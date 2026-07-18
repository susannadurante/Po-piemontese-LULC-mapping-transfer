# Hierarchical L1 -> L2 classification - whole park (Vercelli + Turin sectors).
# Extension principle: extrapolate the MODEL, not the accuracy. The 5 Random Forest
# models trained on Vercelli are applied unchanged to the whole park, using the
# extended predictor stack (output_spettrale_intero_parco/), pixel-aligned to the
# Vercelli grid. On Vercelli pixels the predictors are identical to training
# (|delta|~=0, already certified by the predictor scripts), so the models stay valid.
#
# Differences from the reference-sector classification: paths point to the extended
# folder for composites/predictors/output; DIR_RF stays on the Vercelli models (they
# are extrapolated, not reloaded from elsewhere); the extended topography is already
# decomposed into aspect_sin/aspect_cos (adaptive load: direct if present, else the
# raw aspect decomposed in the loop); a final clip to the study-area perimeter
# produces the deliverable maps. Everything else is the exact mirror of the
# reference classification. Vercelli maps are not touched (output goes to the
# extended folder).
#
# Hierarchical cascade, pixel by pixel:
#   COD_L1 = 1 (water/wetlands)   -> model L2_1 -> 11 / 12 / 13
#   COD_L1 = 2 (river substrate)  -> model L2_2 -> 21 / 22
#   COD_L1 = 3 (herbaceous/shrub) -> direct class 31  (no L2 model)
#   COD_L1 = 4 (tree cover)       -> model L2_4 -> 41 / 42 / 43 / 44 / 45 / 46
#   COD_L1 = 5 (cultivated)       -> model L2_5 -> 51 / 52 / 53 / 54
#   COD_L1 = 6 (impervious)       -> direct class 61
#   COD_L1 = 7 (extraction areas) -> direct class 71
#
# Output (in output_spettrale_intero_parco/08_classification_output/):
#   classificazione_L1.tif / _L2.tif           extended rectangle (INT1U / INT2S)
#   classificazione_L1_clip.tif / _L2_clip.tif clipped to the study area (deliverable)
# Only classification. Post-processing (authoritative 6/7 masks, majority filter)
# and the Area of Applicability are separate downstream steps.

# Packages ----
library(terra)
library(ranger)

# Folders and parameters ----
# Set the project root directory here, or via the POPARK_DATA environment variable
DIR_BASE <- Sys.getenv("POPARK_DATA", unset = "path/to/project/root")
if (!dir.exists(DIR_BASE))
  stop("Invalid DIR_BASE: ", DIR_BASE,
       "
     Set the path at the top of the script, or the POPARK_DATA variable.")

# extended folders (whole park)
DIR_COMP  <- file.path(DIR_BASE, "output_spettrale_intero_parco")                  # 4 extended composites
DIR_STEP3 <- file.path(DIR_BASE, "output_spettrale_intero_parco", "STEP3_output")  # extended indices/topo/GLCM

# models: still the Vercelli ones (NOT retrained)
DIR_RF    <- file.path(DIR_BASE, "output_spettrale", "07_RF_output")

# extended output (dedicated folder: does not touch Vercelli)
DIR_OUT   <- file.path(DIR_BASE, "output_spettrale_intero_parco", "08_classification_output")
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

# perimeter for the final clip (real study area)
PATH_PERIM_CLIP <- file.path(DIR_BASE, "parcopo_areecontigue_dissolto.gpkg")

NUM_CORES <- 6L
options(ranger.num.threads = NUM_CORES)
terraOptions(memfrac = 0.4)   # smaller blocks: more RAM headroom

PATH_L1      <- file.path(DIR_OUT, "classificazione_L1.tif")
PATH_L2      <- file.path(DIR_OUT, "classificazione_L2.tif")
PATH_L1_CLIP <- file.path(DIR_OUT, "classificazione_L1_clip.tif")
PATH_L2_CLIP <- file.path(DIR_OUT, "classificazione_L2_clip.tif")

cat("Hierarchical L1 -> L2 classification - whole park\n")

# Build the predictor stack (mirror of the extraction, Part A) ----
# Same seasons, bands, indices, topographic and GLCM as the extraction, with the
# SAME names. GLCM_COR are excluded: no model uses them (structural NAs).
stagioni <- c("MAM", "JJA", "SON", "DJF")
bande    <- c("B02","B03","B04","B05","B06","B07","B08","B8A","B11","B12")
indici   <- c("NDVI","NDRE","NDWI_McF","NDWI_Gao","MNDWI",
              "NDBI","NBR2","EVI2","BSI","GLI","GVMI")
glcm_metriche <- c("CON","HOM","ENT")   # NO "COR": not used by the models (NA)

nome_composito <- function(s)
  file.path(DIR_COMP,
    paste0("S2_MSIL2A_", s, "2019_T32TLQ-32TLR-32TMQ-32TMR_mediana_10bande.tif"))

# reference raster = EXTENDED MAM composite (defines grid, CRS, resolution). Being
# in the extended folder, r_ref is automatically the extended grid (aligned to
# Vercelli upstream). All other layers are brought onto it.
r_ref <- rast(nome_composito("MAM"))
names(r_ref) <- paste0(bande, "_MAM")
cat(sprintf("Reference grid (extended): %d rows x %d cols, res %.1f m\n",
    nrow(r_ref), ncol(r_ref), res(r_ref)[1]))

layer_list <- list()   # one single-band SpatRaster per predictor
n_resampled <- 0L      # how many layers actually needed resampling

# align to r_ref ONLY if the grid does not ALREADY match (same extent/res/dim/crs)
allinea <- function(r, metodo) {
  if (compareGeom(r, r_ref, stopOnError = FALSE, messages = FALSE)) return(r)
  n_resampled <<- n_resampled + 1L
  resample(r, r_ref, method = metodo)
}

# use the <name>_filled.tif (NAs filled by gap-filling) if present, else the
# original. Makes the stack gap-free before classification.
con_filled <- function(path) {
  f_fill <- sub("\\.tif$", "_filled.tif", path)
  if (file.exists(f_fill)) f_fill else path
}

# Sentinel-2 bands (40 layers) - loaded as DN (x10000); the /10000 division happens
# in the block loop, as at extraction.
for (s in stagioni) {
  f  <- con_filled(nome_composito(s))
  rr <- rast(f); names(rr) <- paste0(bande, "_", s)
  r  <- allinea(rr, "near")
  for (b in bande) { nm <- paste0(b, "_", s); layer_list[[nm]] <- r[[nm]] }
  cat(sprintf("   bands %s loaded\n", s))
}

# spectral indices (44 layers) - already in 0-1, no transform
for (s in stagioni) for (idx in indici) {
  nm <- paste0(idx, "_", s)
  r  <- rast(con_filled(file.path(DIR_STEP3, paste0(nm, ".tif")))); names(r) <- nm
  layer_list[[nm]] <- allinea(r, "near")
}
cat("Spectral indices: 44 layers\n")

# topographic - ADAPTIVE aspect handling.
# In Vercelli aspect was raw (aspect.tif) and decomposed into sin/cos in the loop.
# In the extended folder it is already decomposed into aspect_sin.tif/aspect_cos.tif.
# If the sin/cos files exist -> load them directly (no decomposition); otherwise ->
# load raw aspect.tif and decompose in the loop. Either way the final predictors are
# aspect_sin/aspect_cos, identical to training.
f_asin <- con_filled(file.path(DIR_STEP3, "aspect_sin.tif"))
f_acos <- con_filled(file.path(DIR_STEP3, "aspect_cos.tif"))
ASPECT_PREDECOMP <- file.exists(f_asin) && file.exists(f_acos)

topo <- list(dist_fiume = "dist_fiume.tif", dem = "dem.tif", slope = "slope.tif")
if (ASPECT_PREDECOMP) {
  topo$aspect_sin <- "aspect_sin.tif"
  topo$aspect_cos <- "aspect_cos.tif"
} else {
  topo$aspect_raw <- "aspect.tif"   # fallback: decomposed in the loop
}
for (nm in names(topo)) {
  r <- rast(con_filled(file.path(DIR_STEP3, topo[[nm]]))); names(r) <- nm
  layer_list[[nm]] <- allinea(r, "bilinear")
}
cat(sprintf("Topographic loaded (aspect %s)\n",
    if (ASPECT_PREDECOMP) "already decomposed: sin/cos direct"
    else "raw: aspect_raw -> sin/cos in the loop"))

# GLCM CON/HOM/ENT (12 layers) - use the _v3_filled.tif (NAs filled) if present,
# else the originals _v3.tif.
n_glcm_filled <- 0L
for (s in stagioni) for (m in glcm_metriche) {
  nm     <- paste0("GLCM_", m, "_", s)
  f_fill <- file.path(DIR_STEP3, paste0(nm, "_v3_filled.tif"))
  f_orig <- file.path(DIR_STEP3, paste0(nm, "_v3.tif"))
  f      <- if (file.exists(f_fill)) { n_glcm_filled <- n_glcm_filled + 1L; f_fill } else f_orig
  r <- rast(f); names(r) <- nm
  layer_list[[nm]] <- allinea(r, "near")
}
cat(sprintf("GLCM: 12 layers (%d filled) | layers resampled: %d/%d\n",
    n_glcm_filled, n_resampled, length(layer_list)))

# single predictor stack
all_rast <- terra::rast(layer_list)
cat(sprintf("Predictor stack: %d layers\n", nlyr(all_rast)))

# predictor names AVAILABLE after the loop transforms
BANDE_COLS  <- as.vector(outer(bande, stagioni, function(b, s) paste0(b, "_", s)))
if (ASPECT_PREDECOMP) {
  disponibili <- names(all_rast)                                   # sin/cos already present
} else {
  disponibili <- c(setdiff(names(all_rast), "aspect_raw"), "aspect_sin", "aspect_cos")
}

# Load the models + each model's predictors ----
carica_modello <- function(m) {
  f <- file.path(DIR_RF, m, "rf_model_final.rds")
  if (!file.exists(f)) stop(sprintf("Model not found: %s\n  (the models are the Vercelli ones, 07_RF_output)", f))
  readRDS(f)
}
mod_L1   <- carica_modello("L1")
MODELLI_L2 <- list(
  L2_1 = carica_modello("L2_1"), L2_2 = carica_modello("L2_2"),
  L2_4 = carica_modello("L2_4"), L2_5 = carica_modello("L2_5")
)
preds_L1  <- mod_L1$forest$independent.variable.names
preds_L2  <- lapply(MODELLI_L2, function(mm) mm$forest$independent.variable.names)

# check: every required predictor must be reconstructible from the stack
tutti_pred <- unique(c(preds_L1, unlist(preds_L2)))
mancanti   <- setdiff(tutti_pred, disponibili)
if (length(mancanti) > 0)
  stop("Predictors required by the models but absent from the stack:\n  ",
       paste(mancanti, collapse = ", "))
cat(sprintf("Models loaded. Predictors: L1=%d | L2_1=%d L2_2=%d L2_4=%d L2_5=%d\n",
    length(preds_L1), length(preds_L2$L2_1), length(preds_L2$L2_2),
    length(preds_L2$L2_4), length(preds_L2$L2_5)))

# Helper functions ----
# turn a block's raw values into the predictors expected by the model:
#   bands / 10000 (DN -> reflectance 0-1); (only if raw aspect) aspect_raw -> sin/cos
prepara_blocco <- function(v) {
  v[BANDE_COLS] <- v[BANDE_COLS] / 10000
  if (!ASPECT_PREDECOMP) {
    asp <- v[["aspect_raw"]] * pi / 180
    v[["aspect_sin"]] <- sin(asp)
    v[["aspect_cos"]] <- cos(asp)
    v[["aspect_raw"]] <- NULL
  }
  v
}

# ranger prediction -> integer codes (factor levels "11","12",... -> integers)
predici <- function(model, df) {
  p <- predict(model, data = df, num.threads = NUM_CORES, verbose = FALSE)$predictions
  as.integer(as.character(p))
}

# Block classification (L1, then hierarchical L2) ----
# FIXED-size row blocks, controlled manually (terra's automatic plan ignores the
# extra RAM ranger uses -> it would pick a single block -> out of memory).
# 128 rows (not 256 as in Vercelli): the extended rectangle is ~11,500 px wide, so a
# 256-row block would load ~2.4 GB of predictors alone -> saturation risk on 16 GB.
# At 128 the per-block footprint halves (~1.2 GB).
RIGHE_BLOCCO <- 128L
n_righe   <- nrow(r_ref)
n_blocchi <- ceiling(n_righe / RIGHE_BLOCCO)

out_L1 <- rast(r_ref, nlyrs = 1); names(out_L1) <- "COD_L1"
out_L2 <- rast(r_ref, nlyrs = 1); names(out_L2) <- "COD_L2"

readStart(all_rast)
writeStart(out_L1, PATH_L1, overwrite = TRUE, datatype = "INT1U")
writeStart(out_L2, PATH_L2, overwrite = TRUE, datatype = "INT2S")
cat(sprintf("\nClassifying in %d blocks of %d rows...\n", n_blocchi, RIGHE_BLOCCO))

for (i in seq_len(n_blocchi)) {
  riga0  <- (i - 1L) * RIGHE_BLOCCO + 1L
  nrighe <- min(RIGHE_BLOCCO, n_righe - riga0 + 1L)

  # raw block values: one column per layer, named after the layer
  v <- readValues(all_rast, row = riga0, nrows = nrighe, dataframe = TRUE)
  v <- prepara_blocco(v)
  n <- nrow(v)

  # level 1: macrotype for every pixel with complete predictors
  l1  <- rep(NA_integer_, n)
  ok1 <- stats::complete.cases(v[, preds_L1, drop = FALSE])
  if (any(ok1)) l1[ok1] <- predici(mod_L1, v[ok1, preds_L1, drop = FALSE])

  # level 2: hierarchical refinement
  l2 <- rep(NA_integer_, n)
  # macrotypes without subclasses -> direct L2 class
  l2[which(l1 == 3L)] <- 31L
  l2[which(l1 == 6L)] <- 61L
  l2[which(l1 == 7L)] <- 71L
  # macrotypes with subclasses -> dedicated L2 model, only on that macrotype's pixels
  for (k in names(MODELLI_L2)) {
    kk  <- as.integer(sub("L2_", "", k))                  # "L2_4" -> 4
    sel <- which(l1 == kk & stats::complete.cases(v[, preds_L2[[k]], drop = FALSE]))
    if (length(sel)) l2[sel] <- predici(MODELLI_L2[[k]], v[sel, preds_L2[[k]], drop = FALSE])
  }

  writeValues(out_L1, l1, riga0, nrighe)
  writeValues(out_L2, l2, riga0, nrighe)

  rm(v, l1, l2); if (i %% 5L == 0L) gc(verbose = FALSE)
  cat(sprintf("   block %d/%d (rows %d-%d)\n", i, n_blocchi, riga0, riga0 + nrighe - 1L))
}

writeStop(out_L1); writeStop(out_L2); readStop(all_rast)
cat("Rasters (extended rectangle) written.\n")

# Final clip to the study area (_dissolto) ----
# clip the maps to the real park perimeter: crop (shrinks the extent) + mask (NA
# outside the polygon). The extended-rectangle rasters stay on disk for the Vercelli
# new-vs-old check and for diagnostics.
if (!file.exists(PATH_PERIM_CLIP)) {
  cat(sprintf("\n[WARNING] Clip perimeter not found:\n   %s\n   -> skipping the clip. Fix PATH_PERIM_CLIP and re-run only this section.\n",
      PATH_PERIM_CLIP))
} else {
  cat("\nFinal clip to the study area...\n")
  perim <- vect(PATH_PERIM_CLIP)
  if (!same.crs(perim, r_ref)) perim <- project(perim, crs(r_ref))

  clip_uno <- function(path_in, path_out, dtype) {
    r <- rast(path_in)
    r <- crop(r, perim)
    mask(r, perim, filename = path_out, overwrite = TRUE, datatype = dtype)
  }
  clip_uno(PATH_L1, PATH_L1_CLIP, "INT1U")
  clip_uno(PATH_L2, PATH_L2_CLIP, "INT2S")
  cat("Clipped maps written (_clip).\n")
}

# Summary ----
# pixels per L2 class in the final map (clip if present, else the rectangle)
path_riepilogo <- if (file.exists(PATH_L2_CLIP)) PATH_L2_CLIP else PATH_L2
freq_L2 <- freq(rast(path_riepilogo))
cat("\nClassification summary - whole park\n")
cat(sprintf("   Pixels per L2 class (%s):\n", basename(path_riepilogo)))
for (j in seq_len(nrow(freq_L2)))
  cat(sprintf("     %-4s : %s px\n", freq_L2$value[j],
      format(freq_L2$count[j], big.mark = ".")))
cat(sprintf("\n   Output:\n     %s\n     %s\n     %s\n     %s\n",
    PATH_L1, PATH_L2, PATH_L1_CLIP, PATH_L2_CLIP))
cat("   The Vercelli-sector maps were NOT modified.\n")

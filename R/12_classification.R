# Hierarchical L1 -> L2 classification (applies the trained RF models to the image).
# Applies the 5 trained Random Forest models to the WHOLE image and produces the
# land-cover map, in a two-level cascade:
#   1. the L1 model classifies EVERY pixel into one of the 7 macrotypes (COD_L1);
#   2. where the macrotype has subclasses, the matching L2 model refines the class
#      (COD_L2); where it does not, the L2 class is direct.
# Output: two rasters, the macrotype map (L1) and the final map (L2).
#
# Hierarchical cascade, pixel by pixel:
#   COD_L1 = 1 (water/wetlands)     -> model L2_1 -> 11 / 12 / 13
#   COD_L1 = 2 (river substrate)    -> model L2_2 -> 21 / 22
#   COD_L1 = 3 (herbaceous/shrub)   -> direct class 31  (no L2 model)
#   COD_L1 = 4 (tree cover)         -> model L2_4 -> 41 / 42 / 43 / 44 / 45 / 46
#   COD_L1 = 5 (cultivated)         -> model L2_5 -> 51 / 52 / 53 / 54
#   COD_L1 = 6 (impervious)         -> direct class 61
#   COD_L1 = 7 (extraction areas)   -> direct class 71
#
# Critical requirement: the stack must MATCH the one used at extraction. The models
# were trained on predictors with precise values and names, so here the SAME stack
# is rebuilt over the whole image with the SAME transforms (bands /10000, aspect
# decomposed into sin/cos, same resampling, same layer names). Any value difference
# feeds the model different inputs and silently corrupts the predictions. The stack
# is therefore the exact MIRROR of the extraction (Part A), without the physical
# filter or z-statistics (those cleaned the TRAINING pixels; here everything is
# classified).
#
# Predictor selection is automatic: each model remembers its predictors in
# model$forest$independent.variable.names, so the stack is subset to those names
# per model (works identically for selected-subset and full-set models). GLCM_COR
# are never used (structural NAs) and are not even loaded.
#
# Block processing (16 GB RAM): the image is ~16 M pixels x ~100 layers, too big
# for memory. It is read/predicted/written in fixed-size row BLOCKS, so the RAM per
# block stays bounded. Each pixel is predicted once by L1 and once by its L2 model.
#
# This script does ONLY the classification. Post-processing (majority filter,
# authoritative masks for 6/7) and the Area of Applicability (CAST::aoa) for the
# Turin sector are separate downstream steps.

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
DIR_COMP  <- file.path(DIR_BASE, "output_spettrale")                 # seasonal composites
DIR_STEP3 <- file.path(DIR_BASE, "output_spettrale", "STEP3_output") # indices, topographic, GLCM
DIR_RF    <- file.path(DIR_BASE, "output_spettrale", "07_RF_output") # trained models
DIR_OUT   <- file.path(DIR_BASE, "output_spettrale", "08_classification_output")
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

NUM_CORES <- 6L
options(ranger.num.threads = NUM_CORES)
terraOptions(memfrac = 0.4)   # smaller blocks: more RAM headroom

PATH_L1 <- file.path(DIR_OUT, "classificazione_L1.tif")
PATH_L2 <- file.path(DIR_OUT, "classificazione_L2.tif")

cat("Hierarchical L1 -> L2 classification\n")

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

# reference raster = MAM composite (defines grid, CRS, resolution). All other
# layers are brought onto this grid, as at extraction.
r_ref <- rast(nome_composito("MAM"))
names(r_ref) <- paste0(bande, "_MAM")
cat(sprintf("Reference grid: %d rows x %d cols, res %.1f m\n",
    nrow(r_ref), ncol(r_ref), res(r_ref)[1]))

layer_list <- list()   # one single-band SpatRaster per predictor
n_resampled <- 0L      # how many layers actually needed resampling

# align to r_ref ONLY if the grid does not ALREADY match (same extent/res/dim/crs).
# Predictors are produced on the 10 m composite grid, so they usually match: then
# the resample - very costly over millions of pixels - is skipped entirely and the
# stack is built lazily.
allinea <- function(r, metodo) {
  if (compareGeom(r, r_ref, stopOnError = FALSE, messages = FALSE)) return(r)
  n_resampled <<- n_resampled + 1L
  resample(r, r_ref, method = metodo)
}

# use the <name>_filled.tif (NAs filled by the gap-filling step) if present,
# otherwise the original. Makes the stack gap-free before classification.
con_filled <- function(path) {
  f_fill <- sub("\\.tif$", "_filled.tif", path)
  if (file.exists(f_fill)) f_fill else path
}

# Sentinel-2 bands (40 layers) - loaded as DN (x10000); the /10000 division
# happens in the block loop, as at extraction.
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

# topographic (4 layers): dist_fiume, dem, slope, aspect_raw (bilinear).
# aspect_raw is in degrees; it is decomposed into sin/cos in the block loop.
topo <- list(dist_fiume = "dist_fiume.tif", dem = "dem.tif",
             slope = "slope.tif", aspect_raw = "aspect.tif")
for (nm in names(topo)) {
  r <- rast(con_filled(file.path(DIR_STEP3, topo[[nm]]))); names(r) <- nm
  layer_list[[nm]] <- allinea(r, "bilinear")
}
cat("Topographic: 4 layers (aspect_raw -> sin/cos in the loop)\n")

# GLCM CON/HOM/ENT (12 layers) - use the _v3_filled.tif files (NAs filled) if
# present, otherwise the originals _v3.tif. Texture filling closes the gaps in the
# riparian corridor before classification; the class stays decided by the model on
# real bands/indices.
n_glcm_filled <- 0L
for (s in stagioni) for (m in glcm_metriche) {
  nm     <- paste0("GLCM_", m, "_", s)
  f_fill <- file.path(DIR_STEP3, paste0(nm, "_v3_filled.tif"))
  f_orig <- file.path(DIR_STEP3, paste0(nm, "_v3.tif"))
  f      <- if (file.exists(f_fill)) { n_glcm_filled <- n_glcm_filled + 1L; f_fill } else f_orig
  r <- rast(f); names(r) <- nm
  layer_list[[nm]] <- allinea(r, "near")
}
cat(sprintf("GLCM: 12 layers (%d filled) | total resampling: %d/100\n",
    n_glcm_filled, n_resampled))

# single stack (100 layers: 40 raw bands + 44 indices + 4 topo + 12 GLCM). After
# the loop transforms they become the model's 101 predictors (aspect_raw ->
# aspect_sin + aspect_cos).
all_rast <- terra::rast(layer_list)
cat(sprintf("Predictor stack: %d layers\n", nlyr(all_rast)))

# predictor names AVAILABLE after the loop transforms (aspect_raw -> sin/cos)
BANDE_COLS  <- as.vector(outer(bande, stagioni, function(b, s) paste0(b, "_", s)))
disponibili <- c(setdiff(names(all_rast), "aspect_raw"), "aspect_sin", "aspect_cos")

# Load the models + each model's predictors ----
carica_modello <- function(m) {
  f <- file.path(DIR_RF, m, "rf_model_final.rds")
  if (!file.exists(f)) stop(sprintf("Model not found: %s\n  Train it first for %s.", f, m))
  readRDS(f)
}
mod_L1   <- carica_modello("L1")
MODELLI_L2 <- list(
  L2_1 = carica_modello("L2_1"), L2_2 = carica_modello("L2_2"),
  L2_4 = carica_modello("L2_4"), L2_5 = carica_modello("L2_5")
)
preds_L1  <- mod_L1$forest$independent.variable.names
preds_L2  <- lapply(MODELLI_L2, function(mm) mm$forest$independent.variable.names)

# check: every required predictor must be reconstructible from the stack.
# If a name is missing -> stop immediately (explicit error, not wrong predictions).
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
#   bands / 10000 (DN -> reflectance 0-1); aspect_raw -> aspect_sin, aspect_cos
prepara_blocco <- function(v) {
  v[BANDE_COLS] <- v[BANDE_COLS] / 10000
  asp <- v[["aspect_raw"]] * pi / 180
  v[["aspect_sin"]] <- sin(asp)
  v[["aspect_cos"]] <- cos(asp)
  v[["aspect_raw"]] <- NULL
  v
}

# ranger prediction -> integer codes (factor levels "11","12",... -> integers)
predici <- function(model, df) {
  p <- predict(model, data = df, num.threads = NUM_CORES, verbose = FALSE)$predictions
  as.integer(as.character(p))
}

# Block classification (L1, then hierarchical L2) ----
# FIXED-size row blocks, controlled manually. terra's automatic block plan sizes
# blocks looking ONLY at the rasters, ignoring the extra memory ranger uses to
# predict over millions of pixels: it would pick a single 2541-row block -> out of
# memory. With RIGHE_BLOCCO rows at a time, the RAM per block stays bounded.
RIGHE_BLOCCO <- 256L
n_righe   <- nrow(r_ref)
n_blocchi <- ceiling(n_righe / RIGHE_BLOCCO)

out_L1 <- rast(r_ref, nlyrs = 1); names(out_L1) <- "COD_L1"
out_L2 <- rast(r_ref, nlyrs = 1); names(out_L2) <- "COD_L2"

readStart(all_rast)
# writeStart opens the files for writing; its block plan is IGNORED - the fixed row
# blocks above are used instead. writeValues accepts any start/nrows pair, as long
# as the calls cover all rows in order.
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
cat("Rasters written.\n")

# Summary ----
# pixels per L2 class in the final map (sanity check)
freq_L2 <- freq(rast(PATH_L2))
cat("\nClassification summary - pixels per L2 class:\n")
for (j in seq_len(nrow(freq_L2)))
  cat(sprintf("     %-4s : %s px\n", freq_L2$value[j],
      format(freq_L2$count[j], big.mark = ".")))
cat(sprintf("\n   Output:\n     %s\n     %s\n", PATH_L1, PATH_L2))
cat("   Next steps: post-processing (masks for 6/7, majority filter) + external\n")
cat("   validation + Area of Applicability (CAST::aoa) for the Turin sector.\n")

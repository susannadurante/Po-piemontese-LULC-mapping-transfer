# Export of the predictor stack, for the Area of Applicability analysis.
# The classification script builds the 100-layer stack in memory and discards it
# after predicting: no stack .tif exists on disk. The AoA needs one. This script
# replicates EXACTLY the stack section of the whole-park classification (same
# seasons, bands, indices, topographic layers, GLCM, same _filled files, same
# alignment) and writes to disk ONLY the predictors actually used by the selected
# model. Saving just those, instead of all 100 layers, keeps the file far lighter:
# the AoA does not use the others.
#
# PARAMETRIC: exports the stack of ONE model per run, selected by the MODELLO
#   variable below, and writes it under the name the AoA script expects
#   (stack_predittori_<MODELLO>_intero_parco.tif). Run it once per model whose
#   applicability is of interest, then run the AoA script with the same MODELLO.
#
# Consistency guarantee: the construction logic is identical to the classification.
#   In addition, a final check verifies that every predictor of the model is present
#   in the exported stack, otherwise the AoA would work on inputs different from
#   those used for the classification.
# Run from a terminal (Rscript), never from an interactive IDE.

# Packages ----
library(terra)
library(ranger)

# Model selection ----
MODELLO <- "L2_4"   # the model for this run (L1, L2_1, L2_2, L2_4 or L2_5)

# Folders and parameters ----
# The whole-park seasonal composites and predictors live under
# output_spettrale_intero_parco/ (composites in the root, the rest in STEP3_output/).
# Set the project root directory here, or via the POPARK_DATA environment variable
DIR_BASE <- Sys.getenv("POPARK_DATA", unset = "path/to/project/root")
if (!dir.exists(DIR_BASE))
  stop("Invalid DIR_BASE: ", DIR_BASE,
       "\n     Set the path at the top of the script, or the POPARK_DATA variable.")

DIR_COMP  <- file.path(DIR_BASE, "output_spettrale_intero_parco")
DIR_STEP3 <- file.path(DIR_BASE, "output_spettrale_intero_parco", "STEP3_output")
DIR_RF    <- file.path(DIR_BASE, "output_spettrale", "07_RF_output")   # models stay the ones trained on Vercelli
DIR_OUT   <- file.path(DIR_BASE, "output_spettrale_intero_parco", "08_classification_output")
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

# output: stack of this model's predictors only, over the whole park
PATH_OUT <- file.path(DIR_OUT, sprintf("stack_predittori_%s_intero_parco.tif", MODELLO))

NUM_CORES <- 6L
terraOptions(memfrac = 0.4)

stagioni <- c("MAM", "JJA", "SON", "DJF")
bande    <- c("B02","B03","B04","B05","B06","B07","B08","B8A","B11","B12")
indici   <- c("NDVI","NDRE","NDWI_McF","NDWI_Gao","MNDWI",
              "NDBI","NBR2","EVI2","BSI","GLI","GVMI")
glcm_metriche <- c("CON","HOM","ENT")

nome_composito <- function(s)
  file.path(DIR_COMP,
    paste0("S2_MSIL2A_", s, "2019_T32TLQ-32TLR-32TMQ-32TMR_mediana_10bande.tif"))

cat(sprintf("Export of the predictor stack - model %s\n\n", MODELLO))

# Predictors of the model ----
f_mod <- file.path(DIR_RF, MODELLO, "rf_model_final.rds")
if (!file.exists(f_mod)) stop("Model not found: ", f_mod)
mod   <- readRDS(f_mod)
preds <- mod$forest$independent.variable.names
cat(sprintf("   Predictors of model %s: %d\n", MODELLO, length(preds)))

# aspect: the model uses aspect_sin/aspect_cos; the raw stack holds aspect_raw
usa_aspect <- any(c("aspect_sin", "aspect_cos") %in% preds)

# Reference grid ----
r_ref <- rast(nome_composito("MAM"))
names(r_ref) <- paste0(bande, "_MAM")
cat(sprintf("   Grid: %d x %d, res %.1f m\n", nrow(r_ref), ncol(r_ref), res(r_ref)[1]))

layer_list  <- list()
n_resampled <- 0L

# align to r_ref only if the grid does not already match
allinea <- function(r, metodo) {
  if (compareGeom(r, r_ref, stopOnError = FALSE, messages = FALSE)) return(r)
  n_resampled <<- n_resampled + 1L
  resample(r, r_ref, method = metodo)
}

# use the <name>_filled.tif (NAs filled by gap-filling) if present, else the original
con_filled <- function(path) {
  f_fill <- sub("\\.tif$", "_filled.tif", path)
  if (file.exists(f_fill)) f_fill else path
}

# keep only the layers this model needs: filter by name
serve <- function(nm) nm %in% preds

# Sentinel-2 bands
for (s in stagioni) {
  nomi_b <- paste0(bande, "_", s)
  if (!any(serve(nomi_b))) next
  rr <- rast(con_filled(nome_composito(s))); names(rr) <- nomi_b
  r  <- allinea(rr, "near")
  for (b in bande) { nm <- paste0(b, "_", s); if (serve(nm)) layer_list[[nm]] <- r[[nm]] }
}

# spectral indices
for (s in stagioni) for (idx in indici) {
  nm <- paste0(idx, "_", s); if (!serve(nm)) next
  r <- rast(con_filled(file.path(DIR_STEP3, paste0(nm, ".tif")))); names(r) <- nm
  layer_list[[nm]] <- allinea(r, "near")
}

# topographic (dist_fiume, dem, slope; aspect_raw only if sin/cos are needed)
topo <- list(dist_fiume = "dist_fiume.tif", dem = "dem.tif", slope = "slope.tif")
for (nm in names(topo)) {
  if (!serve(nm)) next
  r <- rast(con_filled(file.path(DIR_STEP3, topo[[nm]]))); names(r) <- nm
  layer_list[[nm]] <- allinea(r, "bilinear")
}
if (usa_aspect) {
  r  <- rast(con_filled(file.path(DIR_STEP3, "aspect.tif"))); names(r) <- "aspect_raw"
  ar <- allinea(r, "bilinear")
  asp <- ar * pi / 180
  if ("aspect_sin" %in% preds) { s1 <- sin(asp); names(s1) <- "aspect_sin"; layer_list[["aspect_sin"]] <- s1 }
  if ("aspect_cos" %in% preds) { c1 <- cos(asp); names(c1) <- "aspect_cos"; layer_list[["aspect_cos"]] <- c1 }
}

# GLCM
for (s in stagioni) for (m in glcm_metriche) {
  nm <- paste0("GLCM_", m, "_", s); if (!serve(nm)) next
  f_fill <- file.path(DIR_STEP3, paste0(nm, "_v3_filled.tif"))
  f_orig <- file.path(DIR_STEP3, paste0(nm, "_v3.tif"))
  f <- if (file.exists(f_fill)) f_fill else f_orig
  r <- rast(f); names(r) <- nm
  layer_list[[nm]] <- allinea(r, "near")
}

# Stack and check ----
# raw bands / 10000, as the classification does inside its block loop, so the stack
# values coincide with those the model was trained on.
BANDE_COLS <- as.vector(outer(bande, stagioni, function(b, s) paste0(b, "_", s)))
for (nm in intersect(BANDE_COLS, names(layer_list)))
  layer_list[[nm]] <- layer_list[[nm]] / 10000

stk <- rast(layer_list)
cat(sprintf("\n   Layers in the stack: %d | resampled: %d\n", nlyr(stk), n_resampled))

manca <- setdiff(preds, names(stk))
if (length(manca) > 0)
  stop("Model predictors ABSENT from the stack: ", paste(manca, collapse = ", "))
stk <- stk[[preds]]   # same order as the model predictors
cat(sprintf("   All predictors of %s present and ordered.\n", MODELLO))

# Write ----
cat(sprintf("\n   Writing %s ...\n", basename(PATH_OUT)))
writeRaster(stk, PATH_OUT, overwrite = TRUE,
            gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES"))
cat("   Done.\n")
cat(sprintf("\n   Now run the Area of Applicability script with MODELLO = \"%s\".\n", MODELLO))

# Spectral indices - whole park (Vercelli + Turin sectors).
# Computes the 11 seasonal spectral indices on the whole-park composites and writes
# them to output_spettrale_intero_parco/STEP3_output/, with the same canonical names
# (<INDEX>_<season>.tif) used for Vercelli.
# Formulas identical to the reference-sector scripts:
#   3 base (normalised difference, raw DN): NDVI, NDWI_McF, NDWI_Gao
#   8 others: MNDWI, BSI, EVI2, NDRE, NDBI, NBR2, GVMI, GLI
#   EVI2 and GVMI use bands in REFLECTANCE 0-1 (/10000); the others use DN.
# Band order in the composite:
#   B02=1 B03=2 B04=3 B05=4 B06=5 B07=6 B08=7 B8A=8 B11=9 B12=10

# Packages ----
library(terra)

# Folders ----
# Set the project root directory here, or via the POPARK_DATA environment variable
DIR_BASE <- Sys.getenv("POPARK_DATA", unset = "path/to/project/root")
if (!dir.exists(DIR_BASE))
  stop("Invalid DIR_BASE: ", DIR_BASE,
       "
     Set the path at the top of the script, or the POPARK_DATA variable.")
DIR_VERCELLI  <- file.path(DIR_BASE, "output_spettrale")                 # for the validation
DIR_STEP3_VC  <- file.path(DIR_VERCELLI, "STEP3_output")
DIR_OUTPUT    <- file.path(DIR_BASE, "output_spettrale_intero_parco")    # whole park
DIR_STEP3     <- file.path(DIR_OUTPUT, "STEP3_output")
dir.create(DIR_STEP3, showWarnings = FALSE, recursive = TRUE)

COMPOSITI <- list(
  MAM = file.path(DIR_OUTPUT, "S2_MSIL2A_MAM2019_T32TLQ-32TLR-32TMQ-32TMR_mediana_10bande.tif"),
  JJA = file.path(DIR_OUTPUT, "S2_MSIL2A_JJA2019_T32TLQ-32TLR-32TMQ-32TMR_mediana_10bande.tif"),
  SON = file.path(DIR_OUTPUT, "S2_MSIL2A_SON2019_T32TLQ-32TLR-32TMQ-32TMR_mediana_10bande.tif"),
  DJF = file.path(DIR_OUTPUT, "S2_MSIL2A_DJF2019_T32TLQ-32TLR-32TMQ-32TMR_mediana_10bande.tif")
)
mancanti <- unlist(COMPOSITI)[!file.exists(unlist(COMPOSITI))]
if (length(mancanti) > 0) stop("Missing composites (generate the 4 extended composites first):\n",
                               paste(" -", mancanti, collapse = "\n"))

# Write helper ----
salva <- function(r, nome, stagione) {
  names(r) <- paste0(nome, "_", stagione)
  writeRaster(r, file.path(DIR_STEP3, paste0(nome, "_", stagione, ".tif")),
              datatype = "FLT4S", overwrite = TRUE,
              gdal = c("COMPRESS=LZW", "TILED=YES"))
  mm <- global(r, c("min", "max"), na.rm = TRUE)
  cat(sprintf("   %-9s min=%.3f  max=%.3f  OK\n", paste0(nome,"_",stagione),
              mm[1,"min"], mm[1,"max"]))
}

cat("11 spectral indices - whole park\n")
cat("Base:  NDVI, NDWI_McF, NDWI_Gao   |   Others: MNDWI, BSI, EVI2, NDRE, NDBI, NBR2, GVMI, GLI\n")
cat("Output:", DIR_STEP3, "\n\n")

# Loop over the 4 seasons ----
for (stagione in names(COMPOSITI)) {
  cat(sprintf("-- Season %s\n", stagione))
  comp <- rast(COMPOSITI[[stagione]])

  B02 <- comp[[1]]; B03 <- comp[[2]]; B04 <- comp[[3]]; B05 <- comp[[4]]
  B08 <- comp[[7]]; B8A <- comp[[8]]; B11 <- comp[[9]]; B12 <- comp[[10]]
  # reflectance 0-1 for indices with absolute constants (EVI2, GVMI)
  B04r <- B04 / 10000;  B08r <- B08 / 10000;  B12r <- B12 / 10000

  # 3 base indices (raw DN, normalised ratios)
  salva((B08 - B04) / (B08 + B04),               "NDVI",     stagione)
  salva((B03 - B08) / (B03 + B08),               "NDWI_McF", stagione)   # McFeeters 1996
  salva((B08 - B11) / (B08 + B11),               "NDWI_Gao", stagione)   # Gao 1996

  # 8 further indices
  salva((B03 - B11) / (B03 + B11),                                   "MNDWI", stagione)
  salva(((B11 + B04) - (B08 + B02)) / ((B11 + B04) + (B08 + B02)),   "BSI",   stagione)
  salva(2.5 * (B08r - B04r) / (B08r + 2.4 * B04r + 1),               "EVI2",  stagione)   # reflectance 0-1
  salva((B8A - B05) / (B8A + B05),                                   "NDRE",  stagione)
  salva((B11 - B08) / (B11 + B08),                                   "NDBI",  stagione)
  salva((B11 - B12) / (B11 + B12),                                   "NBR2",  stagione)
  salva(((B08r + 0.1) - (B12r + 0.02)) / ((B08r + 0.1) + (B12r + 0.02)), "GVMI", stagione) # reflectance 0-1
  salva((2*B03 - B04 - B02) / (2*B03 + B04 + B02),                   "GLI",   stagione)   # Louhaichi 2001

  rm(comp, B02, B03, B04, B05, B08, B8A, B11, B12, B04r, B08r, B12r); gc()
  cat("\n")
}

# Check produced files ----
indici   <- c("NDVI","NDWI_McF","NDWI_Gao","MNDWI","BSI","EVI2","NDRE","NDBI","NBR2","GVMI","GLI")
stagioni <- c("MAM","JJA","SON","DJF")
n_ok <- 0
for (ind in indici) for (s in stagioni) {
  if (file.exists(file.path(DIR_STEP3, paste0(ind,"_",s,".tif")))) n_ok <- n_ok + 1
  else cat(sprintf("   MISSING: %s_%s.tif\n", ind, s))
}
cat(sprintf("   %d/%d rasters produced\n\n", n_ok, length(indici)*length(stagioni)))

# Validation on Vercelli ----
# compare the whole-park indices with the EXISTING Vercelli ones on the overlapping
# window. |delta|~=0 confirms formulas and grid alignment (Vercelli pixels are
# identical, the models stay valid).
cat("Validation on Vercelli (expect |delta|max ~= 0):\n")
for (ind in c("NDVI","EVI2","GVMI")) {       # 1 base + 2 with absolute scaling
  f_new <- file.path(DIR_STEP3,    paste0(ind, "_MAM.tif"))
  f_old <- file.path(DIR_STEP3_VC, paste0(ind, "_MAM.tif"))
  if (file.exists(f_new) && file.exists(f_old)) {
    rn <- rast(f_new); ro <- rast(f_old)
    rc <- crop(rn, ext(ro))
    if (!compareGeom(rc, ro, stopOnError = FALSE)) rc <- resample(rc, ro, method = "near")
    dmax <- as.numeric(global(abs(rc - ro), "max", na.rm = TRUE))
    cat(sprintf("   %s_MAM:  |delta|max = %.6g  %s\n", ind, dmax,
                if (is.finite(dmax) && dmax < 1e-4) "OK (identical)" else "!! check"))
  } else {
    cat(sprintf("   %s_MAM: comparison raster not found, skipping.\n", ind))
  }
}

# GLCM texture with cloud mask - whole park (Vercelli + Turin sectors).
# Recomputes the GLCM textures on the whole-park composites and writes them to
# output_spettrale_intero_parco/STEP3_output/, with the same canonical names
# (GLCM_<code>_<season>_v3.tif) used for Vercelli.
# Identical to the reference-sector step (cloud mask B02>4000, normalisation 0-31,
# 5x5 window, 4 directions), with ONE difference: only CON, HOM, ENT are computed
# (not COR). GLCM_COR has STRUCTURAL NAs on uniform surfaces (sigma=0 in the window
# -> undefined correlation); it is already excluded from the z-statistics and from
# the RF predictors (ranger does not handle NAs). As it enters neither the model
# nor the classification, computing it over the whole park would be pure waste (GLCM
# is the heaviest step after the composites).
# Cloud mask: B02 (blue) > 4000 DN = certain cloud -> B08 set to NA before glcm();
# glcm(na_val=NA) excludes those pixels from the 5x5 window.

# Packages ----
library(terra)
library(raster)
library(glcm)

# Folders ----
# Set the project root directory here, or via the POPARK_DATA environment variable
DIR_BASE <- Sys.getenv("POPARK_DATA", unset = "path/to/project/root")
if (!dir.exists(DIR_BASE))
  stop("Invalid DIR_BASE: ", DIR_BASE,
       "
     Set the path at the top of the script, or the POPARK_DATA variable.")
DIR_OUTPUT <- file.path(DIR_BASE, "output_spettrale_intero_parco")     # whole park
DIR_STEP3  <- file.path(DIR_OUTPUT, "STEP3_output")
dir.create(DIR_STEP3, showWarnings = FALSE, recursive = TRUE)

# extended S2 composites (10 bands: B02,B03,B04,B05,B06,B07,B08,B8A,B11,B12)
COMPOSITI <- list(
  MAM = file.path(DIR_OUTPUT, "S2_MSIL2A_MAM2019_T32TLQ-32TLR-32TMQ-32TMR_mediana_10bande.tif"),
  JJA = file.path(DIR_OUTPUT, "S2_MSIL2A_JJA2019_T32TLQ-32TLR-32TMQ-32TMR_mediana_10bande.tif"),
  SON = file.path(DIR_OUTPUT, "S2_MSIL2A_SON2019_T32TLQ-32TLR-32TMQ-32TMR_mediana_10bande.tif"),
  DJF = file.path(DIR_OUTPUT, "S2_MSIL2A_DJF2019_T32TLQ-32TLR-32TMQ-32TMR_mediana_10bande.tif")
)

# Parameters ----
BANDA_B02 <- 1   # blue - cloud mask
BANDA_B08 <- 7   # NIR - GLCM texture
SOGLIA_B02_DN <- 4000  # DN: B02 > 4000 = certain cloud (reflectance > 0.40)

# only the three metrics used by the model (COR excluded: structural NAs)
STATISTICHE <- c("contrast", "homogeneity", "entropy")
SIGLA_MAP   <- c("contrast"="CON", "homogeneity"="HOM", "entropy"="ENT")

cat("GLCM with cloud mask - whole park\n")
cat("Cloud threshold: B02 >", SOGLIA_B02_DN, "DN (reflectance > 0.40)\n")
cat("Metrics: CON, HOM, ENT  (COR excluded: structural NAs, not used by the RF)\n")
cat("Output:", DIR_STEP3, "\n\n")

mancanti <- unlist(COMPOSITI)[!file.exists(unlist(COMPOSITI))]
if (length(mancanti) > 0) stop("Missing composites (generate the 4 extended composites first):\n",
                               paste(" -", mancanti, collapse = "\n"))
cat("Input composites verified.\n\n")

# Loop over composites ----
for (stagione in names(COMPOSITI)) {
  cat("Season", stagione, "\n")
  path_comp <- COMPOSITI[[stagione]]

  # cloud mask (count, with terra)
  comp_terra <- rast(path_comp)
  b02_terra  <- comp_terra[[BANDA_B02]]
  maschera_nuvole <- b02_terra > SOGLIA_B02_DN
  n_nuvole <- global(maschera_nuvole, fun = "sum", na.rm = TRUE)[1,1]
  n_totale <- ncell(b02_terra) - global(is.na(b02_terra), fun = "sum", na.rm = TRUE)[1,1]
  cat(sprintf("   Cloudy pixels (B02 > %d DN): %s / %s (%.3f%%)\n",
              SOGLIA_B02_DN, format(as.integer(n_nuvole), big.mark = "."),
              format(as.integer(n_totale), big.mark = "."),
              round(n_nuvole / n_totale * 100, 3)))

  # B08 masked and normalised 0-31 (RasterLayer for glcm)
  b08_rl <- raster::raster(path_comp, band = BANDA_B08); names(b08_rl) <- "B08"
  b02_rl <- raster::raster(path_comp, band = BANDA_B02); names(b02_rl) <- "B02"

  b08_masked_rl <- calc(stack(b08_rl, b02_rl), fun = function(x) {
    b08 <- x[, 1]; b02 <- x[, 2]
    b08[!is.na(b02) & b02 > SOGLIA_B02_DN] <- NA     # clouds -> NA
    b08 <- round((b08 / 10000) * 31)                 # normalise 0-31
    b08[!is.na(b08) & b08 < 0]  <- 0L
    b08[!is.na(b08) & b08 > 31] <- 31L
    as.integer(b08)
  })
  names(b08_masked_rl) <- "B08"
  n_na_rl <- cellStats(is.na(b08_masked_rl), sum)
  cat(sprintf("   NA pixels in the final RasterLayer: %s\n",
              format(as.integer(n_na_rl), big.mark = ".")))

  # GLCM (5x5 window, 4 directions) - CON, HOM, ENT
  cat("   Computing GLCM (5x5 window, 4 directions)...\n")
  glcm_out <- glcm::glcm(
    b08_masked_rl,
    window     = c(5, 5),
    shift      = list(c(0,1), c(1,1), c(1,0), c(1,-1)),
    statistics = STATISTICHE,        # only CON, HOM, ENT
    na_val     = NA,
    min_x      = 0,
    max_x      = 31
  )

  # save the GLCM rasters
  for (k in 1:nlayers(glcm_out)) {
    stat_key <- gsub("glcm_", "", names(glcm_out)[k])
    sigla    <- SIGLA_MAP[stat_key]
    layer_sr <- rast(glcm_out[[k]])
    nome_out <- paste0("GLCM_", sigla, "_", stagione, "_v3")
    names(layer_sr) <- nome_out
    writeRaster(layer_sr, file.path(DIR_STEP3, paste0(nome_out, ".tif")),
                datatype = "FLT4S", overwrite = TRUE,
                gdal = c("COMPRESS=LZW", "TILED=YES"))
  }

  # check CON range
  con_sr <- rast(file.path(DIR_STEP3, paste0("GLCM_CON_", stagione, "_v3.tif")))
  mm <- global(con_sr, c("min","max"), na.rm = TRUE)
  cat(sprintf("   GLCM_CON_%s_v3: min=%.3f max=%.3f\n", stagione, mm[1,"min"], mm[1,"max"]))
  cat(sprintf("   GLCM_CON/HOM/ENT_%s_v3.tif saved\n\n", stagione))

  rm(comp_terra, b02_terra, maschera_nuvole, b08_rl, b02_rl, b08_masked_rl, glcm_out, con_sr)
  gc()
}

# Check produced files ----
n_ok <- 0
for (sig in c("CON","HOM","ENT")) for (s in names(COMPOSITI)) {
  if (file.exists(file.path(DIR_STEP3, paste0("GLCM_", sig, "_", s, "_v3.tif")))) n_ok <- n_ok + 1
  else cat(sprintf("   MISSING: GLCM_%s_%s_v3.tif\n", sig, s))
}
cat(sprintf("   %d/12 GLCM rasters produced\n", n_ok))

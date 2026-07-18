# Gap-filling of the predictors - whole park (Vercelli + Turin sectors).
# Spatially interpolates the NA pixels of the extended predictors (composites,
# indices, topography), writing <name>_filled.tif versions. Originals stay intact;
# the classification step uses the _filled versions automatically where present.
#
# What it does and does not do (methodological honesty, unchanged from the
# reference-sector version): it is spatial INTERPOLATION (iterated 3x3 focal mean on
# NA cells only), not data recovery. NAs arise in the compositing (median -> NA
# where no valid observation exists in the season); not even a BAP would recover
# them. For single-year data, spatial interpolation is the correct choice. Nothing
# is retrained: the models stay unchanged.
#
# Differences from the reference-sector version:
#   - paths -> output_spettrale_intero_parco/
#   - MAM is NOT skipped (on the extended extent it has 48 border NAs); processed if it has NAs
#   - aspect_sin / aspect_cos added (produced by the extended topography step)
#   - GLCM skipped here (filled by the GLCM-fill step)
#   - OPTIONAL park-perimeter mask: avoids interpolating large out-of-park NA areas
#     (useful for DJF, ~68k NAs). See MASCHERA_PERIMETRO below.

# Packages ----
library(terra)

# Folders ----
# Set the project root directory here, or via the POPARK_DATA environment variable
DIR_BASE <- Sys.getenv("POPARK_DATA", unset = "path/to/project/root")
if (!dir.exists(DIR_BASE))
  stop("Invalid DIR_BASE: ", DIR_BASE,
       "
     Set the path at the top of the script, or the POPARK_DATA variable.")
DIR_COMP  <- file.path(DIR_BASE, "output_spettrale_intero_parco")
DIR_STEP3 <- file.path(DIR_COMP, "STEP3_output")

# Option: limit filling to the park perimeter (+ margin) ----
# TRUE = set out-of-perimeter NAs aside before interpolating -> no iterations wasted
#        on areas that will be masked anyway after the classification.
# TRUE is recommended given the large out-of-park areas in the rectangle corners.
MASCHERA_PERIMETRO <- FALSE
PERIMETRO   <- file.path(DIR_BASE, "parcopo_areecontigue_buffer.gpkg")
MARGINE_M   <- 300

# Parameters ----
stagioni <- c("MAM", "JJA", "SON", "DJF")
indici   <- c("NDVI","NDRE","NDWI_McF","NDWI_Gao","MNDWI","NDBI","NBR2","EVI2","BSI","GLI","GVMI")
W       <- 3L
MAXITER <- 300L

nome_composito <- function(s)
  file.path(DIR_COMP, paste0("S2_MSIL2A_", s, "2019_T32TLQ-32TLR-32TMQ-32TMR_mediana_10bande.tif"))

# Functions ----
na_tot <- function(r) sum(as.numeric(global(r, "isNA")[, 1]))

# park mask (in the CRS/grid of the given raster), as a projected SpatVector
carica_perimetro <- function() {
  if (!MASCHERA_PERIMETRO) return(NULL)
  if (!file.exists(PERIMETRO)) { cat("   ! perimeter not found, mask disabled\n"); return(NULL) }
  vect(PERIMETRO)
}
PERIM <- carica_perimetro()

# where OUTSIDE the park, NAs stay NA and are not counted/filled, because the mask
# is applied AFTER filling (out-of-park pixels revert to NA)
maschera_parco <- function(r) {
  if (is.null(PERIM)) return(r)
  p <- project(PERIM, crs(r))
  if (MARGINE_M > 0) p <- buffer(p, MARGINE_M)
  mask(r, p)
}

# iterative spatial fill by neighbour mean (NA cells only, per band)
riempi <- function(r) {
  iter <- 0L
  repeat {
    if (na_tot(r) == 0 || iter >= MAXITER) break
    iter <- iter + 1L
    r <- focal(r, w = W, fun = "mean", na.policy = "only", na.rm = TRUE)
  }
  list(r = r, iter = iter, na = na_tot(r))
}

# fill a file ONLY if it has NAs; write <name>_filled.tif. Returns a log line.
processa <- function(f_in, etichetta) {
  if (!file.exists(f_in)) return(sprintf("   %-16s ! file not found", etichetta))
  r   <- rast(f_in)
  # if active, keep only NAs inside the park: crop FIRST, so out-of-park NAs leave
  # the count and are not interpolated
  if (!is.null(PERIM)) r <- maschera_parco(r)
  na0 <- na_tot(r)
  if (na0 == 0) return(sprintf("   %-16s NA=0 -> nothing to do", etichetta))
  out   <- riempi(r)
  if (!is.null(PERIM)) out$r <- maschera_parco(out$r)   # clean any out-of-park interpolated borders
  f_out <- sub("\\.tif$", "_filled.tif", f_in)
  writeRaster(out$r, f_out, overwrite = TRUE)
  sprintf("   %-16s NA %d -> %d  (%d steps) -> %s", etichetta, na0, out$na, out$iter, basename(f_out))
}

# Run ----
cat("Gap-filling predictors (whole park)\n")
cat("Perimeter mask:", if (is.null(PERIM)) "NO (whole rectangle)" else "YES (park + margin)", "\n\n")

cat("Seasonal composites (bands):\n")
for (s in stagioni) cat(processa(nome_composito(s), paste0("composito_", s)), "\n")

cat("\nSpectral indices:\n")
for (s in stagioni) for (idx in indici)
  cat(processa(file.path(DIR_STEP3, paste0(idx, "_", s, ".tif")), paste0(idx, "_", s)), "\n")

cat("\nTopographic:\n")
cat(processa(file.path(DIR_STEP3, "slope.tif"),      "slope"),      "\n")
cat(processa(file.path(DIR_STEP3, "dem.tif"),        "dem"),        "\n")
cat(processa(file.path(DIR_STEP3, "aspect_sin.tif"), "aspect_sin"), "\n")
cat(processa(file.path(DIR_STEP3, "aspect_cos.tif"), "aspect_cos"), "\n")

cat("\nDone. The classification step uses the _filled versions where present.\n")
cat("   GLCM: filled separately by the GLCM-fill step.\n")

# Spatial gap-filling of the predictors.
# What it does and does not do (for methodological honesty): NAs arise in the
# compositing (median -> NA where no valid observation exists in the season). They
# are not recovered from the data here: they are INTERPOLATED spatially. This is
# interpolation, not recovery - even a BAP would not recover those pixels (no valid
# observation, nothing to select). For single-year data, spatial interpolation is
# the correct choice. It does not require re-running the pipeline: the trained
# models stay unchanged; only the classification (which reads the _filled rasters)
# is re-run.
#
# Method: spatial fill by neighbour mean (iterated 3x3 focal), on NA pixels only,
# until none remain. Suited to small, sparse gaps (the case here).
# Fills: seasonal composites (bands), indices, slope. MAM is clean -> skipped. GLCM
# are filled by a separate step and skipped here.
# Output: for each raster with NAs, a <name>_filled.tif (originals kept). The
# classification step uses the _filled versions automatically.

# Packages ----
library(terra)

# Folders ----
# Set the project root directory here, or via the POPARK_DATA environment variable
DIR_BASE <- Sys.getenv("POPARK_DATA", unset = "path/to/project/root")
if (!dir.exists(DIR_BASE))
  stop("Invalid DIR_BASE: ", DIR_BASE,
       "
     Set the path at the top of the script, or the POPARK_DATA variable.")
DIR_COMP  <- file.path(DIR_BASE, "output_spettrale")
DIR_STEP3 <- file.path(DIR_BASE, "output_spettrale", "STEP3_output")

# Parameters ----
stagioni <- c("MAM", "JJA", "SON", "DJF")
indici   <- c("NDVI","NDRE","NDWI_McF","NDWI_Gao","MNDWI","NDBI","NBR2","EVI2","BSI","GLI","GVMI")
W       <- 3L
MAXITER <- 300L

nome_composito <- function(s)
  file.path(DIR_COMP, paste0("S2_MSIL2A_", s, "2019_T32TLQ-32TLR-32TMQ-32TMR_mediana_10bande.tif"))

# Functions ----
# total NAs of a raster (also multiband)
na_tot <- function(r) sum(as.numeric(global(r, "isNA")[, 1]))

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
  na0 <- na_tot(r)
  if (na0 == 0) return(sprintf("   %-16s NA=0 -> nothing to do", etichetta))
  out   <- riempi(r)
  f_out <- sub("\\.tif$", "_filled.tif", f_in)
  writeRaster(out$r, f_out, overwrite = TRUE)
  sprintf("   %-16s NA %d -> %d  (%d steps) -> %s", etichetta, na0, out$na, out$iter, basename(f_out))
}

# Run ----
cat("Gap-filling predictors\n")

cat("Seasonal composites (bands):\n")
for (s in stagioni) cat(processa(nome_composito(s), paste0("composito_", s)), "\n")

cat("\nSpectral indices:\n")
for (s in stagioni) for (idx in indici)
  cat(processa(file.path(DIR_STEP3, paste0(idx, "_", s, ".tif")), paste0(idx, "_", s)), "\n")

cat("\nTopographic:\n")
cat(processa(file.path(DIR_STEP3, "slope.tif"), "slope"), "\n")
cat(processa(file.path(DIR_STEP3, "dem.tif"),   "dem"),   "\n")

cat("\nDone. The classification step uses the _filled rasters where present.\n")

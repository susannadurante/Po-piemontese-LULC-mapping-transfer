# GLCM filling - whole park (Vercelli + Turin sectors).
# Closes the NAs in the GLCM textures BEFORE classification, over the whole park.
# Identical to the reference-sector version (iterated 3x3 focal on NA cells only),
# repointed to the extended folder. Consistent with the predictor gap-filling: it
# fills the whole rectangle (the out-of-park area is removed by the final map clip).
#
# Why: texture has NAs where the cloud/fog mask is wide (up to ~0.7%, more than the
# bands). There the classification would leave pixels unclassified. Filling the GLCM
# lets those pixels be classified - and the CLASS is still decided by the RF on the
# REAL bands/indices present: this fills a SECONDARY predictor, not the final map.
# Method: texture is a CONTINUOUS variable -> neighbour mean (3x3 focal), iterated
# until no NA remains (MAXITER cap).
# Input: output_spettrale_intero_parco/STEP3_output/GLCM_<m>_<s>_v3.tif (12 files).
# Output: .../GLCM_<m>_<s>_v3_filled.tif (originals untouched).

# Packages ----
library(terra)

# Folders ----
# Set the project root directory here, or via the POPARK_DATA environment variable
DIR_BASE <- Sys.getenv("POPARK_DATA", unset = "path/to/project/root")
if (!dir.exists(DIR_BASE))
  stop("Invalid DIR_BASE: ", DIR_BASE,
       "
     Set the path at the top of the script, or the POPARK_DATA variable.")
DIR_STEP3 <- file.path(DIR_BASE, "output_spettrale_intero_parco", "STEP3_output")   # extended folder

# Parameters ----
stagioni      <- c("MAM", "JJA", "SON", "DJF")
glcm_metriche <- c("CON", "HOM", "ENT")     # metrics used by the models (no COR)
W       <- 3L
MAXITER <- 300L

# Function ----
riempi_continuo <- function(r) {
  iter <- 0L
  repeat {
    na <- as.numeric(global(r, "isNA")[1, 1])
    if (na == 0 || iter >= MAXITER) break
    iter <- iter + 1L
    r <- focal(r, w = W, fun = "mean", na.policy = "only", na.rm = TRUE)
  }
  list(r = r, iter = iter, resid = as.numeric(global(r, "isNA")[1, 1]))
}

# Run ----
cat("GLCM filling (whole park)\n")
cat("Folder:", DIR_STEP3, "\n\n")

n_ok <- 0
for (s in stagioni) for (m in glcm_metriche) {
  nm    <- paste0("GLCM_", m, "_", s)
  f_in  <- file.path(DIR_STEP3, paste0(nm, "_v3.tif"))
  f_out <- file.path(DIR_STEP3, paste0(nm, "_v3_filled.tif"))

  if (!file.exists(f_in)) { cat(sprintf("   %-14s ! not found (run the GLCM step first)\n", nm)); next }

  r    <- rast(f_in)
  ncel <- ncell(r)
  na0  <- as.numeric(global(r, "isNA")[1, 1])

  out <- riempi_continuo(r)
  writeRaster(out$r, f_out, overwrite = TRUE)
  n_ok <- n_ok + 1

  cat(sprintf("   %-14s NA: %.3f%% -> %.3f%%  (%d steps)  -> %s\n",
      nm, 100 * na0 / ncel, 100 * out$resid / ncel, out$iter, basename(f_out)))
  rm(r, out); gc(verbose = FALSE)
}

cat(sprintf("\nDone. %d/12 GLCM filled. The classification step uses the _v3_filled.tif if present.\n", n_ok))

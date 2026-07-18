# Topography + river distance - whole park (Vercelli + Turin sectors).
# On the extended grid (the whole-park composites) produces:
#   dist_fiume_rete.tif - distance from the Po axis (official hydrographic network,
#                         ARPA/BDTRE, asse_Po.gpkg)
#   dist_fiume_lcp.tif  - distance from the Po channel isolated from LCP
#                         (clc_3liv==511, connected components, drops minor channels)
#   Vercelli check      - compares both with the old dist_fiume.tif (median diff,
#                         95th pct, % pixels |delta|>10 m + delta raster)
#   dem.tif, slope.tif, aspect_sin.tif, aspect_cos.tif - from DTM_merge.tif
#
# Principle: extrapolate the MODEL, not the accuracy. The Vercelli check tells
# whether the new dist_fiume matches (within a few metres) the one the models were
# trained on -> if so, the models are kept without retraining.
# Prerequisite: the extended composites must already exist (they define the grid).

# Packages ----
library(terra)
library(sf)

# namespace conflicts (as in the rest of the pipeline)
extract <- terra::extract

# Folders ----
# Set the project root directory here, or via the POPARK_DATA environment variable
DIR_BASE <- Sys.getenv("POPARK_DATA", unset = "path/to/project/root")
if (!dir.exists(DIR_BASE))
  stop("Invalid DIR_BASE: ", DIR_BASE,
       "
     Set the path at the top of the script, or the POPARK_DATA variable.")
DIR_EXT  <- file.path(DIR_BASE, "output_spettrale_intero_parco")          # extended folder
DIR_OUT  <- file.path(DIR_EXT, "STEP3_output")
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

PATH_ASSE_PO <- file.path(DIR_BASE, "asse_Po.gpkg")                       # hydrographic network (Po)
PATH_LCP     <- file.path(DIR_BASE, "LCP_clip.gpkg")                      # LCP land cover
PATH_DTM     <- file.path(DIR_BASE, "DTM_merge.tif")                      # extended DTM (reliefs)
PATH_DIST_OLD<- file.path(DIR_BASE, "output_spettrale/STEP3_output/dist_fiume.tif")  # Vercelli

# Parameters ----
AREA_MIN_PO <- 20      # ha: threshold for the Po connected components in LCP (channels < 1 ha)
LCP_CAMPO   <- "clc_3liv"   # class field in LCP
LCP_PO_VAL  <- 511          # "5.1.1 Watercourses"

# Reference grid = an extended composite (extent/resolution/CRS) ----
# look for a whole-park composite; if none, stop with clear instructions
trova_griglia <- function() {
  cand <- list.files(DIR_EXT, "mediana_10bande.*\\.tif$", full.names = TRUE)
  if (!length(cand))
    cand <- list.files(file.path(DIR_BASE, "output_spettrale"),
                       "mediana_10bande_PARCO\\.tif$", full.names = TRUE)
  cand
}
comp <- trova_griglia()
if (!length(comp)) stop(
  "No extended composite found.\n",
  "  Generate the 4 whole-park composites first (they define the grid),\n",
  "  then re-run this script. Looked in:\n   - ", DIR_EXT, "\n   - ",
  file.path(DIR_BASE, "output_spettrale"), " (pattern _PARCO)")
griglia <- rast(comp[1])[[1]]            # one band: only the geometry is needed
griglia[] <- NA
cat("Topography + dist_fiume (whole park)\n")
cat("Reference grid:", basename(comp[1]), "\n")
cat(sprintf("  extent: x[%.0f, %.0f]  y[%.0f, %.0f]  res %.0f m  CRS %s\n",
            ext(griglia)$xmin, ext(griglia)$xmax, ext(griglia)$ymin, ext(griglia)$ymax,
            res(griglia)[1], crs(griglia, describe=TRUE)$code))

# River distance from the hydrographic network (Po axis) ----
cat("\n1. dist_fiume from the hydrographic network (asse_Po.gpkg)\n")
asse <- vect(PATH_ASSE_PO)
asse <- project(asse, crs(griglia))
r_asse <- rasterize(asse, griglia, field = 1, background = NA)
dist_rete <- terra::distance(r_asse)
names(dist_rete) <- "dist_fiume"
writeRaster(dist_rete, file.path(DIR_OUT, "dist_fiume_rete.tif"), overwrite = TRUE)
cat("   saved: dist_fiume_rete.tif\n")

# River distance from LCP (Po isolated from channels: connected components > threshold) ----
cat("\n2. dist_fiume from LCP (Po isolated, clc_3liv==511)\n")
lcp <- st_read(PATH_LCP, quiet = TRUE)
if (!LCP_CAMPO %in% names(lcp)) stop("Field '", LCP_CAMPO, "' not found in LCP. Fields: ",
                                     paste(names(lcp), collapse=", "))
fiume <- lcp[lcp[[LCP_CAMPO]] == LCP_PO_VAL, ]
fiume <- st_transform(fiume, crs(griglia))
# stitch adjacent fragments and keep only the large components (= the Po, not the channels)
fiume_u <- st_union(st_buffer(fiume, 5))
parts   <- st_cast(fiume_u, "POLYGON")
parts   <- st_sf(geometry = parts)
parts$area_ha <- as.numeric(st_area(parts)) / 1e4
po_lcp  <- parts[parts$area_ha >= AREA_MIN_PO, ]
cat(sprintf("   511 components: %d -> Po components (>=%d ha): %d  (%.0f ha)\n",
            nrow(parts), AREA_MIN_PO, nrow(po_lcp), sum(po_lcp$area_ha)))
r_po_lcp  <- rasterize(vect(po_lcp), griglia, field = 1, background = NA)
dist_lcp  <- terra::distance(r_po_lcp)
names(dist_lcp) <- "dist_fiume"
writeRaster(dist_lcp, file.path(DIR_OUT, "dist_fiume_lcp.tif"), overwrite = TRUE)
cat("   saved: dist_fiume_lcp.tif\n")

# Consistency check on Vercelli (d_new - d_old) ----
# decides whether the models stay valid without retraining
cat("\n3. Vercelli check (comparison with the old dist_fiume.tif)\n")
if (file.exists(PATH_DIST_OLD)) {
  d_old <- rast(PATH_DIST_OLD)
  verifica <- function(d_new, etichetta) {
    dn <- resample(crop(d_new, ext(d_old)), d_old, method = "bilinear")
    df <- dn - d_old
    ad <- abs(df)
    med <- as.numeric(global(ad, median, na.rm = TRUE))
    p95 <- as.numeric(global(ad, \(x) quantile(x, 0.95, na.rm = TRUE)))
    p0  <- as.numeric(global(ad > 10, "mean", na.rm = TRUE)) * 100
    cat(sprintf("   [%s vs old]  |delta| median %.1f m | 95th pct %.1f m | pixels |delta|>10 m: %.1f%%\n",
                etichetta, med, p95, p0))
    writeRaster(df, file.path(DIR_OUT, paste0("verifica_diff_", etichetta, ".tif")), overwrite = TRUE)
    invisible(c(med = med, p95 = p95, pct10 = p0))
  }
  v_rete <- verifica(dist_rete, "rete")
  v_lcp  <- verifica(dist_lcp,  "lcp")
  cat("\n   Reading: if |delta| is small almost everywhere (e.g. <10 m on >95% of\n")
  cat("   pixels), the Vercelli dist_fiume does not change -> keep the models without\n")
  cat("   retraining. The 'lcp' version (channel-polygon) should match the old one\n")
  cat("   better (same concept); 'rete' (axis-line) differs inside/near the channel.\n")
  cat("   Difference maps saved (verifica_diff_*.tif) for inspection in QGIS.\n")
} else {
  cat("   ! Old dist_fiume not found:", PATH_DIST_OLD, "-> check skipped.\n")
}

# Topography from DTM (dem, slope, aspect -> sin/cos) ----
cat("\n4. Topography from DTM_merge.tif\n")
dtm <- rast(PATH_DTM)
if (crs(dtm) != crs(griglia)) dtm <- project(dtm, crs(griglia))
dem <- resample(dtm, griglia, method = "bilinear")
names(dem) <- "dem"
writeRaster(dem, file.path(DIR_OUT, "dem.tif"), overwrite = TRUE)

slope <- terrain(dem, "slope", unit = "degrees")
names(slope) <- "slope"
writeRaster(slope, file.path(DIR_OUT, "slope.tif"), overwrite = TRUE)

asp <- terrain(dem, "aspect", unit = "radians")
aspect_sin <- sin(asp); names(aspect_sin) <- "aspect_sin"
aspect_cos <- cos(asp); names(aspect_cos) <- "aspect_cos"
writeRaster(aspect_sin, file.path(DIR_OUT, "aspect_sin.tif"), overwrite = TRUE)
writeRaster(aspect_cos, file.path(DIR_OUT, "aspect_cos.tif"), overwrite = TRUE)
cat("   saved: dem.tif, slope.tif, aspect_sin.tif, aspect_cos.tif\n")

cat("\nDone. Output in:", DIR_OUT, "\n")
cat("Note: the final river distance is produced by the mosaic step (Vercelli\n")
cat("dist_fiume unchanged + Turin from the hydrographic network).\n")

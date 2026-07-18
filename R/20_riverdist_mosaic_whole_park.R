# dist_fiume mosaic (Vercelli unchanged + Turin from the asse_Po network).
# Builds the whole-park dist_fiume WITHOUT altering Vercelli:
#   Vercelli -> the old dist_fiume.tif as is (|delta|=0 -> models unchanged)
#   Turin    -> distance from the Po axis (hydrographic network asse_Po.gpkg), the
#               source most consistent with the old one (|delta| median ~95 m, vs
#               528/898 for the LCP polygons) and continuous along the whole course
#               (182 km).
#   cover(old, network): Vercelli takes precedence; the network only fills Turin.
#
# Safeguards: the old Vercelli dist_fiume.tif is READ-ONLY, never written. Output
#   goes only to output_spettrale_intero_parco/STEP3_output/. The output is named
#   "dist_fiume_mosaico.tif" (NOT "dist_fiume.tif"): the final name is chosen by
#   hand after checking the seam. No existing file is deleted.
#
# Note (declarable in the thesis): in Turin the distance is from the Po axis-line,
# in Vercelli from the channel edge (land cover 5111). It is a slight change of
# criterion on a "soft" predictor (what matters is the near/far order of magnitude
# from the river), acceptable in an extrapolation context. Vercelli, the validated
# part, stays identical to the pixel.

# Packages ----
library(terra)
library(sf)

# Folders ----
# Set the project root directory here, or via the POPARK_DATA environment variable
DIR_BASE <- Sys.getenv("POPARK_DATA", unset = "path/to/project/root")
if (!dir.exists(DIR_BASE))
  stop("Invalid DIR_BASE: ", DIR_BASE,
       "
     Set the path at the top of the script, or the POPARK_DATA variable.")
PATH_ASSE_PO  <- file.path(DIR_BASE, "asse_Po.gpkg")                                   # network (Turin)
PATH_DIST_OLD <- file.path(DIR_BASE, "output_spettrale/STEP3_output/dist_fiume.tif")   # READ-ONLY
DIR_EXT       <- file.path(DIR_BASE, "output_spettrale_intero_parco")
DIR_OUT       <- file.path(DIR_EXT, "STEP3_output")
PATH_OUT      <- file.path(DIR_OUT, "dist_fiume_mosaico.tif")
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

# Safeguard: never write over the old dist_fiume ----
if (normalizePath(PATH_OUT, mustWork = FALSE) == normalizePath(PATH_DIST_OLD, mustWork = FALSE))
  stop("STOP: output equals the old dist_fiume. Aborted for safety.")
if (!file.exists(PATH_DIST_OLD)) stop("Old dist_fiume not found: ", PATH_DIST_OLD)
if (!file.exists(PATH_ASSE_PO))  stop("asse_Po not found: ", PATH_ASSE_PO)

# Reference grid = an extended composite ----
comp <- list.files(DIR_EXT, "mediana_10bande\\.tif$", full.names = TRUE)
if (!length(comp)) stop("No extended composite in ", DIR_EXT)
griglia <- rast(comp[1])[[1]]; griglia[] <- NA
cat("dist_fiume MOSAIC (Vercelli unchanged + Turin from the network)\n")
cat("Grid:", basename(comp[1]), sprintf("(%d x %d, %dm)\n", nrow(griglia), ncol(griglia), res(griglia)[1]))

# 1. Distance from the Po axis (network) over the whole grid ----
cat("\n1. Distance from the Po axis (network) over the whole grid\n")
asse <- project(vect(PATH_ASSE_PO), crs(griglia))
r_asse <- rasterize(asse, griglia, field = 1, background = NA)
dist_rete <- terra::distance(r_asse)
names(dist_rete) <- "dist_fiume"

# 2. Old Vercelli dist_fiume (read-only), aligned to the grid ----
cat("2. Old Vercelli dist_fiume (read-only), aligning to the grid\n")
d_old <- rast(PATH_DIST_OLD)
old_ext <- extend(d_old, griglia)                       # adds NA cells outside Vercelli
if (!compareGeom(old_ext, griglia, stopOnError = FALSE))
  old_ext <- resample(d_old, griglia, method = "near")  # fallback if not perfectly aligned

# 3. Mosaic: Vercelli (old) takes precedence; the network fills Turin ----
cat("3. cover(old_Vercelli, network_Turin)\n")
dist_mosaico <- cover(old_ext, dist_rete)               # x where non-NA, else y
names(dist_mosaico) <- "dist_fiume"
writeRaster(dist_mosaico, PATH_OUT, overwrite = TRUE)   # whole-park folder only
cat("   saved:", PATH_OUT, "\n")

# 4. Check: on Vercelli the mosaic MUST equal the old raster ----
cat("\n4. Check: mosaic vs old on the Vercelli window (must be 0)\n")
m_vc <- crop(dist_mosaico, ext(d_old))
if (!compareGeom(m_vc, d_old, stopOnError = FALSE)) m_vc <- resample(m_vc, d_old, "near")
dmax <- as.numeric(global(abs(m_vc - d_old), "max", na.rm = TRUE))
cat(sprintf("   |delta|max Vercelli = %.6g  %s\n", dmax,
            if (is.finite(dmax) && dmax < 1e-6) "OK (Vercelli identical -> models unchanged)" else "!! WARNING"))

# 5. Seam: size of the jump at the western Vercelli border (info) ----
cat("\n5. Vercelli/Turin seam: size of the gap at the border\n")
xw <- as.numeric(ext(d_old)$xmin)
strip <- ext(xw, xw + 300, ext(d_old)$ymin, ext(d_old)$ymax)   # 300 m strip inside Vercelli
old_s  <- crop(d_old, strip)
rete_s <- resample(crop(dist_rete, strip), old_s, "bilinear")
seam <- as.numeric(global(abs(old_s - rete_s), "mean", na.rm = TRUE))
cat(sprintf("   mean old-vs-network gap near the western border: ~%.0f m\n", seam))
cat("   (the 'axis-line vs channel-polygon' discontinuity; acceptable on a soft predictor)\n")

cat("\nDone. Old Vercelli dist_fiume intact. No overwrite.\n")
cat("   Produced:", PATH_OUT, "\n")
cat("   If check (4) is OK, rename 'dist_fiume_mosaico.tif' to 'dist_fiume.tif' in the\n")
cat("   whole-park folder for the downstream steps.\n")
cat("   Inspect 'dist_fiume_mosaico.tif' in QGIS for a visual seam check.\n")

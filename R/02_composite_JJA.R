# Seasonal median Sentinel-2 L2A composite - summer (JJA, Jun-Aug 2019)
# Reference sector: Parco del Po piemontese, Vercelli-Alessandria (training area).
# Tiles 32TLQ/32TLR/32TMQ/32TMR, 10 bands, clouds/shadows/snow masked via SCL.

# Packages ----
library(terra)
library(gdalcubes)
gdalcubes_options(parallel = 3)   # threads for parallel processing

# Folders ----
# Set the project root directory here, or via the POPARK_DATA environment variable
DIR_BASE <- Sys.getenv("POPARK_DATA", unset = "path/to/project/root")
if (!dir.exists(DIR_BASE))
  stop("Invalid DIR_BASE: ", DIR_BASE,
       "
     Set the path at the top of the script, or the POPARK_DATA variable.")
dir_s2     <- file.path(DIR_BASE, "S2", "JJA")        # Sentinel-2 archives for the season
dir_output <- file.path(DIR_BASE, "output_spettrale")
dir.create(dir_output, showWarnings = FALSE)

# Find Sentinel-2 archives ----
files_zip <- c(
  list.files(dir_s2, pattern = "\\.zip$",  full.names = TRUE, recursive = FALSE),
  list.files(dir_s2, pattern = "\\.SAFE$", full.names = TRUE, recursive = FALSE)
)
cat("Sentinel-2 files found:", length(files_zip), "\n")

# check that the 4 tiles are present
tlq <- sum(grepl("T32TLQ", files_zip))
tlr <- sum(grepl("T32TLR", files_zip))
tmq <- sum(grepl("T32TMQ", files_zip))
tmr <- sum(grepl("T32TMR", files_zip))
cat(sprintf("Tiles 32TLQ:%d  32TLR:%d  32TMQ:%d  32TMR:%d\n", tlq, tlr, tmq, tmr))

# the western tiles (LQ/LR) are required for the area
if (tlq == 0 | tlr == 0) {
  stop("Missing files for tiles 32TLQ and/or 32TLR. Download the images from Copernicus and retry.")
}

# Image collection ----
# load the .db if it exists (fast), otherwise build it from the archives
file_collection <- file.path(dir_output, "S2_JJA_4tile_collection.db")
if (file.exists(file_collection)) {
  cat("Existing collection: loading it...\n")
  col_s2 <- image_collection(file_collection)
} else {
  cat("Collection not found: creating it (may take a few minutes)...\n")
  col_s2 <- create_image_collection(
    files           = files_zip,
    format          = "Sentinel2_L2A",
    unroll_archives = TRUE,
    out_file        = file_collection
  )
}
print(col_s2)

# Cloud mask (SCL) ----
# excluded SCL values (invalid pixels): 0 No data, 1 Saturated, 2 Dark area,
# 3 Cloud shadow, 7 Unclassified, 8-9 Cloud medium/high, 10 Thin cirrus, 11 Snow/Ice
maschera_nuvole <- image_mask("SCL", values = c(0, 1, 2, 3, 7, 8, 9, 10, 11))

# Data cube view ----
# extent = bbox of usosuolo_pulito (Vercelli sector); EPSG:32632 (UTM 32N, metres)
# single temporal window (dt = P3M): median over the whole season
vista <- cube_view(
  srs        = "EPSG:32632",
  extent     = list(
    t0     = "2019-06-01",         # season start
    t1     = "2019-08-31",         # season end
    left   = 428978,         # xmin from st_bbox(usosuolo_pulito)
    right  = 493372,         # xmax
    bottom = 4978856,        # ymin
    top    = 5004263         # ymax
  ),
  dx          = 10,          # 10 m resolution
  dy          = 10,
  dt          = "P3M",
  aggregation = "median",
  resampling  = "bilinear"
)
print(vista)

# Build the composite ----
# temporary folder for the per-band tiffs, cleaned beforehand to avoid leftovers
dir_temp <- file.path(DIR_BASE, "tmp_S2_JJA")
if (dir.exists(dir_temp)) unlink(dir_temp, recursive = TRUE)
dir.create(dir_temp, showWarnings = FALSE)

cat("Building composite...\n")
raster_cube(col_s2, vista, mask = maschera_nuvole) |>
  select_bands(c("B02","B03","B04","B05","B06","B07","B08","B8A","B11","B12")) |>
  write_tif(dir_temp)

# Merge bands into a single GeoTIFF ----
# write_tif produces one file per band: merge them into a 10-band GeoTIFF
files_bande <- list.files(dir_temp, pattern = "\\.tif$", full.names = TRUE)
cat("Band files produced:", length(files_bande), "\n")

composito <- rast(files_bande)
cat("NA pixels per band (expected 0):\n")
print(global(composito, fun = "isNA"))

file_composito <- file.path(dir_output,
  "S2_MSIL2A_JJA2019_T32TLQ-32TLR-32TMQ-32TMR_mediana_10bande.tif")
writeRaster(composito, file_composito, overwrite = TRUE)
cat("Composite saved to:", file_composito, "\n")

unlink(dir_temp, recursive = TRUE)   # remove the temporary folder

# Reload and check statistics ----
# reload from the final file and check per-band stats as a sanity check
composito <- rast(file_composito)
cat("Per-band statistics (scale 0-10000):\n")
print(global(composito, fun = c("min", "mean", "max"), na.rm = TRUE))
# expected: mean B08 (NIR) > mean B04 (red) for summer vegetation

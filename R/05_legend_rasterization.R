# Rasterisation of the legend (COD_L2 + COD_L1) for the Po park LULC mapping.
# Assigns each original CLC class its hierarchical COD_L2 (18 classes) and, by
# derivation, COD_L1 (7 macrotypes); applies the legend to the vector and saves
# the labelled GeoPackage; rasterises COD_L2 and COD_L1 on the 10 m Sentinel-2
# grid. COD_L1 is the first digit of COD_L2, so L1/L2 consistency holds by
# construction and is re-checked with an assertion.
# Rasterisation uses touches=FALSE: each pixel takes the class of the polygon
# covering its centre, avoiding over-representation of small border polygons.
# All classes enter training (no mask); mixed/anomalous pixels are filtered
# downstream by the z-statistics and the physical reflectance filter.

# Packages ----
library(terra)
library(sf)
library(dplyr)

# Folders ----
# Set the project root directory here, or via the POPARK_DATA environment variable
DIR_BASE <- Sys.getenv("POPARK_DATA", unset = "path/to/project/root")
if (!dir.exists(DIR_BASE))
  stop("Invalid DIR_BASE: ", DIR_BASE,
       "
     Set the path at the top of the script, or the POPARK_DATA variable.")
DIR_OUTPUT <- file.path(DIR_BASE, "output_spettrale")
DIR_STEP3  <- file.path(DIR_OUTPUT, "STEP3_output")

PATH_VETTORIALE <- file.path(DIR_BASE, "usosuolo_2019_layout_def_CLEAN.gpkg")
PATH_COMPOSITO  <- file.path(DIR_OUTPUT,
  "S2_MSIL2A_MAM2019_T32TLQ-32TLR-32TMQ-32TMR_mediana_10bande.tif")

if (!dir.exists(DIR_STEP3)) dir.create(DIR_STEP3, recursive = TRUE)

# local tempdir (RAM safety on 16 GB)
tmpdir <- file.path(DIR_STEP3, "temp")
dir.create(tmpdir, showWarnings = FALSE)
terraOptions(memfrac = 0.6, tempdir = tmpdir)

cat("Legend rasterisation. Output:", DIR_STEP3, "\n\n")

for (f in c(PATH_VETTORIALE, PATH_COMPOSITO))
  if (!file.exists(f)) stop("File not found:\n  ", f)
cat("Input files verified.\n\n")

# Legend definition (CLC class -> COD_L2 -> COD_L1) ----
# 'classe_to_l2' maps each original CLC code to its target COD_L2. Entries are
# grouped by COD_L2 (11, 12, 13...) only for readability; order does not matter.
# COD_L1 and the descriptions are derived right after from lookup tables.
classe_to_l2 <- c(
  "5111"=11, "5112"=11, "5113"=11, "5121"=11,
  "4241"=12, "5122"=12, "5123"=12, "1312"=12, "1313"=12,  # 1312/1313 reassigned from 71 to 12 (quarry water bodies = still water)
  "4113"=13, "4242"=13, "4243"=13, "4244"=13,
  "3311"=21, "3312"=21,
  "332"=22,  "335"=22,  "3241"=22,
  "3212"=31, "3242"=31, "3243"=31, "3246"=31,
  "3111"=41, "3112"=41,
  "3113"=42, "3116"=42,
  "3115"=43,
  "224"=44,
  "225"=45,  "3117"=45, "3244"=45,
  "3114"=46,
  "213"=51,
  "226"=52,  "2111"=52, "2112"=52, "2113"=52, "2121"=52,
  "221"=53,  "222"=53,
  "2311"=54,
  "111"=61,  "112"=61,  "121"=61,  "122"=61, "142"=61,
  "1321"=61, "1421"=61, "1422"=61,
  "1311"=71, "1314"=71, "1315"=71, "1322"=71
)

# L1 / L2 descriptions (kept in Italian: written to the legend table and GeoPackage)
desc_l1 <- c(
  "1"="Acque e zone umide", "2"="Substrato fluviale",
  "3"="Vegetazione erbacea e arbustiva", "4"="Copertura arborea",
  "5"="Superfici coltivate", "6"="Superfici impermeabilizzate",
  "7"="Aree estrattive"
)
desc_l2 <- c(
  "11"="Acque interne",
  "12"="Lanche e acque ferme fluviali",
  "13"="Vegetazione igrofila erbacea",
  "21"="Greto e sabbie fluviali",
  "22"="Formazioni pioniere su substrato fluviale",
  "31"="Vegetazione erbacea alta e arbustiva ripariale",
  "41"="Boschi ripariali igrofili",
  "42"="Boschi planiziali di latifoglie mesofite",
  "43"="Robinieti",
  "44"="Pioppeti",
  "45"="Piantagioni arboree da legno",
  "46"="Boschi misti collinari",
  "51"="Risaie",
  "52"="Seminativi e colture erbacee annuali",
  "53"="Colture arboree permanenti",
  "54"="Prati stabili",
  "61"="Tessuto urbano e aree industriali",
  "71"="Aree estrattive e cave"
)

# Build the legend table: one row per CLC class with COD_L2/COD_L1 and
# descriptions. COD_L1 = first digit of COD_L2; descriptions come from the lookups.
tab_classi <- data.frame(
  CLASSE = as.integer(names(classe_to_l2)),
  COD_L2 = as.integer(classe_to_l2),
  stringsAsFactors = FALSE
)
tab_classi$COD_L1    <- tab_classi$COD_L2 %/% 10          # first digit of COD_L2
tab_classi$DESC_L1   <- desc_l1[as.character(tab_classi$COD_L1)]
tab_classi$DESC_L2   <- desc_l2[as.character(tab_classi$COD_L2)]
tab_classi$TIPO      <- "Training"                        # no mask
tab_classi <- tab_classi[, c("CLASSE","COD_L1","DESC_L1",
                             "COD_L2","DESC_L2","TIPO")]

# Legend checks (stop before writing outputs) ----
# Stops if the counts are off (53/18/7), if there are duplicated CLASSE codes,
# missing descriptions, or if L1/L2 consistency is violated.
cat("Legend check:\n")
cat(sprintf("   Total classes : %d  (expected: 53)\n", nrow(tab_classi)))
cat(sprintf("   L2 classes    : %d  (expected: 18)\n", length(unique(tab_classi$COD_L2))))
cat(sprintf("   L1 macrotypes : %d  (expected:  7)\n", length(unique(tab_classi$COD_L1))))
cat(sprintf("   All Training  : %s\n", all(tab_classi$TIPO == "Training")))

if (nrow(tab_classi) != 53)                  stop("Classes != 53")
if (length(unique(tab_classi$COD_L2)) != 18) stop("L2 classes != 18")
if (length(unique(tab_classi$COD_L1)) != 7)  stop("L1 macrotypes != 7")
if (any(duplicated(tab_classi$CLASSE)))      stop("Duplicated CLASSE")
if (any(is.na(tab_classi$DESC_L2)))          stop("Missing DESC_L2")
if (!all(tab_classi$COD_L2 %/% 10 == tab_classi$COD_L1))
  stop("L1/L2 consistency failed")
cat("   No duplicates, L1/L2 consistency OK\n\n")

write.csv(tab_classi, file.path(DIR_STEP3, "tabella_COD_L2.csv"), row.names = FALSE)
cat("   tabella_COD_L2.csv saved\n\n")

# Apply the legend to the vector ----
# Reads the polygons, aligns the CRS to the grid, cleans geometries, joins the
# legend (on CLASSE) and saves the labelled GeoPackage. Each step has a check
# that stops the script if something is off (CRS, geometries, classes).
cat("Updating vector...\n")

uso_sf <- st_read(PATH_VETTORIALE, quiet = TRUE)
cat("   Polygons read:", nrow(uso_sf),
    "| Unique classes (CLASSE):", length(unique(uso_sf$CLASSE)), "\n")

# reference S2 grid (CRS EPSG:32632, 10 m)
griglia <- rast(PATH_COMPOSITO)[[1]]

# align vector CRS to the grid (robust comparison via st_crs)
crs_griglia_sf <- st_crs(crs(griglia, proj = TRUE))
if (!isTRUE(st_crs(uso_sf) == crs_griglia_sf)) {
  cat("   Different CRS -> reprojecting vector...\n")
  uso_sf <- st_transform(uso_sf, crs = crs_griglia_sf)
} else {
  cat("   CRS match\n")
}

# geometry cleaning: drop empty, fix invalid
n_empty <- sum(st_is_empty(uso_sf))
if (n_empty > 0) {
  cat("   !", n_empty, "empty geometries removed\n")
  uso_sf <- uso_sf[!st_is_empty(uso_sf), ]
}
n_invalid <- sum(!st_is_valid(uso_sf))
if (n_invalid > 0) {
  cat("   !", n_invalid, "invalid geometries -> st_make_valid()...\n")
  uso_sf <- st_make_valid(uso_sf)
}
cat("   Polygons after geometry cleaning:", nrow(uso_sf), "\n")

# drop any legend fields left from previous runs
campi_vecchi <- c("COD_L1","DESC_L1","COD_L2","DESC_L2","TIPO","TIPO_ORIG","COD_MASCHERA")
uso_sf <- uso_sf[, !names(uso_sf) %in% campi_vecchi]

# CLASSE as integer (robust join even if the gpkg read it as text or factor)
uso_sf$CLASSE <- as.integer(as.character(uso_sf$CLASSE))

# join the legend on CLASSE
uso_sf <- uso_sf %>% left_join(tab_classi, by = "CLASSE")

# check: no polygon without COD_L2 (i.e. CLASSE missing from the legend)
n_na <- sum(is.na(uso_sf$COD_L2))
if (n_na > 0) {
  mancanti <- sort(unique(uso_sf$CLASSE[is.na(uso_sf$COD_L2)]))
  stop("WARNING: ", n_na, " polygons with CLASSE not in the legend: ",
       paste(mancanti, collapse = ", "))
}
cat("   All polygons have COD_L2 assigned\n")
cat(sprintf("   Training polygons (all classes): %d\n",
    sum(uso_sf$TIPO == "Training")))

st_write(uso_sf, file.path(DIR_STEP3, "usosuolo_2019_COD_L2.gpkg"),
         delete_dsn = TRUE, quiet = TRUE)
cat("   usosuolo_2019_COD_L2.gpkg saved\n\n")

# Rasterisation (COD_L2, then COD_L1 by derivation) ----
# touches=FALSE: each pixel takes the COD_L2 of the polygon covering its centre,
# without over-representing small border polygons. The COD_L1 raster is obtained
# by reclassifying COD_L2 with its first digit.
cat("Rasterising...\n")

rast_l2 <- rasterize(vect(uso_sf), griglia,
                     field = "COD_L2", background = NA, touches = FALSE)

freq_l2     <- as.data.frame(freq(rast_l2))
valori_l2   <- sort(freq_l2$value)
n_px        <- sum(freq_l2$count)
l2_attesi   <- c(11,12,13, 21,22, 31, 41,42,43,44,45,46, 51,52,53,54, 61, 71)
l2_mancanti <- setdiff(l2_attesi, valori_l2)

cat("   COD_L2 in the raster:", paste(valori_l2, collapse=", "), "\n")
cat("   Total pixels:", format(n_px, big.mark="."), "\n")
if (length(l2_mancanti) > 0)
  cat("   ! Missing L2 classes:", paste(l2_mancanti, collapse=", "),
      "- polygons out of extent or too small?\n") else
  cat("   All 18 L2 classes present\n")

writeRaster(rast_l2, file.path(DIR_STEP3, "raster_COD_L2.tif"),
            datatype = "INT2S", overwrite = TRUE,
            gdal = c("COMPRESS=LZW", "TILED=YES"))
cat("   raster_COD_L2.tif written\n")

# raster_COD_L1: reclassify (first digit of COD_L2 = COD_L1)
l2_unici <- sort(unique(tab_classi$COD_L2))
rcl_l1   <- cbind(l2_unici, l2_unici %/% 10)
rast_l1  <- classify(rast_l2, rcl_l1, others = NA)

valori_l1 <- sort(as.data.frame(freq(rast_l1))$value)
cat("   COD_L1 in the raster:", paste(valori_l1, collapse=", "),
    if (length(setdiff(1:7, valori_l1)) == 0) "- all 7 present\n" else "! MISSING\n")

writeRaster(rast_l1, file.path(DIR_STEP3, "raster_COD_L1.tif"),
            datatype = "INT1U", overwrite = TRUE,
            gdal = c("COMPRESS=LZW", "TILED=YES"))
cat("   raster_COD_L1.tif written\n\n")

# Report: pixels per L2 class ----
cat("Pixels per L2 class:\n")
rep_l2 <- merge(
  data.frame(COD_L2 = valori_l2,
             N_pixel = freq_l2$count[match(valori_l2, freq_l2$value)]),
  unique(tab_classi[, c("COD_L2","COD_L1","DESC_L2")]),
  by = "COD_L2"
)
rep_l2 <- rep_l2[order(rep_l2$COD_L2), ]
for (i in seq_len(nrow(rep_l2)))
  cat(sprintf("   L1=%d  L2=%-2d  %-46s %s px\n",
      rep_l2$COD_L1[i], rep_l2$COD_L2[i],
      substr(rep_l2$DESC_L2[i], 1, 46),
      format(rep_l2$N_pixel[i], big.mark = ".")))

cat("\nDone: raster_COD_L2.tif + raster_COD_L1.tif + usosuolo_2019_COD_L2.gpkg\n")

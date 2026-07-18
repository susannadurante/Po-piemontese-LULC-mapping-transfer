# Pixel extraction + z-statistics + pure-pixel selection.
#
# Part A - extracts predictor values for every training pixel: 40 S2 bands
#   (10 x 4 seasons, reflectance 0-1), 44 spectral indices (11 x 4), 5 topographic
#   (dist_fiume, dem, slope, aspect_sin, aspect_cos) and 16 GLCM (4 x 4). All
#   layers are resampled onto r_ref (MAM composite, 10 m) before stacking:
#   method="near" for bands/indices/GLCM, method="bilinear" for continuous
#   surfaces (dem, slope, dist_fiume, aspect_raw).
# Part B - physical filter (B02 > 0.40 reflectance = residual cloud/nodata), then
#   drop_na on the 101 z-statistics predictors only. GLCM_COR (structural NAs on
#   uniform surfaces) are kept in the CSV for traceability but excluded from the
#   z-statistics and from the RF (ranger does not handle NAs).
# Part C - per-class z-statistics: z_total = sqrt(sum(((x-mu)/sigma)^2)) over 101
#   predictors, with an adaptive per-class window mu_z +/- sigma_z (a pixel is
#   pure if z_total falls in that range). Classes with < 50 pixels keep all
#   pixels (too small to estimate the z distribution).
# Part D - writes training_pixels_clean.csv, z_statistics_summary.csv and
#   pixel_per_classe_confronto.csv.
#
# aspect is a circular variable (0 and 360 deg are the same direction), decomposed
# into sin/cos before the z-statistics (Zar 1999).

# Packages ----
library(terra)
library(sf)
library(dplyr)
library(tidyr)

# explicit namespace conflicts: dplyr::select/filter can be masked by terra or
# MASS; terra::extract can be masked by raster or tidyr.
select  <- dplyr::select
filter  <- dplyr::filter
extract <- terra::extract

# allocate up to 80% of RAM for raster operations (reduces crash risk during
# extract() on large stacks). terra's default is 0.6; raising it is safe on a
# dedicated machine.
terraOptions(memfrac = 0.8)

# Folders ----
# Set the project root directory here, or via the POPARK_DATA environment variable
DIR_BASE <- Sys.getenv("POPARK_DATA", unset = "path/to/project/root")
if (!dir.exists(DIR_BASE))
  stop("Invalid DIR_BASE: ", DIR_BASE,
       "
     Set the path at the top of the script, or the POPARK_DATA variable.")
DIR_OUT   <- file.path(DIR_BASE, "output_spettrale", "STEP4_output")
DIR_STEP3 <- file.path(DIR_BASE, "output_spettrale", "STEP3_output")
DIR_COMP  <- file.path(DIR_BASE, "output_spettrale")
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

# training polygons from the legend step (GeoPackage with TIPO, COD_L1, COD_L2).
# COD_L1 and COD_L2 are already fields in the GeoPackage: used directly, no join.
SHP_TRAINING <- file.path(DIR_STEP3, "usosuolo_2019_COD_L2.gpkg")

# Seasons, bands and indices ----
stagioni <- c("MAM", "JJA", "SON", "DJF")
bande    <- c("B02","B03","B04","B05","B06","B07","B08","B8A","B11","B12")

# seasonal composite naming:
# S2_MSIL2A_<season>2019_T32TLQ-32TLR-32TMQ-32TMR_mediana_10bande.tif
nome_composito <- function(s) {
  file.path(DIR_COMP,
    paste0("S2_MSIL2A_", s, "2019_T32TLQ-32TLR-32TMQ-32TMR_mediana_10bande.tif"))
}

# 11 spectral indices x 4 seasons = 44 predictors, named <INDEX>_<season>.tif
indici <- c("NDVI","NDRE","NDWI_McF","NDWI_Gao","MNDWI",
            "NDBI","NBR2","EVI2","BSI","GLI","GVMI")

# 4 GLCM metrics x 4 seasons = 16 predictors, named GLCM_<metric>_<season>_v3.tif
glcm_metriche <- c("CON","HOM","COR","ENT")

# Reference raster ----
# r_ref is the MAM composite at 10 m: it defines the grid onto which every other
# layer is resampled before extraction.
r_ref <- rast(nome_composito("MAM"))
names(r_ref) <- paste0(bande, "_MAM")
cat(sprintf("Reference raster: %d rows x %d cols, res %.1f m\n",
    nrow(r_ref), ncol(r_ref), res(r_ref)[1]))

# Load and resample all rasters ----
# Every layer is resampled onto r_ref before being added to the list, so the
# final stack has identical grid, resolution and extent for each layer.
# Resample strategy: near for bands/indices/GLCM (categorical or texture values,
# not interpolable), bilinear for dem/slope/dist_fiume/aspect_raw (continuous).
cat("\nPart A: loading and resampling predictor rasters\n")

layer_list <- list()  # each element is a single-band SpatRaster

# S2 bands (10 bands x 4 seasons = 40 layers). MAM is already r_ref.
for (s in stagioni) {
  if (s == "MAM") {
    r <- r_ref
  } else {
    r <- rast(nome_composito(s))
    names(r) <- paste0(bande, "_", s)
    r <- resample(r, r_ref, method = "near")
  }
  for (b in bande) {
    nm <- paste0(b, "_", s)
    layer_list[[nm]] <- r[[nm]]
  }
  cat(sprintf("   Bands %s: loaded and resampled (%d bands)\n", s, length(bande)))
}

# spectral indices (11 x 4 = 44 layers), named <INDEX>_<season>.tif in DIR_STEP3
for (s in stagioni) {
  for (idx in indici) {
    nm <- paste0(idx, "_", s)
    r  <- rast(file.path(DIR_STEP3, paste0(idx, "_", s, ".tif")))
    names(r) <- nm
    layer_list[[nm]] <- resample(r, r_ref, method = "near")
  }
}
cat(sprintf("   Indices: %d layers loaded and resampled\n", length(indici) * length(stagioni)))

# topographic: dist_fiume, dem, slope, aspect (4 layers). aspect is loaded as
# aspect_raw and decomposed into sin/cos after extraction. dem/slope come from the
# 5 m Piedmont DTM, already resampled to 10 m; here bilinear is a near-identity but
# kept for consistency. For aspect_raw bilinear on a circular variable is
# theoretically imprecise (values near 0 and 360 deg), but negligible here (near-
# identity resample, few 0/360 pixels in the plain); the sin/cos decomposition
# after extraction resolves the circular issue.
topo_bilinear <- list(
  dist_fiume = file.path(DIR_STEP3, "dist_fiume.tif"),
  dem        = file.path(DIR_STEP3, "dem.tif"),
  slope      = file.path(DIR_STEP3, "slope.tif"),
  aspect_raw = file.path(DIR_STEP3, "aspect.tif")
)
for (nm in names(topo_bilinear)) {
  r <- rast(topo_bilinear[[nm]])
  names(r) <- nm
  layer_list[[nm]] <- resample(r, r_ref, method = "bilinear")
}
cat(sprintf("   Topographic: %d layers loaded and resampled (bilinear)\n", length(topo_bilinear)))

# GLCM (4 metrics x 4 seasons = 16 layers), named GLCM_<metric>_<season>_v3.tif
# (_v3 = computed with a B02>4000 DN cloud mask to exclude contaminated pixels)
for (s in stagioni) {
  for (m in glcm_metriche) {
    nm <- paste0("GLCM_", m, "_", s)
    r  <- rast(file.path(DIR_STEP3, paste0(nm, "_v3.tif")))
    names(r) <- nm
    layer_list[[nm]] <- resample(r, r_ref, method = "near")
  }
}
cat(sprintf("   GLCM: %d layers loaded and resampled\n",
    length(glcm_metriche) * length(stagioni)))

# single stack: all layers share the grid now -> c() is safe. terra::rast()
# forces the conversion to SpatRaster.
all_rast <- terra::rast(layer_list)
cat(sprintf("\n   Total stack: %d layers\n", nlyr(all_rast)))
# expected: 40 bands + 44 indices + 4 topo (aspect_raw included) + 16 GLCM = 104
# (aspect_raw becomes aspect_sin + aspect_cos after extraction -> 105 predictors)

# Load training polygons ----
cat("\nLoading training polygons\n")
poly_train <- st_read(SHP_TRAINING, quiet = TRUE)

# training polygons = those with a defined COD_L2 (non-NA). This is the substantive
# criterion, more robust than filtering on the textual TIPO field. In the final
# legend all classes have COD_L2 defined, so no valid polygon is lost here.
poly_train <- poly_train %>% filter(!is.na(COD_L2))

# reproject the vector to the raster CRS (EPSG:32632 for the Po)
poly_train <- st_transform(poly_train, crs(r_ref))

cat(sprintf("   Training polygons: %d\n", nrow(poly_train)))
cat(sprintf("   Unique L2 classes: %d\n", length(unique(poly_train$COD_L2))))

# Pixel extraction ----
# batch extraction per polygon block: avoids loading the whole stack into RAM in a
# single call (which hangs R on large areas with many layers). Slower for big
# classes but keeps RAM bounded.
poly_vect <- vect(poly_train)

# add a row ID to the vector for the later join
poly_train <- poly_train %>% mutate(ID = row_number())
poly_vect  <- vect(poly_train)

cat("   Extracting (batches of polygons)...\n")

batch_size <- 200
n_poly     <- nrow(poly_train)
n_batch    <- ceiling(n_poly / batch_size)
batch_list <- list()

for (i in seq_len(n_batch)) {
  idx_start <- (i - 1) * batch_size + 1
  idx_end   <- min(i * batch_size, n_poly)
  poly_sub  <- poly_vect[idx_start:idx_end, ]

  df_b <- extract(all_rast, poly_sub, xy = TRUE, ID = TRUE)
  # ID in df_b is relative to the batch (1-based in the sub) -> map to the global ID
  df_b$ID <- df_b$ID + idx_start - 1
  batch_list[[i]] <- df_b

  if (i %% 10 == 0 || i == n_batch)
    cat(sprintf("   Batch %d/%d done (%d polygons)\n", i, n_batch, idx_end))
}

df_raw <- dplyr::bind_rows(batch_list)

# attach COD_L2 and COD_L1 to each pixel via the polygon ID. Codes are taken
# directly from the GeoPackage fields (not from a join on COD_L2, which would fail
# on COD_L2=NA). In the final legend every class has COD_L2 defined.
poly_df <- poly_train %>%
  st_drop_geometry() %>%
  select(ID, COD_L1, COD_L2)

df_raw <- df_raw %>%
  left_join(poly_df, by = "ID") %>%
  rename(ID_PIXEL = ID)

cat(sprintf("   Total pixels extracted: %s\n", format(nrow(df_raw), big.mark = ".")))
cat(sprintf("   Data frame columns: %d\n", ncol(df_raw)))

# Convert bands to reflectance 0-1 ----
# Sentinel-2 composites are in DN scaled x10000. z-statistics and the RF work on
# reflectance 0-1, for consistency with the spectral indices (already in 0-1).
cat("\nConverting S2 bands to reflectance 0-1 (/10000)\n")

# column names in the right order: band outer, season inner (B02_MAM, B02_JJA...)
cols_bande <- as.vector(outer(bande, stagioni, paste, sep = "_"))
stopifnot(all(cols_bande %in% names(df_raw)))

df_raw <- df_raw %>%
  mutate(across(all_of(cols_bande), ~ . / 10000))

# Decompose aspect into sin/cos ----
# aspect (slope orientation) is circular: 0 and 360 deg are the same direction but
# arithmetically far apart. Raw aspect in the z-statistics would bias north-facing
# slopes (values near 0/360). sin/cos decomposition solves it (Zar 1999):
#   aspect_sin = sin(aspect * pi/180)  -> N-S component
#   aspect_cos = cos(aspect * pi/180)  -> E-W component
cat("Decomposing aspect into sin/cos (circular variable)\n")

df_raw <- df_raw %>%
  mutate(
    asp_rad    = aspect_raw * pi / 180,
    aspect_sin = sin(asp_rad),
    aspect_cos = cos(asp_rad)
  ) %>%
  select(-aspect_raw, -asp_rad)
# after this: aspect_raw removed, aspect_sin/aspect_cos added
# -> topographic predictors: dist_fiume, dem, slope, aspect_sin, aspect_cos (5)

# Part B - physical filter + drop_na ----
cat("\nPart B: physical filter + drop_na\n")

# physical filter: B02 (blue) reflectance > 0.40 in any season flags residual
# clouds or nodata; exclude the pixel if anomalous in AT LEAST ONE season.
cols_b02  <- paste0("B02_", stagioni)
n_pre_fis <- nrow(df_raw)
df_train  <- df_raw %>%
  filter(if_all(all_of(cols_b02), ~ . <= 0.40))
n_post_fis <- nrow(df_train)
cat(sprintf("   Physical filter B02<=0.40: %s -> %s pixels (-%s anomalies)\n",
    format(n_pre_fis,  big.mark = "."),
    format(n_post_fis, big.mark = "."),
    format(n_pre_fis - n_post_fis, big.mark = ".")))

# drop NAs on the z-statistics predictors only. drop_na(all_of(cols_zstat)) rather
# than na.omit(df_train): na.omit would drop pixels for NAs in GLCM_COR, which do
# not enter the z-statistics, losing valid pixels (especially for small classes).
# GLCM_COR NAs stay in the data frame; they are excluded from the RF at training,
# since ranger does not handle missing values.
cols_glcm_cor_b2 <- paste0("GLCM_COR_", stagioni)
cols_id_b2       <- c("ID_PIXEL", "x", "y", "COD_L2", "COD_L1")
cols_zstat_b2    <- setdiff(names(df_train), c(cols_id_b2, cols_glcm_cor_b2))

n_pre_na  <- nrow(df_train)
df_train  <- df_train %>% drop_na(all_of(cols_zstat_b2))
n_post_na <- nrow(df_train)
cat(sprintf("   drop_na (101 pred.):    %s -> %s pixels (-%s with NA)\n\n",
    format(n_pre_na,  big.mark = "."),
    format(n_post_na, big.mark = "."),
    format(n_pre_na - n_post_na, big.mark = ".")))

# Save the pre-z-statistics extraction (input for predictor selection) ----
# full data frame BEFORE the z-stat purity filter: all reference pixels with all
# predictors. It is the input for predictor selection (ROC-AUC + Boruta + Kendall),
# which works on COLUMNS (features) and is independent of the z-stat filter, which
# works on ROWS (pixels).
write.csv(df_train,
          file.path(DIR_OUT, "predittori_estratti_completi.csv"),
          row.names = FALSE)
cat(sprintf("Saved predittori_estratti_completi.csv: %s pixels x %d columns (Boruta input)\n",
            format(nrow(df_train), big.mark = "."), ncol(df_train)))

# Part C - z-statistics with adaptive mu +/- sigma threshold ----
cat("Part C: z-statistics - adaptive per-class mu +/- sigma threshold\n\n")

# columns that identify the pixel but are not predictors (extract with xy=TRUE
# returns lowercase "x" and "y")
cols_id <- c("ID_PIXEL", "x", "y", "COD_L2", "COD_L1")

# GLCM_COR (texture correlation) is excluded from the z-statistics: on uniform
# surfaces (water bodies, continuous impervious areas) correlation is constant ->
# SD = 0 -> z undefined. The 4 GLCM_COR columns stay in the CSV for traceability
# only; with structural NAs they are also excluded from the RF (ranger).
cols_glcm_cor <- paste0("GLCM_COR_", stagioni)

# z-statistics predictors: everything except IDs and GLCM_COR
# expected: 110 total - 5 IDs - 4 GLCM_COR = 101 predictors
cols_zstat <- setdiff(names(df_train), c(cols_id, cols_glcm_cor))
N_PRED     <- length(cols_zstat)
cat(sprintf("   Predictors in the z-statistics: %d\n", N_PRED))
cat(sprintf("   (105 total predictors - 4 GLCM_COR = %d)\n\n", N_PRED))

# class-by-class loop
classi_L2    <- sort(unique(df_train$COD_L2))
summary_list <- list()
clean_list   <- list()

for (cl in classi_L2) {

  df_cl <- df_train %>% filter(COD_L2 == cl)
  n_cl  <- nrow(df_cl)

  # classes with fewer than 50 pixels: too small to estimate the z distribution
  # reliably -> keep all pixels.
  if (n_cl < 50) {
    cat(sprintf("   L2=%s: %d pixels < 50 -> keep all\n", cl, n_cl))
    clean_list[[as.character(cl)]]   <- df_cl
    summary_list[[as.character(cl)]] <- data.frame(
      COD_L2 = cl, n_pre = n_cl, n_post = n_cl,
      mu_z = NA, sigma_z = NA, t1 = NA, t2 = NA,
      note = "< 50 pixels, all kept"
    )
    next
  }

  # predictor matrix for this class (pixels x predictors)
  mat <- as.matrix(df_cl[, cols_zstat])

  # per-predictor mean and SD over the class
  mu_pred <- colMeans(mat, na.rm = TRUE)
  sd_pred <- apply(mat, 2, sd, na.rm = TRUE)

  # predictors with SD = 0 are constant in the class: the division would give NaN;
  # the contribution to the z distance is forced to 0 (via is.nan -> 0 below).
  # The predictor is not dropped: if it does not vary, it does not discriminate.
  sd_pred[sd_pred == 0] <- NA  # -> NaN after sweep -> then -> 0

  # z-statistics per pixel:
  #   1. centre each predictor: (x_p - mu_p)
  #   2. standardise by SD:      / sigma_p
  #   3. square and sum over predictors: sum z^2
  #   4. square root: simplified (diagonal) Mahalanobis distance
  z_mat           <- sweep(mat, 2, mu_pred, "-")
  z_mat           <- sweep(z_mat, 2, sd_pred, "/")
  z_mat[is.nan(z_mat)] <- 0        # SD=0: contribution forced to 0
  z_totale        <- sqrt(rowSums(z_mat^2, na.rm = TRUE))

  # adaptive threshold from the class's own z distribution: mu_z is the centre,
  # sigma_z its spread. Keep pixels in (mu_z - sigma_z, mu_z + sigma_z), i.e.
  # those spectrally close to the class centre. Homogeneous classes -> small
  # sigma_z -> narrow window; heterogeneous classes -> wide window.
  mu_z    <- mean(z_totale)
  sigma_z <- sd(z_totale)
  t1      <- mu_z - sigma_z
  t2      <- mu_z + sigma_z

  idx_puri <- which(z_totale > t1 & z_totale < t2)
  df_puri  <- df_cl[idx_puri, ]

  cat(sprintf("   L2=%s: %d -> %d pure pixels (%.1f%%)  [t1=%.2f, t2=%.2f]\n",
      cl, n_cl, nrow(df_puri),
      100 * nrow(df_puri) / n_cl, t1, t2))

  clean_list[[as.character(cl)]]   <- df_puri
  summary_list[[as.character(cl)]] <- data.frame(
    COD_L2  = cl,
    n_pre   = n_cl,
    n_post  = nrow(df_puri),
    mu_z    = round(mu_z,    4),
    sigma_z = round(sigma_z, 4),
    t1      = round(t1,      4),
    t2      = round(t2,      4),
    note    = ""
  )
}

# merge all classes into single data frames
df_clean   <- dplyr::bind_rows(clean_list)
df_summary <- dplyr::bind_rows(summary_list)

# Part D - output ----
cat("\nPart D: writing output\n")

# training_pixels_clean.csv: all pure pixels of all L2 classes in one CSV. Not
# split by model here (the split is done by the train/validation step). COD_L1 is
# on every row, so the downstream per-model filter is trivial.
write.csv(df_clean,
          file.path(DIR_OUT, "training_pixels_clean.csv"),
          row.names = FALSE)
cat(sprintf("   training_pixels_clean.csv:      %s pixels x %d columns\n",
    format(nrow(df_clean), big.mark = "."), ncol(df_clean)))

# z_statistics_summary.csv: per-class thresholds and pixel counts before/after.
# Useful to document the selection in the thesis methodology.
write.csv(df_summary,
          file.path(DIR_OUT, "z_statistics_summary.csv"),
          row.names = FALSE)
cat("   z_statistics_summary.csv:       written\n")

# pixel_per_classe_confronto.csv: n_pre / n_post / n_removed / % removed per class.
# Shows at a glance which classes lose the most pixels.
df_confronto <- df_summary %>%
  select(COD_L2, n_pre, n_post) %>%
  mutate(
    n_rimossi   = n_pre - n_post,
    pct_rimossi = round(100 * (n_pre - n_post) / n_pre, 1)
  )
write.csv(df_confronto,
          file.path(DIR_OUT, "pixel_per_classe_confronto.csv"),
          row.names = FALSE)
cat("   pixel_per_classe_confronto.csv: written\n")

# Final summary ----
cat("\nExtraction + z-statistics summary:\n")
cat(sprintf("   Total pixels extracted:        %s\n", format(n_pre_fis,      big.mark=".")))
cat(sprintf("   After physical filter B02<=0.40: %s\n", format(n_post_fis,   big.mark=".")))
cat(sprintf("   After drop_na (101 pred.):     %s\n", format(n_post_na,      big.mark=".")))
cat(sprintf("   After z-statistics (output):   %s\n", format(nrow(df_clean), big.mark=".")))
cat(sprintf("   z-statistics predictors:       %d  (101 = 105 - 4 GLCM_COR)\n", N_PRED))
cat(sprintf("   L2 classes in the clean training: %d\n", length(unique(df_clean$COD_L2))))
cat("   Threshold: adaptive per-class mu +/- sigma\n")
cat("   GLCM_COR: excluded from z-statistics and the RF; kept in the CSV for completeness\n")

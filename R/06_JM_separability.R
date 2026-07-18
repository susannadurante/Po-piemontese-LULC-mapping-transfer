# Jeffries-Matusita separability on the final legend: L2 (18 classes) and
# L1 (7 macrotypes). Full 40x40 covariance (10 bands x 4 seasons), after
# Richards & Jia (2006).
#
# Methodological choices:
#   1. JM is computed on the final-legend GROUPS (COD_L2 and COD_L1), not on the
#      raw CLC codes: it measures separability between the classes actually mapped.
#   2. The 40 features are extracted in a SINGLE PASS from a 40-band stack. A
#      single pass keeps pixels aligned across seasons (no merge by PX_ID),
#      avoiding corrupted inter-seasonal covariances.
#   3. Reads usosuolo_2019_COD_L2.gpkg produced by the legend rasterisation; run
#      it AFTER the legend step.
#
# Reference: Richards, J.A. & Jia, X. (2006). Remote Sensing Digital Image
# Analysis, Springer - JM with the full covariance matrix.
# RAM: the extraction loads ~1.4 M pixels x 40 features (~0.5 GB); run in a clean
# R session.

# Packages ----
library(terra)

# Folders ----
# Set the project root directory here, or via the POPARK_DATA environment variable
DIR_BASE <- Sys.getenv("POPARK_DATA", unset = "path/to/project/root")
if (!dir.exists(DIR_BASE))
  stop("Invalid DIR_BASE: ", DIR_BASE,
       "
     Set the path at the top of the script, or the POPARK_DATA variable.")
DIR_OUTPUT <- file.path(DIR_BASE, "output_spettrale")
DIR_STEP3  <- file.path(DIR_OUTPUT, "STEP3_output")

# final legend = gpkg from the legend step (COD_L2 and COD_L1 per polygon)
PATH_GPKG <- file.path(DIR_STEP3, "usosuolo_2019_COD_L2.gpkg")

BANDE    <- c("B02","B03","B04","B05","B06","B07","B08","B8A","B11","B12")
STAGIONI <- c("MAM","JJA","SON","DJF")
compositi <- setNames(
  file.path(DIR_OUTPUT,
    sprintf("S2_MSIL2A_%s2019_T32TLQ-32TLR-32TMQ-32TMR_mediana_10bande.tif", STAGIONI)),
  STAGIONI)

# 40 features ordered by season, all bands (B02_MAM..B12_MAM, B02_JJA..)
FEATURE_NAMES <- unlist(lapply(STAGIONI, function(s) paste0(BANDE, "_", s)))

terraOptions(memfrac = 0.6)
options(stringsAsFactors = FALSE)

cat("JM separability (L2 + L1) on the final legend\n\n")

for (f in c(PATH_GPKG, compositi))
  if (!file.exists(f)) stop("File not found:\n  ", f)
cat("Input files verified.\n\n")

# descriptions (from the verified legend, only for the output labels; Italian)
desc_l2 <- c(
  "11"="Acque interne","12"="Lanche e acque ferme fluviali",
  "13"="Vegetazione igrofila erbacea","21"="Greto e sabbie fluviali",
  "22"="Formazioni pioniere su substrato fluviale",
  "31"="Vegetazione erbacea alta e arbustiva ripariale",
  "41"="Boschi ripariali igrofili","42"="Boschi planiziali di latifoglie mesofite",
  "43"="Robinieti","44"="Pioppeti","45"="Piantagioni arboree da legno",
  "46"="Boschi misti collinari","51"="Risaie",
  "52"="Seminativi e colture erbacee annuali","53"="Colture arboree permanenti",
  "54"="Prati stabili","61"="Tessuto urbano e aree industriali",
  "71"="Aree estrattive e cave")
desc_l1 <- c(
  "1"="Acque e zone umide","2"="Substrato fluviale",
  "3"="Vegetazione erbacea e arbustiva","4"="Copertura arborea",
  "5"="Superfici coltivate","6"="Superfici impermeabilizzate","7"="Aree estrattive")

# 40-band stack (10 bands x 4 seasons) ----
cat("Building the 40-band stack...\n")
rasters <- lapply(compositi[STAGIONI], rast)

# check that the 4 composites share the same grid (required for stacking)
for (s in STAGIONI[-1]) {
  if (nlyr(rasters[[s]]) != length(BANDE))
    stop("Composite ", s, " has ", nlyr(rasters[[s]]), " bands (expected 10)")
  if (!compareGeom(rasters[[STAGIONI[1]]], rasters[[s]], stopOnError = FALSE))
    stop("Composite ", s, " does not share the grid of ", STAGIONI[1])
}
for (s in STAGIONI) names(rasters[[s]]) <- paste0(BANDE, "_", s)

# explicit c(): do.call(c, list) does not always dispatch terra's S4 method
# (would return a list instead of a SpatRaster). Order MAM,JJA,SON,DJF.
comp40 <- c(rasters[["MAM"]], rasters[["JJA"]], rasters[["SON"]], rasters[["DJF"]])
if (nlyr(comp40) != 40) stop("Stack != 40 bands: ", nlyr(comp40))
if (!identical(names(comp40), FEATURE_NAMES))
  stop("Stack band order differs from FEATURE_NAMES")
cat("   Stack:", nlyr(comp40), "bands, order verified\n\n")

# Extract the 40 features in a single pass ----
# each row = one pixel with all 40 features already aligned (no merge by ID)
cat("Reading legend gpkg + extracting pixels...\n")
v <- vect(PATH_GPKG)
cat("   Polygons:", nrow(v), "\n")
cat("   Extracting 40 bands (may take a few minutes)...\n")

df <- terra::extract(comp40, v, ID = TRUE, touches = FALSE)

# map polygon ID -> COD_L2 / COD_L1 from the gpkg
attr_v     <- as.data.frame(v)
df$COD_L2  <- attr_v$COD_L2[df$ID]
df$COD_L1  <- attr_v$COD_L1[df$ID]
df$ID      <- NULL

# reflectance 0-1
df[FEATURE_NAMES] <- df[FEATURE_NAMES] / 10000

# keep only pixels with ALL 40 features valid (NA removal in one shot)
ok <- stats::complete.cases(df[, FEATURE_NAMES]) & !is.na(df$COD_L2)
df <- df[ok, ]
gc(verbose = FALSE)
cat("   Valid pixels (40 complete features):", format(nrow(df), big.mark="."), "\n\n")

# Functions ----
# JM with full covariance (Richards & Jia 2006):
#   B  = (1/8)(mu1-mu2)^T Sm^-1 (mu1-mu2) + (1/2)[ ln|Sm| - (1/2)(ln|S1|+ln|S2|) ]
#   JM = 2 (1 - exp(-B)),  range 0 (identical) ... 2 (well separated)
# Note: standard JM (Richards & Jia) is 2(1-exp(-B)), saturating at 2. The form
#   sqrt(2(1-exp(-B))) is NOT correct: it saturates at sqrt(2)~=1.414 and would
#   make the 1.8/1.9 thresholds unreachable.
calcola_jm <- function(mu1, mu2, s1, s2) {
  sm  <- (s1 + s2) / 2
  inv <- tryCatch(solve(sm), error = function(e) NULL)
  if (is.null(inv)) return(NA_real_)
  d   <- mu1 - mu2
  t1  <- (1/8) * as.numeric(t(d) %*% inv %*% d)
  ld_m <- as.numeric(determinant(sm, logarithm = TRUE)$modulus)
  ld_1 <- as.numeric(determinant(s1, logarithm = TRUE)$modulus)
  ld_2 <- as.numeric(determinant(s2, logarithm = TRUE)$modulus)
  t2  <- 0.5 * (ld_m - 0.5 * (ld_1 + ld_2))
  JM  <- 2 * (1 - exp(-(t1 + t2)))
  max(0, min(2, as.numeric(JM)))
}

# mean (40) + covariance (40x40) for a group, with a ridge if near-singular
mu_sigma <- function(M) {
  mu <- colMeans(M)
  sg <- stats::cov(M)
  if (tryCatch(kappa(sg, exact = FALSE) > 1e10, error = function(e) TRUE))
    sg <- sg + diag(1e-6, ncol(sg))
  list(mu = mu, sigma = sg, n = nrow(M))
}

# JM matrix + pair table for one level (grouping column)
jm_livello <- function(df, col, desc, etichetta) {
  cat(sprintf("-- %s - means/covariances per group...\n", etichetta))
  gruppi <- sort(unique(df[[col]]))
  ms <- list()
  for (g in gruppi) {
    M <- as.matrix(df[df[[col]] == g, FEATURE_NAMES])
    if (nrow(M) < ncol(M) + 1)
      cat(sprintf("   ! group %s: %d pixels (< 41) - unstable covariance\n",
                  g, nrow(M)))
    ms[[as.character(g)]] <- mu_sigma(M)
  }
  ng    <- length(gruppi)
  mat   <- matrix(NA_real_, ng, ng,
                  dimnames = list(as.character(gruppi), as.character(gruppi)))
  n_cop <- ng * (ng - 1) / 2
  A <- B <- integer(n_cop); JMv <- numeric(n_cop); nA <- nB <- integer(n_cop)
  k <- 1L
  for (i in 1:(ng - 1)) for (j in (i + 1):ng) {
    gi <- as.character(gruppi[i]); gj <- as.character(gruppi[j])
    jm <- calcola_jm(ms[[gi]]$mu, ms[[gj]]$mu, ms[[gi]]$sigma, ms[[gj]]$sigma)
    mat[gi, gj] <- jm; mat[gj, gi] <- jm
    A[k] <- gruppi[i]; B[k] <- gruppi[j]; JMv[k] <- jm
    nA[k] <- ms[[gi]]$n; nB[k] <- ms[[gj]]$n
    k <- k + 1L
  }
  # diagonal left NA: the separability of a class from itself is not a
  # meaningful pair (useful for the thesis heatmap)
  coppie <- data.frame(
    COD_A = A, DESC_A = unname(desc[as.character(A)]),
    COD_B = B, DESC_B = unname(desc[as.character(B)]),
    JM = round(JMv, 4), n_A = nA, n_B = nB)
  coppie <- coppie[order(coppie$JM), ]
  list(mat = mat, coppie = coppie)
}

# Compute L2 (18 classes) and L1 (7 macrotypes) ----
res_l2 <- jm_livello(df, "COD_L2", desc_l2, "L2 (18 classes)")
res_l1 <- jm_livello(df, "COD_L1", desc_l1, "L1 (7 macrotypes)")

# Save ----
write.csv2(res_l2$coppie, file.path(DIR_STEP3, "JM_L2_coppie.csv"),  row.names = FALSE)
write.csv2(as.data.frame(res_l2$mat), file.path(DIR_STEP3, "JM_L2_matrice.csv"), row.names = TRUE)
write.csv2(res_l1$coppie, file.path(DIR_STEP3, "JM_L1_coppie.csv"),  row.names = FALSE)
write.csv2(as.data.frame(res_l1$mat), file.path(DIR_STEP3, "JM_L1_matrice.csv"), row.names = TRUE)
cat("\nSaved to STEP3_output/:\n")
cat("   JM_L2_coppie.csv | JM_L2_matrice.csv (18x18)\n")
cat("   JM_L1_coppie.csv | JM_L1_matrice.csv (7x7)\n\n")

# Console summary ----
cat("L2 - 10 least separable pairs (lowest JM):\n")
print(head(res_l2$coppie, 10), row.names = FALSE)
cat("\nL1 - least separable pairs:\n")
print(head(res_l1$coppie, 10), row.names = FALSE)
cat("\nIndicative thresholds: JM<1.0 nearly inseparable | 1.0-1.8 borderline | >1.9 well separated\n")

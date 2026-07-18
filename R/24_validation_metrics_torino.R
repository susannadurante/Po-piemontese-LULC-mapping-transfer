# Full validation metrics (Turin sector) + comparison with the Vercelli k-fold.
# The RF training script computes Richiardi's full set on Vercelli (OA, Kappa, UA,
# PA, SPEC, F1, BA, WA, TSS, G-Mean). To place "Vercelli k-fold" and "Turin
# expert-judgment" side by side in the thesis, the columns must match: this reuses
# calcola_metriche taken VERBATIM from the training script, so both sides of the
# comparison come from the same code.
#
# Matrix orientation (note): in the training script the matrix is rows = PREDICTED,
# cols = TRUE. Earlier validation matrices used the opposite. Here the training
# orientation (pred, true) is adopted for consistency. PA and UA keep their usual
# definitions: PA = TP/true_k, UA = TP/predicted_k. State the orientation in the
# thesis caption if the tables are placed side by side.
#
# Olofsson does NOT apply: the area-adjusted estimators need a probability sample.
# The Turin sample is expert-judgment: no unbiased estimates, no intervals. The
# metrics here are descriptive of the sample and likely OPTIMISTIC (expert judgment
# tends to pick pure examples, which the model gets right).

# Packages ----
suppressPackageStartupMessages({ library(terra); library(sf) })

# Folders ----
# Set the project root directory here, or via the POPARK_DATA environment variable
DIR_BASE <- Sys.getenv("POPARK_DATA", unset = "path/to/project/root")
if (!dir.exists(DIR_BASE))
  stop("Invalid DIR_BASE: ", DIR_BASE,
       "
     Set the path at the top of the script, or the POPARK_DATA variable.")
DIR_VAL   <- file.path(DIR_BASE, "validazione")
DIR_CLASS <- file.path(DIR_BASE, "output_spettrale_intero_parco", "08_classification_output")
DIR_RF    <- file.path(DIR_BASE, "output_spettrale", "07_RF_output")

PATH_PT_L1  <- file.path(DIR_VAL, "validazione_punti_L1_torino.gpkg")
PATH_PT_L2  <- file.path(DIR_VAL, "validazione_punti_L2_torino.gpkg")
PATH_RAS_L1 <- file.path(DIR_CLASS, "classificazione_L1_clip.tif")
PATH_RAS_L2 <- file.path(DIR_CLASS, "classificazione_L2_clip.tif")

CL_L1 <- 1:7
CL_L2 <- c(11,12,13, 21,22, 31, 41,42,43,44,45,46, 51,52,53,54, 61, 71)

# aggregation scheme B2 (within the macrotype)
AGG <- c("11"=11L, "12"=12L, "13"=12L, "21"=20L, "22"=20L, "31"=31L,
         "41"=41L, "42"=42L, "43"=43L, "46"=46L, "44"=40L, "45"=40L,
         "51"=51L, "52"=52L, "53"=53L, "54"=54L, "61"=61L, "71"=71L)
CL_L2A <- c(11, 12, 20, 31, 41, 42, 43, 46, 40, 51, 52, 53, 54, 61, 71)

# Metric function (copied verbatim from the RF training script) ----
# Matrix oriented rows = PREDICTED, cols = TRUE. Verified there against
# caret::confusionMatrix to 1e-9.
calcola_metriche <- function(pred, true, all_levels) {
  cm     <- table(factor(pred, levels=all_levels),
                  factor(true, levels=all_levels))
  N      <- sum(cm); diagv <- diag(cm)
  pred_k <- rowSums(cm); true_k <- colSums(cm)

  p_o   <- sum(diagv) / N
  p_e   <- sum(pred_k * true_k) / N^2
  Kappa <- (p_o - p_e) / (1 - p_e)

  UA   <- diagv / pmax(pred_k, 1)
  PA   <- diagv / pmax(true_k, 1)
  FP   <- pred_k - diagv; FN <- true_k - diagv
  TN   <- N - diagv - FP - FN
  SPEC <- TN / pmax(TN + FP, 1)
  F1   <- ifelse((UA + PA) > 0, 2 * UA * PA / (UA + PA), 0)
  BA_k  <- (PA + SPEC) / 2
  TSS_k <- PA + SPEC - 1

  BA  <- mean(BA_k,  na.rm = TRUE)
  TSS <- mean(TSS_k, na.rm = TRUE)
  WA  <- sum((true_k / N) * UA, na.rm = TRUE)
  GM  <- prod(PA[true_k > 0])^(1 / sum(true_k > 0))

  list(cm_table = cm,
       OA    = round(100 * p_o, 2),
       Kappa = round(Kappa, 4),
       BA    = round(100 * BA, 2),
       WA    = round(100 * WA, 2),
       TSS   = round(TSS, 4),
       GMean = round(GM, 4),
       UA    = round(100 * UA,   2),
       PA    = round(100 * PA,   2),
       SPEC  = round(100 * SPEC, 2),
       F1    = round(100 * F1,   2),
       BA_k  = round(100 * BA_k, 2),
       TSS_k = round(TSS_k, 4),
       n_true = as.integer(true_k),
       n_pred = as.integer(pred_k))
}

# Extraction ----
leggi <- function(path_pt, campo, path_ras, etichetta) {
  if (!file.exists(path_pt))  stop("Points not found (", etichetta, "): ", path_pt)
  if (!file.exists(path_ras)) {
    cat("[!] Raster not found: ", path_ras, "\n")
    if (dir.exists(dirname(path_ras)))
      for (f in list.files(dirname(path_ras), pattern="\\.tif$")) cat("      - ", f, "\n")
    stop("Fix the raster path ", etichetta)
  }
  pt  <- st_read(path_pt, quiet = TRUE)
  ras <- rast(path_ras)
  if (st_crs(pt)$wkt != crs(ras)) pt <- st_transform(pt, crs(ras))
  list(ref  = as.integer(pt[[campo]]),
       pred = terra::extract(ras, vect(pt))[[2]])
}

cat("Full metrics - Turin validation (expert-judgment)\n\n")

d1 <- leggi(PATH_PT_L1, "classe_L1", PATH_RAS_L1, "L1")
d2 <- leggi(PATH_PT_L2, "classe_L2", PATH_RAS_L2, "L2")

# Three evaluations ----
valuta <- function(ref, pred, livelli, nome) {
  ok <- !is.na(ref) & !is.na(pred) & ref %in% livelli & pred %in% livelli
  r  <- calcola_metriche(pred[ok], ref[ok], livelli)
  cat(sprintf("\n== %s  (n = %d) ==\n", nome, sum(ok)))
  cat(sprintf("   OA %.2f%%  |  BA %.2f%%  |  WA %.2f%%  |  TSS %.4f  |  G-Mean %.4f  |  Kappa %.4f\n",
      r$OA, r$BA, r$WA, r$TSS, r$GMean, r$Kappa))
  cat(sprintf("\n   %-6s %7s %7s %8s %8s %8s %8s %8s\n",
      "classe","n_vero","n_pred","PA%","UA%","SPEC%","F1%","TSS"))
  for (i in seq_along(livelli))
    cat(sprintf("   %-6s %7d %7d %8.1f %8.1f %8.1f %8.1f %8.3f\n",
        livelli[i], r$n_true[i], r$n_pred[i],
        r$PA[i], r$UA[i], r$SPEC[i], r$F1[i], r$TSS_k[i]))
  r
}

r_L1  <- valuta(d1$ref, d1$pred, CL_L1, "L1 - 7 macrotypes")
r_L2  <- valuta(d2$ref, d2$pred, CL_L2, "L2 full - 18 classes")

ric <- function(x) { y <- AGG[as.character(x)]
  if (any(!is.na(x) & is.na(y))) stop("Codes outside the aggregation scheme")
  as.integer(y) }
r_L2A <- valuta(ric(d2$ref), ric(d2$pred), CL_L2A, "L2 aggregated (B2) - 15 classes")

# Comparison with the Vercelli k-fold ----
# read from the kfold_summary.csv files produced by the training script: no numbers
# retyped by hand, no risk of reporting stale values in the thesis.
cat("\n\n== Comparison: Vercelli (k-fold, internal) vs Turin (extrapolation) ==\n\n")

leggi_kfold <- function(m) {
  f <- file.path(DIR_RF, m, "kfold_summary.csv")
  if (!file.exists(f)) return(NULL)
  read.csv(f)
}

cat(sprintf("   %-8s %10s %10s %10s %10s\n", "model", "OA%", "BA%", "WA%", "Kappa"))
cat("   -- Vercelli (k-fold pooled) --\n")
trovati <- 0L
for (m in c("L1","L2_1","L2_2","L2_4","L2_5")) {
  k <- leggi_kfold(m)
  if (is.null(k)) { cat(sprintf("   %-8s   [kfold_summary.csv not found]\n", m)); next }
  trovati <- trovati + 1L
  g <- function(nm) if (nm %in% names(k)) k[[nm]][1] else NA
  cat(sprintf("   %-8s %10s %10s %10s %10s\n", m,
      format(g("OA")), format(g("BA")), format(g("WA")), format(g("Kappa"))))
}
if (trovati == 0L)
  cat("   [!] No kfold_summary.csv found in ", DIR_RF, "\n", sep = "")

cat("   -- Turin (expert-judgment) --\n")
cat(sprintf("   %-8s %10.2f %10.2f %10.2f %10.4f\n", "L1",       r_L1$OA,  r_L1$BA,  r_L1$WA,  r_L1$Kappa))
cat(sprintf("   %-8s %10.2f %10.2f %10.2f %10.4f\n", "L2 full",  r_L2$OA,  r_L2$BA,  r_L2$WA,  r_L2$Kappa))
cat(sprintf("   %-8s %10.2f %10.2f %10.2f %10.4f\n", "L2 aggr",  r_L2A$OA, r_L2A$BA, r_L2A$WA, r_L2A$Kappa))

cat("\n   Reading:\n")
cat("   The k-fold estimates the MODEL on the labelled area (interpolation).\n")
cat("   The Turin points estimate the MAP under extrapolation. The gap between the\n")
cat("   two is the transfer cost, and it is the central result of the work.\n")
cat("   The Turin metrics are descriptive and likely optimistic: the sample is\n")
cat("   expert-judgment, not probabilistic -> no Olofsson, no confidence intervals.\n")

# Save ----
salva <- function(r, livelli, tag) {
  write.csv(as.data.frame.matrix(r$cm_table),
            file.path(DIR_VAL, sprintf("matrice_%s_predXvero.csv", tag)))
  write.csv(data.frame(
    classe = livelli, n_vero = r$n_true, n_pred = r$n_pred,
    PA = r$PA, UA = r$UA, SPEC = r$SPEC, F1 = r$F1, BA = r$BA_k, TSS = r$TSS_k
  ), file.path(DIR_VAL, sprintf("metriche_%s_complete.csv", tag)), row.names = FALSE)
}
salva(r_L1,  CL_L1,  "L1")
salva(r_L2,  CL_L2,  "L2")
salva(r_L2A, CL_L2A, "L2_aggregata")

glob <- data.frame(
  livello = c("L1", "L2_completo", "L2_aggregato"),
  n     = c(sum(r_L1$cm_table), sum(r_L2$cm_table), sum(r_L2A$cm_table)),
  OA    = c(r_L1$OA, r_L2$OA, r_L2A$OA),
  BA    = c(r_L1$BA, r_L2$BA, r_L2A$BA),
  WA    = c(r_L1$WA, r_L2$WA, r_L2A$WA),
  TSS   = c(r_L1$TSS, r_L2$TSS, r_L2A$TSS),
  GMean = c(r_L1$GMean, r_L2$GMean, r_L2A$GMean),
  Kappa = c(r_L1$Kappa, r_L2$Kappa, r_L2A$Kappa)
)
write.csv(glob, file.path(DIR_VAL, "metriche_globali_torino.csv"), row.names = FALSE)

cat("\nSaved in validazione/:\n")
cat("     matrice_L1_predXvero.csv, matrice_L2_predXvero.csv, matrice_L2_aggregata_predXvero.csv\n")
cat("     metriche_L1_complete.csv, metriche_L2_complete.csv, metriche_L2_aggregata_complete.csv\n")
cat("     metriche_globali_torino.csv\n")

# Random Forest: tuning + k-fold validation + final model.
# For one model of the hierarchy (L1 or an L2_k), in a single pass:
#   1. hyperparameter TUNING (mtry, min.node.size) on a subsample;
#   2. internal VALIDATION by class-stratified k-fold cross-validation at pixel
#      level -> estimates the model accuracy;
#   3. FINAL MODEL trained and saved, for the classification step.
#
# PARAMETRIC: this single script trains and validates ONE model per run, selected
#   by the MODELLO variable below. Run it 5 times, setting MODELLO to L1, L2_1,
#   L2_2, L2_4, L2_5 in turn, to reproduce all five models.
#
# How the k-fold works (per pixel): recombines train_pixel.csv + validation_pixel.csv
#   = all the model's pixels, then reassigns k stratified folds. Each fold is the
#   hold-out in turn while the others train: every pixel is predicted EXACTLY ONCE.
#   Two views follow: a pooled confusion matrix over all pixels (basis for
#   Olofsson), and per-fold OA -> mean and variability band across folds.
#
# What it estimates (and not): this is INTERNAL validation on the labelled area; it
#   estimates the quality of the MODEL. It does NOT validate the Turin-sector map,
#   which is extrapolation: that needs the independent photo-interpreted set and the
#   Area of Applicability map (CAST::aoa), downstream of the classification.
#
# Predictor selection (optional, automatic): if predittori_selezionati_<MODELLO>.csv
#   (from the selection step) exists, the model uses only those predictors;
#   otherwise all of them. GLCM_COR (structural NAs) are always excluded, since
#   ranger does not handle missing values.
# Balancing: baseline = standard RF, no weights. Optional hook USA_PESI=TRUE ->
#   case.weights = 1/n computed ONLY on each fold's training set.
#
# Run one model at a time in a fresh R session (16 GB RAM), from the terminal.

# Packages ----
library(data.table)
library(ranger)
# NO caret: folds and metrics computed by hand (data.table) to save RAM

# Model selection ----
MODELLO <- "L1"   # the model for this run (set to L1, L2_1, L2_2, L2_4 or L2_5)

# per-model configuration: target + optional floor on min.node.size
CONFIG <- list(
  L1   = list(target = "COD_L1", node_floor = 1L),
  L2_1 = list(target = "COD_L2", node_floor = 1L),
  L2_2 = list(target = "COD_L2", node_floor = 1L),
  L2_4 = list(target = "COD_L2", node_floor = 1L),
  L2_5 = list(target = "COD_L2", node_floor = 5L)  # L2=53 ratio 71:1 -> avoid Empty node
)
stopifnot(MODELLO %in% names(CONFIG))
TARGET     <- CONFIG[[MODELLO]]$target
NODE_FLOOR <- CONFIG[[MODELLO]]$node_floor

# Folders and parameters ----
# Set the project root directory here, or via the POPARK_DATA environment variable
DIR_BASE <- Sys.getenv("POPARK_DATA", unset = "path/to/project/root")
if (!dir.exists(DIR_BASE))
  stop("Invalid DIR_BASE: ", DIR_BASE,
       "
     Set the path at the top of the script, or the POPARK_DATA variable.")
DIR_IN   <- file.path(DIR_BASE, "output_spettrale", "05_06_output", MODELLO)
DIR_OUT  <- file.path(DIR_BASE, "output_spettrale", "07_RF_output", MODELLO)
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

SEED        <- 50L
NUM_CORES   <- 6L     # 6 threads: RAM is free, this is the main speed lever (the
                       # crashes came from min.node.size=1, not from the cores).
                       # Drop to 4L on a machine with few cores.
NUM_TREES   <- 200L
K_FOLD      <- 5L      # number of folds (5 = RAM / band-stability compromise)
PROP_TUNING <- 0.03    # 3% stratified for the tuning (once only)
PROP_TRAIN  <- 0.80    # fraction of pixels used to train (each fold and the final
                       # model). 0.80 = as much data as possible while keeping RAM
                       # in check; maximises accuracy, especially on minority
                       # classes. Lower it on smaller machines.
MTRY_MAX    <- 34L     # cap on the tuned mtry, for fold-training speed. The grid
                       # explores up to ~p/3; with ~101 predictors the max is ~34,
                       # so this safety cap rarely bites. Set NULL to disable.
USA_PESI    <- FALSE   # TRUE = case.weights 1/n on the fold's training only

options(ranger.num.threads = NUM_CORES)
set.seed(SEED)

COLS_ESCLUSE <- c("ID_PIXEL", "x", "y", "COD_L1", "COD_L2")

cat(sprintf("RF k-fold | model=%s | target=%s | k=%d | seed=%d | weights=%s\n",
    MODELLO, TARGET, K_FOLD, SEED, USA_PESI))

# Olofsson 2014 (area-adjusted accuracy) ----
# Reference: Olofsson P. et al. (2014) Remote Sens. Environ. 148:42-57.
# Note: with the k-fold the validation comes from the folds (internal CV), not from
# an independent probability sample -> read the estimates as relative indicators.
calcola_olofsson <- function(cm_table, pixel_area_m2 = 100) {
  nclass <- nrow(cm_table)
  conf   <- 1.96
  maparea <- rowSums(cm_table) * pixel_area_m2 / 10000
  A    <- sum(maparea)
  W_i  <- maparea / A
  n_i  <- rowSums(cm_table)
  p    <- sweep(cm_table, 1, n_i, "/")
  p    <- sweep(p, 1, W_i, "*")
  p[is.nan(p)] <- 0
  p_area    <- colSums(p) * A
  p_area_CI <- conf * A * sqrt(colSums((W_i * p - p^2) / pmax(n_i-1, 1)))
  OA    <- sum(diag(p))
  PA    <- diag(p) / colSums(p)
  UA    <- diag(p) / rowSums(p)
  F1    <- ifelse((PA+UA)>0, 2*PA*UA/(PA+UA), 0)
  OA_CI <- conf * sqrt(sum(W_i^2 * UA * (1-UA) / pmax(n_i-1, 1)))
  UA_CI <- conf * sqrt(UA * (1-UA) / pmax(n_i-1, 1))
  N_j   <- sapply(seq_len(nclass), function(x) sum(maparea/n_i*cm_table[,x], na.rm=TRUE))
  tmp   <- sapply(seq_len(nclass), function(x)
    sum(maparea[-x]^2 * cm_table[-x,x] / n_i[-x] *
        (1 - cm_table[-x,x]/n_i[-x]) / pmax(n_i[-x]-1, 1), na.rm=TRUE))
  PA_CI <- conf * sqrt(1/N_j^2 *
    (maparea^2*(1-PA)^2*UA*(1-UA)/pmax(n_i-1,1) + PA^2*tmp))
  data.frame(
    classe=rownames(cm_table),
    area_ha=round(p_area,2), area_ha_CI=round(p_area_CI,2),
    PA=round(PA*100,2), PA_CI=round(PA_CI*100,2),
    UA=round(UA*100,2), UA_CI=round(UA_CI*100,2),
    F1=round(F1*100,2), OA=round(OA*100,2), OA_CI=round(OA_CI*100,2)
  )
}

# Manual metrics (no caret) ----
# OA, Kappa, per-class UA/PA/F1/SPEC, plus aggregate indices for imbalanced data
# (BA, WA, TSS, G-Mean), computed by hand from the confusion matrix. Verified
# numerically == caret::confusionMatrix (OA, Kappa, PA=Sensitivity, UA=Pos Pred
# Value, Specificity, Balanced Accuracy, TSS) to 1e-9.
#   UA  = precision = TP/(predicted k) ; PA = recall/sensitivity = TP/(true k)
#   SPEC= TN/(TN+FP) ; F1 = 2*UA*PA/(UA+PA)
#   BA  = mean_k (PA_k+SPEC_k)/2  (balanced accuracy, macro)
#   TSS = mean_k (PA_k+SPEC_k-1)  (true skill statistic, Allouche 2006, macro)
#   WA  = sum_k (true_freq_k/N)*UA_k  (weighted accuracy)
#   GM  = (prod_k PA_k)^(1/K)  (G-mean of the recalls)
# Units: OA,UA,PA,F1,SPEC,BA,WA,BA_k in % ; Kappa,TSS,TSS_k,GMean in [0,1]/[-1,1].
# The matrix is oriented rows=PREDICTED, cols=TRUE (as in the rest of the script).
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
       TSS_k = round(TSS_k, 4))
}

# Stratified per-class fold assignment (no caret) ----
# round-robin within each class after a shuffle -> each class is spread ~evenly
# over all folds, and each fold keeps the natural (imbalanced) proportions.
assegna_fold <- function(y, K, seed) {
  set.seed(seed)
  fold <- integer(length(y))
  for (cl in unique(y)) {
    idx <- which(y == cl)
    idx <- sample(idx)                                 # shuffle within the class
    fold[idx] <- rep_len(seq_len(K), length(idx))      # round-robin
  }
  fold
}

# Load and recombine (Split A train+val = all pixels) ----
f_train <- file.path(DIR_IN, "train_pixel.csv")
f_val   <- file.path(DIR_IN, "validation_pixel.csv")
if (!file.exists(f_train) || !file.exists(f_val))
  stop("Missing train_pixel.csv / validation_pixel.csv in ", DIR_IN)

full <- rbind(fread(f_train), fread(f_val))
setDT(full)
full[[TARGET]] <- as.factor(full[[TARGET]])
all_levels <- levels(full[[TARGET]])
cat(sprintf("   Total pixels (recombined): %s | classes: %d\n",
    format(nrow(full), big.mark="."), length(all_levels)))

cols_pred <- setdiff(names(full), COLS_ESCLUSE)
n_tot_pred <- length(cols_pred)
# predictor selection: if DIR_IN has predittori_selezionati_<MODELLO>.csv (from the
# selection step), the model uses ONLY those predictors. With RICHIEDI_SELEZIONE=TRUE,
# if the file is missing the script STOPS (avoids accidentally training on all 101);
# set it to FALSE only to deliberately produce the all-predictors baseline.
RICHIEDI_SELEZIONE <- TRUE
f_sel <- file.path(DIR_IN, paste0("predittori_selezionati_", MODELLO, ".csv"))
if (file.exists(f_sel)) {
  sel <- fread(f_sel)$predittore
  cols_pred <- intersect(cols_pred, sel)
  if (length(cols_pred) == 0)
    stop("The selection file has no valid predictors for this model.")
  cat(sprintf("   Predictor selection ACTIVE: %d of %d used\n",
      length(cols_pred), n_tot_pred))
} else if (RICHIEDI_SELEZIONE) {
  stop(sprintf(paste0("%s NOT found in:\n     %s\n",
       "   -> Run the selection step first, or check DIR_IN.\n",
       "   -> To deliberately train on ALL predictors: RICHIEDI_SELEZIONE <- FALSE"),
       basename(f_sel), DIR_IN))
} else {
  cat(sprintf("   Predictors: %d (all, BASELINE without selection)\n", n_tot_pred))
}

# exclude predictors with NAs: ranger does NOT handle missing values and stops with
# "Missing data in columns: ...". GLCM_COR have STRUCTURAL NAs (sigma=0 in the local
# window on uniform surfaces -> undefined correlation); they are already excluded
# from the z-statistics but remain in the CSV. Here they are excluded from the RF
# predictors: consistent with the z-statistics, without losing pixels or imputing.
na_cols <- cols_pred[vapply(cols_pred, function(cc) anyNA(full[[cc]]), logical(1))]
if (length(na_cols) > 0) {
  cat(sprintf("   Predictors excluded for NAs (not handled by ranger): %s\n",
      paste(na_cols, collapse = ", ")))
  cols_pred <- setdiff(cols_pred, na_cols)
  cat(sprintf("   Predictors used by the RF: %d\n", length(cols_pred)))
}

# rare-class check against k:
#  - n < K  -> STOP: with round-robin some folds would be WITHOUT that class
#  - n < 2K -> warning: class present everywhere but per-fold metrics unstable
tab0 <- table(full[[TARGET]])
rotte <- names(tab0)[tab0 < K_FOLD]
if (length(rotte) > 0)
  stop(sprintf("Classes with fewer than K=%d pixels (empty folds): %s - lower K or merge.",
       K_FOLD, paste(rotte, collapse=", ")))
instabili <- names(tab0)[tab0 < K_FOLD * 2]
if (length(instabili) > 0)
  warning(sprintf("Classes with few pixels for k=%d (unstable per-fold metrics): %s",
          K_FOLD, paste(instabili, collapse=", ")))

fold_id <- assegna_fold(full[[TARGET]], K_FOLD, SEED)
cat("   Pixel distribution per fold x class:\n")
print(table(fold = fold_id, classe = full[[TARGET]]))

# Tuning (once, on a stratified subsample) ----
set.seed(SEED)
idx_tune <- full[, .(ri = .I[sample(.N, max(1L, round(.N * PROP_TUNING)))]),
                 by = TARGET]$ri
train_tune <- full[idx_tune]
cat(sprintf("   Tuning on %s pixels (%.0f%%)...\n",
    format(nrow(train_tune), big.mark="."), PROP_TUNING * 100))

# mtry grid ADAPTIVE to the number of predictors p: with high p (~90-105) a fixed
# grid up to ~sqrt(p)~=10 misses the optimum (often well beyond 10). Explore from a
# low fraction up to ~p/3.
p_pred    <- length(cols_pred)
mtry_vals <- sort(unique(pmax(1L, c(floor(sqrt(p_pred)/2), floor(sqrt(p_pred)),
                                    round(p_pred/5), round(p_pred/3)))))
mtry_vals <- mtry_vals[mtry_vals <= p_pred]
hyper_grid <- expand.grid(mtry = mtry_vals,
                          min_node_size = c(5L, 10L),   # NOT 1: with node=1 trees
                          OOB = NA_real_)               # grow maximal -> forest
                                                        # ~5-6 GB -> out-of-memory.
                                                        # node>=5 cuts memory and
                                                        # time ~5x, OA ~unchanged.
cat(sprintf("   mtry grid (p=%d): %s\n", p_pred,
    paste(mtry_vals, collapse=", ")))
set.seed(SEED)
cat(sprintf("   Starting tuning: %d combinations on %d cores (a few minutes)...\n",
    nrow(hyper_grid), NUM_CORES))
for (i in seq_len(nrow(hyper_grid))) {
  rf <- ranger(x = train_tune[, ..cols_pred], y = train_tune[[TARGET]],
               num.trees = 50L, mtry = hyper_grid$mtry[i],
               min.node.size = hyper_grid$min_node_size[i],
               num.threads = NUM_CORES, seed = SEED,
               write.forest = FALSE, verbose = FALSE)
  hyper_grid$OOB[i] <- rf$prediction.error
  cat(sprintf("      [%d/%d] mtry=%2d node=%2d  ->  OOB=%.4f\n",
      i, nrow(hyper_grid), hyper_grid$mtry[i],
      hyper_grid$min_node_size[i], hyper_grid$OOB[i]))
  flush.console()
  rm(rf); gc(verbose = FALSE)
}
hyper_grid <- hyper_grid[order(hyper_grid$OOB), ]
best_mtry <- hyper_grid$mtry[1]
if (!is.null(MTRY_MAX) && best_mtry > MTRY_MAX) {
  cat(sprintf("   optimal mtry=%d capped to MTRY_MAX=%d (speed)\n",
      best_mtry, MTRY_MAX))
  best_mtry <- MTRY_MAX
}
best_node <- max(hyper_grid$min_node_size[1], NODE_FLOOR)
cat(sprintf("   Optimal: mtry=%d, min.node.size=%d (floor=%d), OOB=%.4f\n",
    best_mtry, best_node, NODE_FLOOR, hyper_grid$OOB[1]))
fwrite(hyper_grid, file.path(DIR_OUT, "tuning.csv"))
rm(train_tune, idx_tune); gc(verbose = FALSE)

# K-fold loop ----
pred_oof    <- character(nrow(full))      # out-of-fold predictions (1 per pixel)
oa_per_fold <- numeric(K_FOLD)
ka_per_fold <- numeric(K_FOLD)

# helper: normalised case.weights 1/n, computed on the given training set
calcola_pesi <- function(dt) {
  n_cl <- dt[, .N, by = TARGET]
  w    <- setNames(1 / n_cl$N, as.character(n_cl[[TARGET]]))
  w    <- w / sum(w)
  as.numeric(w[as.character(dt[[TARGET]])])
}

for (j in seq_len(K_FOLD)) {
  idx_va <- which(fold_id == j)
  idx_tr <- which(fold_id != j)

  # stratified subsampling of the fold's training set (index-first -> avoids big copies)
  set.seed(SEED + j)
  helper <- data.table(g = full[[TARGET]][idx_tr], pos = idx_tr)
  keep <- helper[, .(ri = pos[sample(.N, max(1L, round(.N * PROP_TRAIN)))]),
                 by = g]$ri
  rm(helper)
  tr_dt <- full[keep]

  cw <- if (USA_PESI) calcola_pesi(tr_dt) else NULL

  set.seed(SEED + j)
  rf <- ranger(x = tr_dt[, ..cols_pred], y = tr_dt[[TARGET]],
               num.trees = NUM_TREES, mtry = best_mtry,
               min.node.size = best_node, case.weights = cw,
               importance = "none", num.threads = NUM_CORES,
               seed = SEED, verbose = FALSE)
  rm(tr_dt); gc(verbose = FALSE)

  va_dt <- full[idx_va]
  pr <- predict(rf, data = va_dt[, ..cols_pred],
                num.threads = NUM_CORES)$predictions
  rm(rf); gc(verbose = FALSE)

  pred_oof[idx_va] <- as.character(pr)
  m <- calcola_metriche(pr, va_dt[[TARGET]], all_levels)
  oa_per_fold[j] <- m$OA
  ka_per_fold[j] <- m$Kappa
  cat(sprintf("   Fold %d/%d: OA=%.1f%%  Kappa=%.3f  (train=%s px, val=%s px)\n",
      j, K_FOLD, m$OA, m$Kappa,
      format(length(keep), big.mark="."), format(length(idx_va), big.mark=".")))
  rm(va_dt, pr, m, keep); gc(verbose = FALSE)
}

# Aggregate metrics ----
# (a) pooled: each pixel predicted once -> single CM over all pixels
met_pool <- calcola_metriche(pred_oof, as.character(full[[TARGET]]), all_levels)
cat(sprintf("\n   POOLED (out-of-fold): OA=%.1f%%  Kappa=%.3f  BA=%.1f%%  WA=%.1f%%  TSS=%.3f  G-Mean=%.3f\n",
    met_pool$OA, met_pool$Kappa, met_pool$BA, met_pool$WA, met_pool$TSS, met_pool$GMean))

# (b) between-fold variability. Note: the k estimates are NOT independent (the
#     training sets overlap) -> the t formula underestimates the true variance and
#     is NOT a rigorous frequentist 95% CI (Bengio & Grandvalet 2004, JMLR 5:1089-
#     1105). Reported as an indicative EMPIRICAL band of between-fold variability.
oa_mean <- mean(oa_per_fold); oa_sd <- sd(oa_per_fold)
oa_band <- qt(0.975, K_FOLD - 1) * oa_sd / sqrt(K_FOLD)   # band width (t, k-1 df)
cat(sprintf("   CV OA: mean=%.1f%%  sd=%.2f  fold-var. band (indicative)=[%.1f, %.1f]\n",
    oa_mean, oa_sd, oa_mean - oa_band, oa_mean + oa_band))

fwrite(data.frame(fold = seq_len(K_FOLD), OA = oa_per_fold, Kappa = ka_per_fold),
       file.path(DIR_OUT, "kfold_per_fold.csv"))
fwrite(data.frame(OA_media = round(oa_mean,2), OA_sd = round(oa_sd,2),
                  banda_inf = round(oa_mean - oa_band,2),
                  banda_sup = round(oa_mean + oa_band,2),
                  OA_pooled = met_pool$OA, Kappa_pooled = met_pool$Kappa,
                  BA_pooled = met_pool$BA, WA_pooled = met_pool$WA,
                  TSS_pooled = met_pool$TSS, GMean_pooled = met_pool$GMean),
       file.path(DIR_OUT, "kfold_summary.csv"))
fwrite(data.frame(classe = all_levels, UA = met_pool$UA, PA = met_pool$PA,
                  F1 = met_pool$F1, SPEC = met_pool$SPEC,
                  BA = met_pool$BA_k, TSS = met_pool$TSS_k),
       file.path(DIR_OUT, "validation_kfold.csv"))
fwrite(as.data.frame(met_pool$cm_table),
       file.path(DIR_OUT, "confusion_matrix_kfold.csv"))

olf <- tryCatch(calcola_olofsson(as.matrix(met_pool$cm_table)),
                error = function(e) { warning(paste("Olofsson:", e$message)); NULL })
if (!is.null(olf)) {
  fwrite(olf, file.path(DIR_OUT, "olofsson_kfold.csv"))
  cat(sprintf("   Olofsson area-adjusted OA=%.1f%% +/- %.1f%%\n", olf$OA[1], olf$OA_CI[1]))
}
rm(pred_oof); gc(verbose = FALSE)

# Final model for classification ----
# the k-fold ESTIMATES the accuracy; the map is produced with this model. It is
# trained on the SAME PROP_TRAIN fraction (0.80) as the folds, so the k-fold
# estimate is FAITHFUL to this model. If the folds used less data than the final
# model, the CV would underestimate it; using the same fraction, it does not. If
# you raise PROP_TRAIN only for the final model, state that the k-fold estimate is
# conservative.
cat("\nFinal model (80% of the pixels) for the classification...\n")
set.seed(SEED)
idx_fin <- full[, .(ri = .I[sample(.N, max(1L, round(.N * PROP_TRAIN)))]),
                by = TARGET]$ri
fin_dt <- full[idx_fin]
cw_fin <- if (USA_PESI) calcola_pesi(fin_dt) else NULL
set.seed(SEED)
rf_final <- ranger(x = fin_dt[, ..cols_pred], y = fin_dt[[TARGET]],
                   num.trees = NUM_TREES, mtry = best_mtry,
                   min.node.size = best_node, case.weights = cw_fin,
                   importance = "none", num.threads = NUM_CORES,
                   seed = SEED, verbose = FALSE)
saveRDS(rf_final, file.path(DIR_OUT, "rf_model_final.rds"))
cat(sprintf("   Final model OOB: %.4f | saved: rf_model_final.rds\n",
    rf_final$prediction.error))
rm(fin_dt, idx_fin, rf_final); gc(verbose = FALSE)

# Variable importance (separate model, on a subsample for RAM) ----
# this was the heaviest step (it ran on the whole dataset). On a stratified
# subsample the importance ranking stays stable and RAM holds.
cat("\nVariable importance (subsample)...\n")
set.seed(SEED)
idx_vi <- full[, .(ri = .I[sample(.N, max(1L, round(.N * PROP_TRAIN)))]),
               by = TARGET]$ri
vi_dt  <- full[idx_vi]
set.seed(SEED)
rf_vi <- ranger(x = vi_dt[, ..cols_pred], y = vi_dt[[TARGET]],
                num.trees = 200L, mtry = best_mtry, importance = "impurity",
                num.threads = NUM_CORES, seed = SEED, verbose = FALSE)
vi <- sort(rf_vi$variable.importance, decreasing = TRUE)
fwrite(data.frame(predittore = names(vi), importanza = round(vi, 2),
                  rank = seq_along(vi)),
       file.path(DIR_OUT, "variable_importance.csv"))
cat(sprintf("   Top 5: %s\n", paste(names(vi)[1:5], collapse = ", ")))
rm(rf_vi, vi_dt, idx_vi); gc(verbose = FALSE)

# Summary ----
cat("\nRF k-fold summary -", MODELLO, "\n")
cat(sprintf("   CV OA: %.1f%% +/- %.1f (fold-var. band, indicative)  |  Pooled OA: %.1f%%  Kappa: %.3f\n",
    oa_mean, oa_band, met_pool$OA, met_pool$Kappa))
cat(sprintf("   Pooled BA: %.1f%%  WA: %.1f%%  TSS: %.3f  G-Mean: %.3f\n",
    met_pool$BA, met_pool$WA, met_pool$TSS, met_pool$GMean))
cat(sprintf("   Output in: %s\n", DIR_OUT))
cat("   Files: kfold_per_fold.csv, kfold_summary.csv, validation_kfold.csv,\n")
cat("          confusion_matrix_kfold.csv, olofsson_kfold.csv, tuning.csv,\n")
cat("          variable_importance.csv, rf_model_final.rds\n")
cat("   Note: INTERNAL estimate (labelled area). The Turin sector is extrapolation\n")
cat("         -> validate with the independent photo-interpreted set + Area of\n")
cat("         Applicability map (CAST::aoa) downstream of the classification.\n")

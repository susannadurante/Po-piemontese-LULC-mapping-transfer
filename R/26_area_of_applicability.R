# Area of Applicability (AoA) - Meyer & Pebesma (2021), package CAST.
# A Random Forest ALWAYS predicts a class, even for pixels that resemble nothing in
# the training set. The AoA quantifies WHERE the prediction is reliable:
#   DI (Dissimilarity Index): distance of each pixel from the nearest training point
#     in predictor space, WEIGHTED by the importance of each predictor.
#   threshold: maximum DI observed in the training set under cross-validation.
#   AOA: 1 where DI <= threshold (inside the domain), 0 where DI > threshold
#     (extrapolation).
#
# Hypothesis under test: the Turin floodplain falls OUTSIDE the AoA, because its
# 'dem' values (~192 m) are new with respect to the Vercelli training (~92 m in the
# floodplain). If confirmed, this is a third independent line of evidence for the
# altimetric artefact - geometric in nature and A PRIORI, since it does not depend
# on the validation points.
#
# PARAMETRIC: computes the AoA for ONE model per run, selected by the MODELLO
#   variable below. Run it once per model whose applicability is of interest
#   (L1, L2_1, L2_2, L2_4, L2_5). Each run needs the matching predictor stack
#   exported by the stack-export script.
#
# Scale and memory (16 GB): aoa() over ~16 Mpx x ~68 predictors is out of reach in
#   one call. Strategy: trainDI() ONCE on the training set (light: it computes
#   weights and threshold), then aoa() applied to row BLOCKS of the raster, reusing
#   that trainDI. Memory stays bounded and the result is identical.
# Run from a terminal (Rscript), never from an interactive IDE.

# Packages ----
suppressPackageStartupMessages({
  library(terra); library(ranger); library(data.table)
})
if (!requireNamespace("CAST", quietly = TRUE))
  stop("Package CAST missing. Install with:\n",
       "  install.packages(c('CAST','caret'))")
library(CAST)

# Model selection ----
MODELLO <- "L2_4"   # the model for this run (L1, L2_1, L2_2, L2_4 or L2_5)

# target column: the L1 model predicts the macrotype, the L2_k models the L2 class.
# Set explicitly (not by pattern matching on the column names), so the wrong column
# can never be picked when both COD_L1 and COD_L2 are present in the training table.
TARGET <- if (MODELLO == "L1") "COD_L1" else "COD_L2"

# Folders and parameters ----
# Set the project root directory here, or via the POPARK_DATA environment variable
DIR_BASE <- Sys.getenv("POPARK_DATA", unset = "path/to/project/root")
if (!dir.exists(DIR_BASE))
  stop("Invalid DIR_BASE: ", DIR_BASE,
       "\n     Set the path at the top of the script, or the POPARK_DATA variable.")

DIR_RF    <- file.path(DIR_BASE, "output_spettrale", "07_RF_output")
DIR_TRAIN <- file.path(DIR_BASE, "output_spettrale", "05_06_output")
DIR_CLASS <- file.path(DIR_BASE, "output_spettrale_intero_parco", "08_classification_output")
DIR_OUT   <- file.path(DIR_CLASS, "AOA")
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)

# whole-park predictor stack, written by the stack-export script
PATH_STACK <- file.path(DIR_CLASS, sprintf("stack_predittori_%s_intero_parco.tif", MODELLO))

# training pixels of this model (from the train/validation split)
PATH_TRAIN <- file.path(DIR_TRAIN, MODELLO, "train_pixel.csv")

RIGHE_BLOCCO <- 200L     # rows per block (lower it if RAM is tight)
NUM_CORES    <- 6L
SEED         <- 50L

# trainDI computes the dissimilarity INTERNAL to the training set: it is O(n^2).
# With ~266,000 pixels the computation is intractable (days). The threshold and the
# dissimilarity patterns are estimated in a statistically equivalent way on a
# class-stratified sample. 15-20 thousand pixels -> minutes.
N_TRAIN_MAX <- 15000L

# Model and predictors ----
cat(sprintf("Area of Applicability - model %s\n\n", MODELLO))

f_mod <- file.path(DIR_RF, MODELLO, "rf_model_final.rds")
if (!file.exists(f_mod)) stop("Model not found: ", f_mod)
mod   <- readRDS(f_mod)
preds <- mod$forest$independent.variable.names
cat(sprintf("   Predictors of model %s: %d\n", MODELLO, length(preds)))
cat(sprintf("   'dem' among the predictors: %s\n",
    if ("dem" %in% preds) "YES" else "NO -- the AoA will not capture the elevation effect"))

# Training data ----
if (!file.exists(PATH_TRAIN)) stop("Training not found: ", PATH_TRAIN)
tr <- fread(PATH_TRAIN)

manca_tr <- setdiff(preds, names(tr))
if (length(manca_tr) > 0)
  stop("Predictors absent from the training set: ", paste(manca_tr, collapse = ", "))
if (!TARGET %in% names(tr))
  stop("Target column '", TARGET, "' absent from the training set.")

# stratified subsample of the reference set (proportional quota per class)
set.seed(SEED)
if (nrow(tr) > N_TRAIN_MAX) {
  cl  <- tr[[TARGET]]
  idx <- unlist(lapply(split(seq_len(nrow(tr)), cl), function(ii) {
    # as.numeric guards against integer overflow in the product below
    n_k <- max(1L, round(as.numeric(N_TRAIN_MAX) * length(ii) / nrow(tr)))
    if (length(ii) <= n_k) ii else sample(ii, n_k)
  }))
  tr_ref <- tr[idx, ]
  cat(sprintf("   Training: %s pixels -> stratified sample of %s (for trainDI)\n",
      format(nrow(tr), big.mark = "."), format(length(idx), big.mark = ".")))
} else {
  tr_ref <- tr
  cat(sprintf("   Training pixels: %s\n", format(nrow(tr), big.mark = ".")))
}

# predictor matrix, rows ALIGNED with tr_ref (used again for the importance fallback)
train_pred <- as.data.frame(tr_ref[, ..preds])

# trainDI: weights and threshold (computed once, light) ----
# Weights come from the model importance. The final model was trained with
# importance='none', so the importance is read from variable_importance.csv if the
# training script saved it, otherwise recomputed with a quick ranger run.
vi_file <- file.path(DIR_RF, MODELLO, "variable_importance.csv")
if (file.exists(vi_file)) {
  vi_df <- read.csv(vi_file)
  names(vi_df)[1:2] <- c("var", "imp")       # expected: predictor, importance
  w <- setNames(vi_df$imp[match(preds, vi_df$var)], preds)
  w[is.na(w)] <- 0
  cat("   Importance read from variable_importance.csv\n")
} else {
  cat("   variable_importance.csv absent: recomputing importance (ranger impurity)...\n")
  set.seed(SEED)
  # sample INDICES of train_pred, so predictors and labels stay aligned: both are
  # indexed on tr_ref. Indexing the full 'tr' by the row names of a subsample would
  # mismatch the labels, since tr_ref is already a subsample of tr.
  idx_vi <- sample(nrow(train_pred), min(50000L, nrow(train_pred)))
  rf_vi  <- ranger(x = train_pred[idx_vi, , drop = FALSE],
                   y = as.factor(tr_ref[[TARGET]][idx_vi]),
                   num.trees = 200L, importance = "impurity",
                   num.threads = NUM_CORES, seed = SEED, verbose = FALSE)
  w <- rf_vi$variable.importance[preds]; w[is.na(w)] <- 0
}

cat("   Computing trainDI (applicability threshold)...\n")
tdi <- trainDI(train = train_pred, variables = preds, weight = as.data.frame(t(w)))
cat(sprintf("   AoA threshold (DI): %.4f\n", tdi$threshold))
saveRDS(tdi, file.path(DIR_OUT, sprintf("trainDI_%s.rds", MODELLO)))

# Predictor stack ----
if (!file.exists(PATH_STACK))
  stop("Predictor stack not found: ", PATH_STACK,
       "\n  Export it first with the stack-export script, for this same model.")

stk <- rast(PATH_STACK)
manca_stk <- setdiff(preds, names(stk))
if (length(manca_stk) > 0)
  stop("Predictors absent from the stack: ", paste(manca_stk, collapse = ", "))
stk <- stk[[preds]]
cat(sprintf("   Stack: %d x %d px, %d layers\n", nrow(stk), ncol(stk), nlyr(stk)))

# aoa() in row blocks ----
out_DI  <- rast(stk, nlyrs = 1); names(out_DI)  <- "DI"
out_AOA <- rast(stk, nlyrs = 1); names(out_AOA) <- "AOA"
f_DI  <- file.path(DIR_OUT, sprintf("DI_%s.tif",  MODELLO))
f_AOA <- file.path(DIR_OUT, sprintf("AOA_%s.tif", MODELLO))

nr     <- nrow(stk)
starts <- seq(1, nr, by = RIGHE_BLOCCO)
cat(sprintf("   aoa() over %d blocks of %d rows...\n", length(starts), RIGHE_BLOCCO))

invisible(writeStart(out_DI,  f_DI,  overwrite = TRUE))
invisible(writeStart(out_AOA, f_AOA, overwrite = TRUE))
readStart(stk)     # opens the read connection (required for block-wise readValues)
on.exit(try(readStop(stk), silent = TRUE), add = TRUE)

t0 <- Sys.time()
for (i in seq_along(starts)) {
  r0 <- starts[i]; nrows <- min(RIGHE_BLOCCO, nr - r0 + 1)
  v  <- terra::readValues(stk, row = r0, nrows = nrows, dataframe = TRUE)
  ok <- complete.cases(v)
  DI <- rep(NA_real_, nrow(v))
  AO <- rep(NA_real_, nrow(v))
  if (any(ok)) {
    a <- aoa(newdata = v[ok, , drop = FALSE], trainDI = tdi, verbose = FALSE)
    DI[ok] <- a$DI
    AO[ok] <- a$AOA
  }
  writeValues(out_DI,  DI, r0, nrows)
  writeValues(out_AOA, AO, r0, nrows)
  if (i %% 5 == 0 || i == 1 || i == length(starts))
    cat(sprintf("     block %d/%d  (%.0f%%)\n", i, length(starts), 100 * i / length(starts)))
}
readStop(stk)
writeStop(out_DI); writeStop(out_AOA)
cat(sprintf("   aoa() time: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# Statistics ----
aoa_r  <- rast(f_AOA)
tab    <- freq(aoa_r)
dentro <- tab$count[tab$value == 1]; fuori <- tab$count[tab$value == 0]
dentro <- if (length(dentro) == 0) 0 else dentro
fuori  <- if (length(fuori)  == 0) 0 else fuori
tot    <- dentro + fuori

cat("\nResult\n")
cat(sprintf("   Inside AoA  : %s px (%.1f%%)\n", format(dentro, big.mark = "."), 100 * dentro / tot))
cat(sprintf("   OUTSIDE AoA : %s px (%.1f%%)  <- extrapolation\n",
    format(fuori, big.mark = "."), 100 * fuori / tot))
cat(sprintf("\n   Open AOA_%s.tif in QGIS on top of the classification: the areas at 0 are\n", MODELLO))
cat("   those where the model extrapolates. Check whether the Turin floodplain is among them.\n")
cat(sprintf("\n   Saved in %s:\n     DI_%s.tif, AOA_%s.tif, trainDI_%s.rds\n",
    DIR_OUT, MODELLO, MODELLO, MODELLO))

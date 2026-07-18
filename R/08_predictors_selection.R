# Predictor selection (ROC-AUC + Boruta + Kendall).
# Reduces the 101 predictors to an informative, non-redundant subset SEPARATELY
# for each model (L1, L2_1, L2_2, L2_4, L2_5): selection is a property of the
# model, so different targets give different predictors. Three filters in
# sequence:
#   1. ROC-AUC - drops predictors that separate NO class (one-vs-rest AUC below
#      threshold for every class).
#   2. Boruta  - statistical test: compares each predictor with shadow copies
#      (randomly shuffled = pure noise) and keeps only those that significantly
#      beat the noise ("Confirmed").
#   3. Kendall - among highly correlated predictors, removes one (redundancy).
# The final number of predictors is not imposed: it emerges from the data.
# GLCM_COR predictors are excluded upfront (structural NAs, outside the RF).
# Selection runs on a per-class stratified subsample (keeps Boruta/Kendall
# tractable); reproducible via a fixed seed.

# Packages ----
library(data.table)
library(dplyr)

# Boruta/pROC/caret are extra dependencies: check they are installed.
.pkg_mancanti <- character(0)
for (pkg in c("Boruta", "pROC", "caret"))
  if (!requireNamespace(pkg, quietly = TRUE)) .pkg_mancanti <- c(.pkg_mancanti, pkg)
if (length(.pkg_mancanti) > 0)
  stop("Missing packages: ", paste(.pkg_mancanti, collapse = ", "),
       "\n  Install with: install.packages(c(",
       paste(sprintf('"%s"', .pkg_mancanti), collapse = ", "), "))")
suppressPackageStartupMessages({ library(Boruta); library(pROC); library(caret) })

# explicit namespace conflicts
select <- dplyr::select
filter <- dplyr::filter

# Parameters ----
SEED           <- 50L
AUC_MIN        <- 0.75    # filter 1: keep a predictor if it separates at least one class at AUC >= 0.75
KENDALL_CUTOFF <- 0.90    # filter 3: remove one predictor from each pair with |tau| > cutoff
BORUTA_MAXRUNS <- 60L     # Boruta iterations (higher = more robust but slower)
SUB_PER_CLASSE <- 2000L   # per-class stratified subsample (AUC + Boruta)
N_KENDALL      <- 3000L   # subsample for the Kendall correlation matrix
NUM_CORES      <- 6L

set.seed(SEED)
cat(sprintf("Parameters: AUC>=%.2f, Kendall|tau|>%.2f, Boruta maxRuns=%d, sub/class=%d\n",
    AUC_MIN, KENDALL_CUTOFF, BORUTA_MAXRUNS, SUB_PER_CLASSE))

# Folders ----
# Set the project root directory here, or via the POPARK_DATA environment variable
DIR_BASE <- Sys.getenv("POPARK_DATA", unset = "path/to/project/root")
if (!dir.exists(DIR_BASE))
  stop("Invalid DIR_BASE: ", DIR_BASE,
       "
     Set the path at the top of the script, or the POPARK_DATA variable.")
DIR_STEP4 <- file.path(DIR_BASE, "output_spettrale", "STEP4_output")
DIR_0506  <- file.path(DIR_BASE, "output_spettrale", "05_06_output")
PATH_IN   <- file.path(DIR_STEP4, "predittori_estratti_completi.csv")
if (!file.exists(PATH_IN)) stop("Input not found: ", PATH_IN)

# Model definition (target + pixel filter) ----
# L1 uses all pixels (target = macrotype COD_L1). Each L2_k uses only pixels of
# its own macrotype (COD_L1 == k) and targets the L2 class (COD_L2).
modelli <- list(
  L1   = list(target = "COD_L1", l1 = NULL),
  L2_1 = list(target = "COD_L2", l1 = 1L),
  L2_2 = list(target = "COD_L2", l1 = 2L),
  L2_4 = list(target = "COD_L2", l1 = 4L),
  L2_5 = list(target = "COD_L2", l1 = 5L)
)

# Load data and define candidates ----
cat("\nLoading", basename(PATH_IN), "...\n")
dt <- fread(PATH_IN)

COLS_ID       <- c("ID_PIXEL", "x", "y", "COD_L1", "COD_L2")
COLS_GLCM_COR <- grep("^GLCM_COR_", names(dt), value = TRUE)
# candidates = all predictors except IDs and GLCM_COR (structural NAs, outside the RF)
candidati <- setdiff(names(dt), c(COLS_ID, COLS_GLCM_COR))
cat(sprintf("   Pixels: %s | candidates: %d (excluded %d GLCM_COR)\n",
    format(nrow(dt), big.mark = "."), length(candidati), length(COLS_GLCM_COR)))

# Functions ----
# stratified subsample: up to n_per pixels per class (reproducible)
subcampiona <- function(d, target, n_per) {
  d <- as.data.table(d)
  d[, .SD[sample(.N, min(.N, n_per))], by = target]
}

# discriminating power of a predictor = maximum one-vs-rest AUC across classes.
# For each class k the AUC of separating "k vs rest" is computed; the maximum is
# taken, since a predictor is useful if it separates even a single class well.
# (max(a, 1-a) makes the measure independent of the relationship direction.)
auc_max <- function(x, y) {
  cls <- unique(y)
  aucs <- vapply(cls, function(k) {
    a <- tryCatch(
      suppressMessages(as.numeric(pROC::auc(response  = as.integer(y == k),
                                            predictor = x, quiet = TRUE))),
      error = function(e) 0.5)
    max(a, 1 - a)
  }, numeric(1))
  max(aucs, na.rm = TRUE)
}

# Loop over models ----
summary_rows <- list()

for (m in names(modelli)) {
  cfg <- modelli[[m]]
  cat(sprintf("\n-- Model %s (target=%s)\n", m, cfg$target))

  # subset of the model's pixels + valid target
  sub <- if (is.null(cfg$l1)) dt else dt[COD_L1 == cfg$l1]
  sub <- sub[!is.na(get(cfg$target))]
  if (length(unique(sub[[cfg$target]])) < 2) { cat("   <2 classes -> skip.\n"); next }
  cat(sprintf("   Model pixels: %s | classes: %d\n",
      format(nrow(sub), big.mark = "."), length(unique(sub[[cfg$target]]))))

  # per-class stratified subsample (for AUC + Boruta)
  set.seed(SEED)
  subs <- subcampiona(sub[, c(cfg$target, candidati), with = FALSE],
                      cfg$target, SUB_PER_CLASSE)
  y <- as.factor(subs[[cfg$target]])
  X <- as.data.frame(subs[, ..candidati])
  cat(sprintf("   Subsample: %s pixels\n", format(nrow(subs), big.mark = ".")))

  n_auc <- NA_integer_; n_boruta <- NA_integer_

  # filter 1 - ROC-AUC (one-vs-rest, maximum across classes)
  auc_vals <- vapply(candidati, function(p) auc_max(X[[p]], y), numeric(1))
  keep_auc <- names(auc_vals)[auc_vals >= AUC_MIN]
  n_auc    <- length(keep_auc)
  cat(sprintf("   [1] ROC-AUC >= %.2f : %d / %d predictors\n",
      AUC_MIN, n_auc, length(candidati)))

  if (n_auc < 2) {
    final <- keep_auc
    cat("   <2 predictors after AUC -> skip Boruta/Kendall.\n")
  } else {
    # filter 2 - Boruta (comparison with shadow features)
    set.seed(SEED)
    bor  <- Boruta(x = X[, keep_auc], y = y,
                   maxRuns = BORUTA_MAXRUNS, num.threads = NUM_CORES)
    bor  <- tryCatch(TentativeRoughFix(bor), error = function(e) bor)
    conf <- getSelectedAttributes(bor, withTentative = FALSE)
    n_boruta <- length(conf)
    cat(sprintf("   [2] Boruta Confirmed : %d / %d\n", n_boruta, n_auc))

    if (n_boruta < 2) {
      final <- conf
    } else {
      # filter 3 - Kendall (redundancy removal)
      set.seed(SEED)
      idx_k <- sample(nrow(X), min(nrow(X), N_KENDALL))
      cor_k <- cor(X[idx_k, conf], method = "kendall")
      cor_k[is.na(cor_k)] <- 0   # undefined correlation (constant predictor) -> 0
      rm_corr <- caret::findCorrelation(cor_k, cutoff = KENDALL_CUTOFF, names = TRUE)
      final   <- setdiff(conf, rm_corr)
      cat(sprintf("   [3] Kendall (|tau|>%.2f) removes %d -> final: %d\n",
          KENDALL_CUTOFF, length(rm_corr), length(final)))
    }
  }

  # write the model output (column 'predittore' = format read by the training script)
  dir_m <- file.path(DIR_0506, m)
  dir.create(dir_m, recursive = TRUE, showWarnings = FALSE)
  fwrite(data.table(predittore = final),
         file.path(dir_m, "predittori_selezionati.csv"))
  cat(sprintf("   -> %s/predittori_selezionati.csv (%d predictors)\n", m, length(final)))

  summary_rows[[m]] <- data.frame(
    modello       = m,
    n_candidati   = length(candidati),
    n_post_auc    = n_auc,
    n_post_boruta = n_boruta,
    n_finale      = length(final)
  )
}

# Summary ----
summ <- do.call(rbind, summary_rows)
fwrite(summ, file.path(DIR_0506, "selezione_summary.csv"))

cat("\nSelection summary:\n")
print(summ, row.names = FALSE)
cat(sprintf("\n   selezione_summary.csv -> %s\n", DIR_0506))
cat("   Per model: predittori_selezionati.csv in 05_06_output/<model>/\n")
cat("   The training script reads it automatically.\n")

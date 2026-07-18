# Train / validation split (per pixel, 70/30).
# Reads the pure pixels (training_pixels_clean.csv, output of the extraction step)
# and splits them into training (70%) and validation (30%), separately for the 5
# RF models of the hierarchy. The split is stratified by class: each class keeps
# the same 70/30 proportion in both partitions.
#
# The 5 hierarchy models:
#   L1   : all pixels, target = macrotype COD_L1 (7 macrotypes)
#   L2_1 : only COD_L1=1 (water/wetlands),    target = COD_L2 -> classes 11/12/13
#   L2_2 : only COD_L1=2 (river substrate),   target = COD_L2 -> classes 21/22
#   L2_4 : only COD_L1=4 (tree cover),        target = COD_L2 -> classes 41..46
#   L2_5 : only COD_L1=5 (cultivated),        target = COD_L2 -> classes 51..54
#   Macrotypes 3, 6, 7 have a single L2 class (31, 61, 71): no L2 model needed,
#   they are classified directly by the L1 model.
#
# Why a per-pixel split (not per polygon): a per-pixel split does not separate
# train and validation spatially (neighbouring pixels of the same polygon can end
# up in both -> slightly optimistic estimate from spatial autocorrelation,
# Roberts et al. 2017). Spatial reliability is handled downstream by (a) the AOA,
# area of applicability (Meyer & Pebesma 2021) and (b) external validation on
# independent photo-interpreted points. A per-polygon split is not produced here.

# Packages ----
library(data.table)
library(caret)
# (dplyr not needed: the per-model filter uses eval() on quote() expressions)

# Folders ----
# Set the project root directory here, or via the POPARK_DATA environment variable
DIR_BASE <- Sys.getenv("POPARK_DATA", unset = "path/to/project/root")
if (!dir.exists(DIR_BASE))
  stop("Invalid DIR_BASE: ", DIR_BASE,
       "
     Set the path at the top of the script, or the POPARK_DATA variable.")
DIR_IN   <- file.path(DIR_BASE, "output_spettrale", "STEP4_output")
DIR_OUT  <- file.path(DIR_BASE, "output_spettrale", "05_06_output")
CSV_IN   <- file.path(DIR_IN, "training_pixels_clean.csv")

if (!file.exists(CSV_IN))
  stop("training_pixels_clean.csv not found. Run the z-statistics extraction step first.")

dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

SEED       <- 50L
PROP_TRAIN <- 0.70
set.seed(SEED)

# Load CSV ----
cat("Loading training_pixels_clean.csv\n")
df <- as.data.frame(fread(CSV_IN, stringsAsFactors = FALSE))
cat(sprintf("   %s pixels x %d columns\n", format(nrow(df), big.mark="."), ncol(df)))

stopifnot(all(c("ID_PIXEL", "COD_L1", "COD_L2") %in% names(df)))

# hierarchical labels as integers (consistent with the legend codes)
df$COD_L1 <- as.integer(df$COD_L1)
df$COD_L2 <- as.integer(df$COD_L2)

cat(sprintf("   Unique COD_L1: %s\n", paste(sort(unique(df$COD_L1)), collapse=", ")))
cat(sprintf("   Unique COD_L2: %s\n", paste(sort(unique(df$COD_L2)), collapse=", ")))

# Model definition ----
# Each model is defined by name, target variable and a data filter. Filters are
# dynamic (on COD_L1/COD_L2): any final-legend class is captured automatically in
# the right model, without hand-written class lists.
modelli <- list(
  list(nome="L1",   target="COD_L1", filtro=quote(!is.na(COD_L1))),
  list(nome="L2_1", target="COD_L2", filtro=quote(COD_L1==1 & !is.na(COD_L2))),
  list(nome="L2_2", target="COD_L2", filtro=quote(COD_L1==2 & !is.na(COD_L2))),
  list(nome="L2_4", target="COD_L2", filtro=quote(COD_L1==4 & !is.na(COD_L2))),
  list(nome="L2_5", target="COD_L2", filtro=quote(COD_L1==5 & !is.na(COD_L2)))
)

# Loop over models: filter + 70/30 per-pixel split + write ----
summary_list <- list()

for (m in modelli) {

  dir_m <- file.path(DIR_OUT, m$nome)
  dir.create(dir_m, showWarnings = FALSE, recursive = TRUE)

  # filter the model rows (eval of the quote() expression in the df context)
  keep <- eval(m$filtro, envir = df)
  df_m <- df[keep, , drop = FALSE]

  target_vals <- df_m[[m$target]]

  cat(sprintf("\n-- Model %s: %s pixels | target=%s | classes: %s\n",
      m$nome, format(nrow(df_m), big.mark="."), m$target,
      paste(sort(unique(target_vals)), collapse=", ")))

  # 70/30 split stratified by class (per pixel). createDataPartition guarantees
  # each class is present in train and validation at the requested proportions.
  set.seed(SEED)
  idx_train <- caret::createDataPartition(
    factor(target_vals), p = PROP_TRAIN, list = FALSE
  )[, 1]

  df_tr <- df_m[ idx_train, , drop = FALSE]
  df_va <- df_m[-idx_train, , drop = FALSE]

  fwrite(df_tr, file.path(dir_m, "train_pixel.csv"))
  fwrite(df_va, file.path(dir_m, "validation_pixel.csv"))

  cat(sprintf("   train_pixel.csv:      %s pixels\n", format(nrow(df_tr), big.mark=".")))
  cat(sprintf("   validation_pixel.csv: %s pixels\n", format(nrow(df_va), big.mark=".")))

  # per-class diagnostics (n train / n val / % train)
  cls <- sort(unique(target_vals))
  for (c in cls) {
    n_tr <- sum(df_tr[[m$target]] == c)
    n_va <- sum(df_va[[m$target]] == c)
    n_to <- n_tr + n_va
    summary_list[[length(summary_list) + 1]] <- data.frame(
      modello   = m$nome,
      target    = m$target,
      classe    = c,
      n_train   = n_tr,
      n_val     = n_va,
      n_tot     = n_to,
      pct_train = round(100 * n_tr / n_to, 1)
    )
  }
}

# Split diagnostics ----
split_summary <- do.call(rbind, summary_list)
fwrite(split_summary, file.path(DIR_OUT, "split_summary.csv"))

cat("\nTrain/validation split summary:\n")
cat(sprintf("   Models produced: %s\n",
    paste(sapply(modelli, function(x) x$nome), collapse=", ")))
cat(sprintf("   Output in: %s\n", DIR_OUT))
cat("   Per model: train_pixel.csv + validation_pixel.csv\n")
cat("   Diagnostics: split_summary.csv\n")
cat("   Split: per pixel only (70/30 stratified). Split B not produced.\n")

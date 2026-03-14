#Runs external cross validation of MintTea run


options(stringsAsFactors = FALSE)

setwd("C:/PhD/Canterbury dataset/DIABLO")


#Runs on the species results change to genus for genus runs
minttea_dir <- "mixomics_outputs_2block_species_final/minttea_outputs"
out_rds <- "mixomics_outputs_2block_species_final/minttea_input_microbiome_metabolome.rds"
selected_minttea_setting <- "keep_20//des_0.3//nrep_10//nfol_5//ncom_5//edge_0.9"
selected_minttea_modules <- c("module1", "module2", "module3")
repeats <- 100
folds <- 5
seed <- 20260306

minttea_result_rds <- file.path(minttea_dir, "minttea_result.rds")
out_results_csv <- file.path(minttea_dir, "selected_setting_repeated_external_cv_results.csv")
out_summary_csv <- file.path(minttea_dir, "selected_setting_repeated_external_cv_summary.csv")

if (!file.exists(minttea_result_rds)) stop("Missing MintTea RDS: ", minttea_result_rds)
if (!file.exists(out_rds)) stop("Missing input RDS: ", out_rds)

minttea_obj <- readRDS(minttea_result_rds)
out_df <- readRDS(out_rds)

if (!(selected_minttea_setting %in% names(minttea_obj))) {
  stop("Selected setting not found: ", selected_minttea_setting)
}

compute_auc_rank <- function(y_true_binary01, y_score) {
  y <- as.integer(y_true_binary01)
  s <- suppressWarnings(as.numeric(y_score))
  ok <- is.finite(s) & !is.na(y)
  y <- y[ok]
  s <- s[ok]
  if (length(y) < 2) return(NA_real_)
  n_pos <- sum(y == 1L)
  n_neg <- sum(y == 0L)
  if (n_pos == 0L || n_neg == 0L) return(NA_real_)
  r <- rank(s, ties.method = "average")
  (sum(r[y == 1L]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

make_stratified_fold_ids <- function(y, k, seed = NULL) {
  y <- as.character(y)
  out <- integer(length(y))
  lev <- unique(y)
  if (!is.null(seed)) set.seed(seed)
  for (cl in lev) {
    idx <- which(y == cl)
    idx <- sample(idx, length(idx), replace = FALSE)
    out[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  out
}

run_repeated_external_cv_for_module <- function(proc_data, feature_names, repeats, folds, seed) {
  y_chr <- as.character(proc_data$disease_state)
  y <- ifelse(y_chr == "disease", 1L, 0L)

  min_class_n <- min(table(y))
  k_use <- min(as.integer(folds), as.integer(min_class_n))
  if (k_use < 2L) stop("Too few samples per class for CV.")

  X <- as.matrix(proc_data[, feature_names, drop = FALSE])
  X <- apply(X, 2, function(v) suppressWarnings(as.numeric(v)))
  X <- as.matrix(X)

  rows <- list()
  rr <- 1L

  for (rep_i in seq_len(as.integer(repeats))) {
    fold_id <- make_stratified_fold_ids(y, k = k_use, seed = as.integer(seed) + rep_i)
    for (fold_i in seq_len(k_use)) {
      test_idx <- which(fold_id == fold_i)
      train_idx <- which(fold_id != fold_i)
      y_train <- y[train_idx]
      y_test <- y[test_idx]
      if (length(unique(y_train)) < 2 || length(unique(y_test)) < 2) next

      x_train <- X[train_idx, , drop = FALSE]
      x_test <- X[test_idx, , drop = FALSE]

      pca <- try(stats::prcomp(x_train, center = TRUE, scale. = FALSE), silent = TRUE)
      if (inherits(pca, "try-error") || ncol(pca$x) < 1) next

      train_pc1 <- as.numeric(pca$x[, 1])
      test_pc <- try(stats::predict(pca, newdata = x_test), silent = TRUE)
      if (inherits(test_pc, "try-error") || ncol(test_pc) < 1) next
      test_pc1 <- as.numeric(test_pc[, 1])

      fit <- try(stats::glm(y_train ~ train_pc1, family = stats::binomial()), silent = TRUE)
      if (inherits(fit, "try-error")) next

      prob_test <- try(stats::predict(fit, newdata = data.frame(train_pc1 = test_pc1), type = "response"), silent = TRUE)
      if (inherits(prob_test, "try-error")) next

      rows[[rr]] <- data.frame(
        repeat_id = rep_i,
        fold_id = fold_i,
        n_train = length(train_idx),
        n_test = length(test_idx),
        auc = compute_auc_rank(y_test, prob_test),
        stringsAsFactors = FALSE
      )
      rr <- rr + 1L
    }
  }

  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

setting_obj <- minttea_obj[[selected_minttea_setting]]
available_modules <- names(setting_obj)
module_names <- intersect(selected_minttea_modules, available_modules)
module_names <- module_names[grepl("^module", module_names)]
if (length(module_names) == 0) stop("No selected modules found under selected setting.")

cv_rows <- list()
cv_summ <- list()
kk <- 1L
ss <- 1L

for (m in module_names) {
  mod_obj <- setting_obj[[m]]
  feats <- as.character(mod_obj$features)
  feats <- feats[feats %in% colnames(out_df)]
  if (length(feats) < 2) next

  cv_df <- run_repeated_external_cv_for_module(out_df, feats, repeats, folds, seed)
  if (nrow(cv_df) == 0) next

  cv_df$settings <- selected_minttea_setting
  cv_df$module <- m
  cv_df$module_size <- length(feats)
  cv_rows[[kk]] <- cv_df
  kk <- kk + 1L

  fold_auc <- cv_df$auc[is.finite(cv_df$auc)]
  cv_summ[[ss]] <- data.frame(
    settings = selected_minttea_setting,
    module = m,
    module_size = length(feats),
    n_rows = nrow(cv_df),
    auc_mean = if (length(fold_auc) > 0) mean(fold_auc) else NA_real_,
    auc_sd = if (length(fold_auc) > 1) stats::sd(fold_auc) else NA_real_,
    auc_q025 = if (length(fold_auc) > 0) as.numeric(stats::quantile(fold_auc, 0.025, na.rm = TRUE)) else NA_real_,
    auc_q50 = if (length(fold_auc) > 0) as.numeric(stats::quantile(fold_auc, 0.5, na.rm = TRUE)) else NA_real_,
    auc_q975 = if (length(fold_auc) > 0) as.numeric(stats::quantile(fold_auc, 0.975, na.rm = TRUE)) else NA_real_,
    stringsAsFactors = FALSE
  )
  ss <- ss + 1L
}

if (length(cv_rows) == 0) stop("No valid repeated external CV outputs were produced.")

cv_out <- do.call(rbind, cv_rows)
cv_sum <- do.call(rbind, cv_summ)
cv_sum <- cv_sum[order(cv_sum$auc_mean, decreasing = TRUE, na.last = TRUE), , drop = FALSE]

utils::write.csv(cv_out, out_results_csv, row.names = FALSE)
utils::write.csv(cv_sum, out_summary_csv, row.names = FALSE)

cat("Repeated external CV results CSV:", out_results_csv, "\n")
cat("Repeated external CV summary CSV:", out_summary_csv, "\n")
cat("Top module by repeated external CV AUC:\n")
print(utils::head(cv_sum, 1))

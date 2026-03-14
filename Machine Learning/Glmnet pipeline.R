# TWO-BLOCK GROUP DISCRIMINATION (Late fusion)
#   Multinomial elastic net (glmnet) on each block separately,
#   then fuse predicted class probabilities.

# Inputs:
#   Preferred: aligned_blocks.rds
#     (created by Aligned block generation.R)

# Outputs (late_fusion_glmnet/):
#   - LF_perf_summary.txt
#   - LF_perf_by_repeat.csv
#   - LF_confusion_matrix_pooled.csv
#   - LF_alpha_choices.csv
#   - LF_perf_compact_table.csv
#   - LF_feature_stability_microbiome.csv
#   - LF_feature_stability_metabolome.csv
#   - LF_microbiome_core_features.csv
#   - LF_metabolome_core_features.csv
#   - LF_driver_features_top.csv
#   - LF_confusion_metrics_pooled.csv
#   - OOF score distribution plot (metabolome)
#   - LF_oof_predictions.csv
#   - ROC/PR curve SVGs (per-block + combined)
#   - LF_permutation_summary.csv + LF_permutation_results.rds


# Notes:
#   - This uses repeated stratified M-fold CV.
#   - Hyperparameters (alpha, lambda) are tuned INSIDE each outer training split.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(glmnet)
  library(pROC)
  library(ggplot2)
  library(PRROC)
})

setwd("C:/PhD/Canterbury dataset/DIABLO")

# USER CONFIG

base_out_dir   <- "mixomics_outputs_2block"
aligned_rds    <- file.path(base_out_dir, "aligned_blocks.rds")
out_dir        <- file.path(base_out_dir, "late_fusion_glmnet")

# Outer CV
outer_folds    <- 3          # recommend 3 or 4 with n=28 and 8 treated
outer_repeats  <- 200        # increase for stability
outer_seed     <- 1

# Inner CV (for alpha/lambda selection)
alpha_grid     <- c(0, 0.25, 0.5, 0.75, 1)  # 0=ridge, 1=lasso
inner_folds    <- 3

# Fusion strategy
fusion_mode    <- "inner_accuracy"     # "equal" or "inner_accuracy"
fusion_policy  <- "metabolome_primary" # "always" or "metabolome_primary"

# Threshold tuning (binary only)
positive_class   <- "CeD"
threshold_metric <- "balanced_accuracy"  # "balanced_accuracy" or "youden"
tune_threshold   <- TRUE

# Variance filtering 
# Values are SD cutoffs; set to 0 to disable for a block.
variance_filter_sd <- list(
  Microbiome = 0,
  Metabolome = 0.05
)

# Feature stability reporting
core_feature_prop_cutoff <- 0.6   # e.g., show metabolites selected in >=60% of refits
driver_top_n <- 10

# Permutation testing + plotting
n_permutations <- 500
perm_seed <- 123
roc_plot_dims <- list(width = 6, height = 5)
pr_plot_dims  <- list(width = 6, height = 5)


# HELPERS


stop_quiet <- function(...) stop(paste0(...), call. = FALSE)

make_stratified_fold_ids <- function(y, k, seed) {
  # Returns integer vector fold_id in 1..k; stratified by class.
  set.seed(seed)
  y <- as.factor(y)
  n <- length(y)
  fold_id <- integer(n)

  for (lev in levels(y)) {
    idx <- which(y == lev)
    idx <- sample(idx, length(idx), replace = FALSE)
    # round-robin assignment
    fold_id[idx] <- rep(1:k, length.out = length(idx))
  }
  fold_id
}

scale_train_test <- function(X_train, X_test) {
  mu <- colMeans(X_train)
  sdv <- apply(X_train, 2, sd)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  Xtr <- sweep(X_train, 2, mu, "-")
  Xtr <- sweep(Xtr, 2, sdv, "/")
  Xte <- sweep(X_test, 2, mu, "-")
  Xte <- sweep(Xte, 2, sdv, "/")
  list(X_train = Xtr, X_test = Xte, mu = mu, sd = sdv)
}

as_prob_matrix <- function(pred_obj, classes) {
  # glmnet multinomial predict(type="response") returns an array:
  #   n x K x 1 (for one lambda)
  if (is.null(pred_obj)) return(NULL)
  arr <- pred_obj
  if (is.vector(arr)) {
    # should not happen for multinomial, but keep safe
    mat <- matrix(arr, ncol = length(classes))
    colnames(mat) <- classes
    return(mat)
  }
  if (length(dim(arr)) == 3) {
    mat <- arr[, , 1, drop = FALSE]
    # convert to 2D
    mat <- matrix(mat, nrow = dim(arr)[1], ncol = dim(arr)[2])
    colnames(mat) <- dimnames(arr)[[2]]
    return(mat)
  }
  if (length(dim(arr)) == 2) {
    mat <- arr
    colnames(mat) <- colnames(mat) %||% classes
    return(mat)
  }
  NULL
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

resolve_positive_class <- function(classes, preferred) {
  if (preferred %in% classes) return(preferred)
  if (length(classes) >= 2) return(classes[2])
  classes[1]
}

compute_class_weights <- function(y_train) {
  y_train <- as.factor(y_train)
  tab <- table(y_train)
  w <- 1 / as.numeric(tab)
  names(w) <- names(tab)
  w <- w / mean(w)
  w[as.character(y_train)]
}

apply_variance_filter <- function(X_train_full, X_test_full, tol) {
  if (is.null(tol) || tol <= 0 || ncol(X_train_full) == 0) {
    keep_idx <- seq_len(ncol(X_train_full))
  } else {
    sds <- apply(X_train_full, 2, sd)
    keep_idx <- which(is.finite(sds) & sds > tol)
    if (length(keep_idx) == 0) keep_idx <- seq_len(ncol(X_train_full))
  }
  list(
    X_train = if (length(keep_idx) == 0) X_train_full else X_train_full[, keep_idx, drop = FALSE],
    X_test = if (length(keep_idx) == 0) X_test_full else X_test_full[, keep_idx, drop = FALSE],
    keep_idx = keep_idx
  )
}

extract_beta_vector <- function(cvfit, pos_class, feature_names) {
  cf <- tryCatch(coef(cvfit, s = "lambda.min"), error = function(e) NULL)
  if (is.null(cf)) return(rep(0, length(feature_names)))
  mat <- if (is.list(cf)) cf[[pos_class]] %||% cf[[1]] else cf
  if (is.null(mat)) return(rep(0, length(feature_names)))
  mat <- as.matrix(mat)
  rn <- rownames(mat)
  vals <- as.numeric(mat)
  keep <- rn != "(Intercept)"
  rn <- rn[keep]
  vals <- vals[keep]
  beta <- setNames(vals, rn)
  out <- beta[feature_names]
  out[is.na(out)] <- 0
  out
}

summarize_beta_matrix <- function(B) {
  if (is.null(B) || nrow(B) == 0) {
    return(data.frame(feature = character(0), mean_beta = numeric(0), mean_abs_beta = numeric(0)))
  }
  valid <- !is.na(B)
  denom <- rowSums(valid)
  B2 <- B
  B2[!valid] <- 0
  mean_beta <- ifelse(denom > 0, rowSums(B2) / denom, 0)
  absB <- abs(B2)
  mean_abs_beta <- ifelse(denom > 0, rowSums(absB) / denom, 0)
  data.frame(
    feature = rownames(B),
    mean_beta = mean_beta,
    mean_abs_beta = mean_abs_beta,
    stringsAsFactors = FALSE
  )
}

compute_binary_stats <- function(y_true, y_pred, positive_class) {
  y_true <- as.factor(y_true)
  y_pred <- as.factor(y_pred)

  ok <- !is.na(y_true) & !is.na(y_pred)
  y_true <- y_true[ok]
  y_pred <- y_pred[ok]
  if (length(unique(y_true)) < 2) {
    return(list(tp = NA_integer_, tn = NA_integer_, fp = NA_integer_, fn = NA_integer_,
                sensitivity = NA_real_, specificity = NA_real_, balanced_accuracy = NA_real_))
  }

  neg_class <- setdiff(levels(y_true), positive_class)
  if (length(neg_class) != 1) {
    return(list(tp = NA_integer_, tn = NA_integer_, fp = NA_integer_, fn = NA_integer_,
                sensitivity = NA_real_, specificity = NA_real_, balanced_accuracy = NA_real_))
  }

  tp <- sum(y_true == positive_class & y_pred == positive_class)
  tn <- sum(y_true == neg_class & y_pred == neg_class)
  fp <- sum(y_true == neg_class & y_pred == positive_class)
  fn <- sum(y_true == positive_class & y_pred == neg_class)

  sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  balanced_accuracy <- mean(c(sensitivity, specificity), na.rm = TRUE)

  list(tp = tp, tn = tn, fp = fp, fn = fn,
       sensitivity = sensitivity, specificity = specificity,
       balanced_accuracy = balanced_accuracy)
}

permutation_test_auc <- function(y_true, prob_pos, positive_class, n_perm, seed) {
  if (n_perm <= 0) return(NULL)
  ok <- is.finite(prob_pos) & !is.na(y_true)
  y_true <- y_true[ok]
  prob_pos <- prob_pos[ok]
  if (length(unique(y_true)) < 2) return(NULL)
  obs <- compute_auc(y_true, prob_pos, positive_class)
  if (!is.finite(obs)) return(NULL)
  set.seed(seed)
  null_vals <- replicate(n_perm, {
    y_perm <- sample(y_true)
    compute_auc(y_perm, prob_pos, positive_class)
  })
  p_val <- (1 + sum(null_vals >= obs)) / (n_perm + 1)
  list(observed = obs, null = null_vals, p_value = p_val)
}

permutation_test_ber <- function(y_true, y_pred, classes, n_perm, seed) {
  if (n_perm <= 0) return(NULL)
  ok <- !is.na(y_true) & !is.na(y_pred)
  y_true <- y_true[ok]
  y_pred <- y_pred[ok]
  if (length(unique(y_true)) < 2) return(NULL)
  obs <- compute_metrics(y_true, y_pred, classes)$ber
  if (!is.finite(obs)) return(NULL)
  set.seed(seed)
  null_vals <- replicate(n_perm, {
    y_perm <- sample(y_true)
    compute_metrics(y_perm, y_pred, classes)$ber
  })
  p_val <- (1 + sum(null_vals <= obs)) / (n_perm + 1)
  list(observed = obs, null = null_vals, p_value = p_val)
}

fit_multinom_enet_inner <- function(X_train, y_train, alpha_grid, inner_folds, weights) {
  # Fits cv.glmnet for each alpha; selects alpha+lambda with min CV misclass.
  y_train <- as.factor(y_train)
  best <- list(alpha = NA_real_, cv = NULL, cvm = Inf)

  for (a in alpha_grid) {
    cvfit <- tryCatch({
      cv.glmnet(
        x = X_train,
        y = y_train,
        family = "multinomial",
        alpha = a,
        nfolds = inner_folds,
        type.measure = "class",
        standardize = FALSE,
        weights = weights
      )
    }, error = function(e) NULL)

    if (is.null(cvfit)) next

    # cvfit$cvm is mean misclass; choose lambda.min
    cvm_min <- min(cvfit$cvm, na.rm = TRUE)
    if (is.finite(cvm_min) && cvm_min < best$cvm) {
      best <- list(alpha = a, cv = cvfit, cvm = cvm_min)
    }
  }

  if (is.null(best$cv)) return(NULL)
  best$lambda_min <- best$cv$lambda.min
  best
}

predict_probs <- function(cvfit, X_test, classes) {
  p <- tryCatch({
    predict(cvfit, newx = X_test, s = "lambda.min", type = "response")
  }, error = function(e) NULL)
  as_prob_matrix(p, classes)
}

prob_to_class <- function(P) {
  if (is.null(P) || nrow(P) == 0) return(character(0))
  colnames(P)[max.col(P, ties.method = "first")]
}

prob_to_class_threshold <- function(P, positive_class, threshold, use_threshold) {
  if (is.null(P) || nrow(P) == 0) return(character(0))
  if (!isTRUE(use_threshold) || is.na(threshold) || !(positive_class %in% colnames(P))) {
    return(prob_to_class(P))
  }
  pos <- P[, positive_class]
  neg <- setdiff(colnames(P), positive_class)
  if (length(neg) != 1) return(prob_to_class(P))
  pred <- ifelse(pos >= threshold, positive_class, neg)
  as.character(pred)
}

compute_metrics <- function(y_true, y_pred, classes) {
  y_true <- factor(y_true, levels = classes)
  y_pred <- factor(y_pred, levels = classes)

  ok <- !is.na(y_true) & !is.na(y_pred)
  y_true <- y_true[ok]
  y_pred <- y_pred[ok]

  cm <- table(y_true, y_pred)
  overall_err <- mean(y_true != y_pred)

  per_class_err <- sapply(classes, function(cl) {
    idx <- which(y_true == cl)
    if (length(idx) == 0) return(NA_real_)
    mean(y_pred[idx] != cl)
  })
  ber <- mean(per_class_err, na.rm = TRUE)

  list(overall_err = overall_err, ber = ber, per_class_err = per_class_err, cm = cm)
}

extract_prob_pos <- function(P, positive_class) {
  if (is.null(P) || nrow(P) == 0) return(numeric(0))
  if (!(positive_class %in% colnames(P))) return(rep(NA_real_, nrow(P)))
  as.numeric(P[, positive_class])
}

inner_cv_probabilities <- function(X_train, y_train, alpha, lambda, inner_folds, classes, positive_class) {
  y_train <- as.factor(y_train)
  fold_id <- make_stratified_fold_ids(y_train, k = inner_folds, seed = 1001)
  prob_pos <- rep(NA_real_, length(y_train))

  for (f in seq_len(inner_folds)) {
    te <- which(fold_id == f)
    tr <- which(fold_id != f)

    w_tr <- compute_class_weights(y_train[tr])
    fit <- tryCatch({
      glmnet(
        x = X_train[tr, , drop = FALSE],
        y = y_train[tr],
        family = "multinomial",
        alpha = alpha,
        lambda = lambda,
        standardize = FALSE,
        weights = w_tr
      )
    }, error = function(e) NULL)

    if (is.null(fit)) next

    P <- tryCatch({
      predict(fit, newx = X_train[te, , drop = FALSE], s = lambda, type = "response")
    }, error = function(e) NULL)

    Pm <- as_prob_matrix(P, classes)
    if (is.null(Pm)) next
    prob_pos[te] <- extract_prob_pos(Pm, positive_class)
  }

  prob_pos
}

choose_threshold <- function(y_true, prob_pos, positive_class, metric = "balanced_accuracy") {
  ok <- is.finite(prob_pos) & !is.na(y_true)
  y_true <- as.factor(y_true[ok])
  prob_pos <- prob_pos[ok]
  if (length(unique(y_true)) < 2) return(list(threshold = NA_real_, score = NA_real_))

  cand <- unique(prob_pos)
  if (length(cand) > 200) cand <- seq(0.01, 0.99, by = 0.01)

  best_t <- NA_real_
  best_s <- -Inf
  neg_class <- setdiff(levels(y_true), positive_class)
  if (length(neg_class) != 1) return(list(threshold = NA_real_, score = NA_real_))

  for (t in cand) {
    pred <- ifelse(prob_pos >= t, positive_class, neg_class)
    pred <- factor(pred, levels = levels(y_true))
    y_true_f <- factor(y_true, levels = levels(y_true))

    tp <- sum(pred == positive_class & y_true_f == positive_class)
    tn <- sum(pred == neg_class & y_true_f == neg_class)
    fp <- sum(pred == positive_class & y_true_f == neg_class)
    fn <- sum(pred == neg_class & y_true_f == positive_class)

    sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
    if (!is.finite(sens) || !is.finite(spec)) next

    score <- if (metric == "youden") (sens + spec - 1) else (sens + spec) / 2
    if (score > best_s) {
      best_s <- score
      best_t <- t
    }
  }

  list(threshold = best_t, score = best_s)
}

compute_auc <- function(y_true, prob_pos, positive_class) {
  ok <- is.finite(prob_pos) & !is.na(y_true)
  y_true <- as.factor(y_true[ok])
  prob_pos <- prob_pos[ok]
  if (length(unique(y_true)) < 2) return(NA_real_)
  pos <- y_true == positive_class
  n_pos <- sum(pos)
  n_neg <- sum(!pos)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  r <- rank(prob_pos, ties.method = "average")
  auc <- (sum(r[pos]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
  as.numeric(auc)
}

compute_pr_auc <- function(y_true, prob_pos, positive_class) {
  ok <- is.finite(prob_pos) & !is.na(y_true)
  y_true <- as.factor(y_true[ok])
  prob_pos <- prob_pos[ok]
  if (length(unique(y_true)) < 2) return(NA_real_)

  ord <- order(prob_pos, decreasing = TRUE)
  y_sorted <- y_true[ord] == positive_class

  tp <- cumsum(y_sorted)
  fp <- cumsum(!y_sorted)
  prec <- tp / pmax(1, tp + fp)
  rec <- tp / max(1, sum(y_sorted))

  rec <- c(0, rec)
  prec <- c(1, prec)

  sum((rec[-1] - rec[-length(rec)]) * prec[-1])
}

extract_selected_features <- function(cvfit, classes) {
  # Returns character vector of features with any non-zero coef in any class
  # (excluding intercept)
  cf <- tryCatch(coef(cvfit, s = "lambda.min"), error = function(e) NULL)
  if (is.null(cf)) return(character(0))

  feats <- character(0)
  # coef() returns a list of sparse matrices, one per class
  for (cl in names(cf)) {
    m <- cf[[cl]]
    rn <- rownames(m)
    v <- as.numeric(m)
    nz <- which(v != 0)
    if (length(nz) == 0) next
    f <- rn[nz]
    f <- setdiff(f, "(Intercept)")
    feats <- union(feats, f)
  }
  feats
}


# LOAD DATA

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
sink(file.path(out_dir, "LF_run_log.txt"), split = TRUE)

cat("=== Late fusion multinomial elastic net (glmnet) ===\n")
cat("aligned_rds:", aligned_rds, "\n")

if (!file.exists(aligned_rds)) {
  stop_quiet(
    "Missing ", aligned_rds, ". Run two_block_run_final_models.R first to create aligned_blocks.rds, ",
    "or modify this script to read your matrices directly."
  )
}

ab <- readRDS(aligned_rds)
if (!all(c("Microbiome", "Metabolome", "Group") %in% names(ab))) {
  stop_quiet("aligned_blocks.rds must contain Microbiome, Metabolome, Group.")
}

X_micro <- ab$Microbiome
X_metab <- ab$Metabolome
y <- droplevels(as.factor(ab$Group))

if (nrow(X_micro) != length(y) || nrow(X_metab) != length(y)) {
  stop_quiet("Row counts of blocks must match length(Group).")
}
if (!identical(rownames(X_micro), rownames(X_metab))) {
  stop_quiet("Row names of Microbiome and Metabolome blocks must match.")
}

if (is.null(rownames(X_micro)) || any(rownames(X_micro) == "")) {
  assigned_ids <- paste0("Sample_", seq_len(nrow(X_micro)))
  rownames(X_micro) <- assigned_ids
  rownames(X_metab) <- assigned_ids
}

sample_ids <- rownames(X_micro)

classes <- levels(y)
pos_class <- resolve_positive_class(classes, positive_class)
chance_acc <- 1 / max(1, length(classes))

if (length(classes) != 2) {
  tune_threshold <- FALSE
  cat("Note: threshold tuning only supported for binary classification. Using max-prob class.\n")
}
cat("Group counts:\n")
print(table(y))


# OUTER CV LOOP

all_preds_micro <- character(0)
all_preds_metab <- character(0)
all_preds_fused <- character(0)
all_true <- character(0)
all_prob_micro <- numeric(0)
all_prob_metab <- numeric(0)
all_prob_fused <- numeric(0)

oof_rows <- list()

n_models_total <- outer_repeats * outer_folds
micro_beta_mat <- matrix(NA_real_, nrow = ncol(X_micro), ncol = n_models_total,
                         dimnames = list(colnames(X_micro), NULL))
metab_beta_mat <- matrix(NA_real_, nrow = ncol(X_metab), ncol = n_models_total,
                         dimnames = list(colnames(X_metab), NULL))

pos_label_clean <- gsub("[^A-Za-z0-9]+", "", pos_class)
prob_col_names <- list(
  micro = paste0("p_", pos_label_clean, "_microbiome"),
  metab = paste0("p_", pos_label_clean, "_metabolome"),
  fused = paste0("p_", pos_label_clean, "_fused")
)

# store per-repeat summaries
rep_rows <- list()

# alpha choices
alpha_rows <- list()

# threshold choices
threshold_rows <- list()

# feature stability counters
micro_feat_counts <- setNames(integer(ncol(X_micro)), colnames(X_micro))
metab_feat_counts <- setNames(integer(ncol(X_metab)), colnames(X_metab))

total_models_fit_micro <- 0L
total_models_fit_metab <- 0L

for (r in seq_len(outer_repeats)) {
  fold_id <- make_stratified_fold_ids(y, k = outer_folds, seed = outer_seed + r)

  rep_true <- character(0)
  rep_pred_micro <- character(0)
  rep_pred_metab <- character(0)
  rep_pred_fused <- character(0)
  rep_prob_micro <- numeric(0)
  rep_prob_metab <- numeric(0)
  rep_prob_fused <- numeric(0)

  for (f in seq_len(outer_folds)) {
    te <- which(fold_id == f)
    tr <- which(fold_id != f)

    model_idx <- (r - 1) * outer_folds + f
    if (!is.null(micro_beta_mat)) micro_beta_mat[, model_idx] <- NA_real_
    if (!is.null(metab_beta_mat)) metab_beta_mat[, model_idx] <- NA_real_

    y_tr <- y[tr]
    y_te <- y[te]
    w_tr <- compute_class_weights(y_tr)

    # --- Microbiome model ---
    micro_train_full <- X_micro[tr, , drop = FALSE]
    micro_test_full <- X_micro[te, , drop = FALSE]
    vf_m <- apply_variance_filter(micro_train_full, micro_test_full, variance_filter_sd$Microbiome)
    sc_m <- scale_train_test(vf_m$X_train, vf_m$X_test)
    micro_keep_idx <- vf_m$keep_idx
    micro_feat_names <- if (length(micro_keep_idx) > 0) colnames(X_micro)[micro_keep_idx] else character(0)
    fit_m <- fit_multinom_enet_inner(sc_m$X_train, y_tr, alpha_grid = alpha_grid,
                                     inner_folds = inner_folds, weights = w_tr)

    Pm <- NULL
    w_m <- NA_real_
    th_m <- NA_real_
    prob_in_m <- NULL
    if (!is.null(fit_m)) {
      Pm <- predict_probs(fit_m$cv, sc_m$X_test, classes)
      total_models_fit_micro <- total_models_fit_micro + 1L
      sel_m <- extract_selected_features(fit_m$cv, classes)
      micro_feat_counts[sel_m] <- micro_feat_counts[sel_m] + 1L
      w_m <- 1 - fit_m$cvm

      if (ncol(X_micro) > 0) {
        beta_vec <- rep(0, ncol(X_micro))
        if (length(micro_feat_names) > 0) {
          beta_local <- extract_beta_vector(fit_m$cv, pos_class, micro_feat_names)
          beta_vec[micro_keep_idx] <- beta_local
        }
        micro_beta_mat[, model_idx] <- beta_vec
      }

      if (isTRUE(tune_threshold)) {
        prob_in_m <- inner_cv_probabilities(sc_m$X_train, y_tr, fit_m$alpha, fit_m$lambda_min,
                                            inner_folds, classes, pos_class)
        th_m <- choose_threshold(y_tr, prob_in_m, pos_class, metric = threshold_metric)$threshold
      }
    }

    # --- Metabolome model ---
    metab_train_full <- X_metab[tr, , drop = FALSE]
    metab_test_full <- X_metab[te, , drop = FALSE]
    vf_t <- apply_variance_filter(metab_train_full, metab_test_full, variance_filter_sd$Metabolome)
    sc_t <- scale_train_test(vf_t$X_train, vf_t$X_test)
    metab_keep_idx <- vf_t$keep_idx
    metab_feat_names <- if (length(metab_keep_idx) > 0) colnames(X_metab)[metab_keep_idx] else character(0)
    fit_t <- fit_multinom_enet_inner(sc_t$X_train, y_tr, alpha_grid = alpha_grid,
                                     inner_folds = inner_folds, weights = w_tr)

    Pt <- NULL
    w_t <- NA_real_
    th_t <- NA_real_
    prob_in_t <- NULL
    if (!is.null(fit_t)) {
      Pt <- predict_probs(fit_t$cv, sc_t$X_test, classes)
      total_models_fit_metab <- total_models_fit_metab + 1L
      sel_t <- extract_selected_features(fit_t$cv, classes)
      metab_feat_counts[sel_t] <- metab_feat_counts[sel_t] + 1L
      w_t <- 1 - fit_t$cvm

      if (ncol(X_metab) > 0) {
        beta_vec <- rep(0, ncol(X_metab))
        if (length(metab_feat_names) > 0) {
          beta_local <- extract_beta_vector(fit_t$cv, pos_class, metab_feat_names)
          beta_vec[metab_keep_idx] <- beta_local
        }
        metab_beta_mat[, model_idx] <- beta_vec
      }

      if (isTRUE(tune_threshold)) {
        prob_in_t <- inner_cv_probabilities(sc_t$X_train, y_tr, fit_t$alpha, fit_t$lambda_min,
                                            inner_folds, classes, pos_class)
        th_t <- choose_threshold(y_tr, prob_in_t, pos_class, metric = threshold_metric)$threshold
      }
    }

    # record alpha choices
    alpha_rows[[length(alpha_rows) + 1]] <- data.frame(
      repeat_id = r,
      fold = f,
      block = c("Microbiome", "Metabolome"),
      alpha = c(ifelse(is.null(fit_m), NA_real_, fit_m$alpha), ifelse(is.null(fit_t), NA_real_, fit_t$alpha)),
      inner_misclass = c(ifelse(is.null(fit_m), NA_real_, fit_m$cvm), ifelse(is.null(fit_t), NA_real_, fit_t$cvm)),
      stringsAsFactors = FALSE
    )

    # --- Fuse probabilities ---
    # If one block failed, fall back to the other.
    use_metab_only <- isTRUE(fusion_policy == "metabolome_primary") && !is.null(Pt) &&
      (is.null(Pm) || !is.finite(w_m) || (is.finite(w_t) && w_t >= w_m) ||
         (is.finite(w_m) && w_m <= chance_acc))

    if (is.null(Pm) && is.null(Pt)) {
      # worst case: uniform
      Pf <- matrix(1 / length(classes), nrow = length(te), ncol = length(classes))
      colnames(Pf) <- classes
      th_f <- NA_real_
    } else if (is.null(Pm)) {
      Pf <- Pt
      th_f <- th_t
    } else if (is.null(Pt)) {
      Pf <- Pm
      th_f <- th_m
    } else if (use_metab_only) {
      Pf <- Pt
      th_f <- th_t
    } else {
      if (fusion_mode == "inner_accuracy" && is.finite(w_m) && is.finite(w_t) && (w_m + w_t) > 0) {
        wm <- if (w_m <= chance_acc) 0 else max(w_m, 1e-6)
        wt <- max(w_t, 1e-6)
        if (wm == 0 && wt > 0) {
          Pf <- Pt
        } else {
          Pf <- (wm * Pm + wt * Pt) / (wm + wt)
        }
      } else {
        Pf <- (Pm + Pt) / 2
      }
      th_f <- NA_real_
      if (isTRUE(tune_threshold) && !is.null(prob_in_m) && !is.null(prob_in_t)) {
        if (fusion_mode == "inner_accuracy" && is.finite(w_m) && is.finite(w_t) && (w_m + w_t) > 0) {
          wm <- if (w_m <= chance_acc) 0 else max(w_m, 1e-6)
          wt <- max(w_t, 1e-6)
          if (wm == 0 && wt > 0) {
            prob_in_f <- prob_in_t
          } else {
            prob_in_f <- (wm * prob_in_m + wt * prob_in_t) / (wm + wt)
          }
        } else {
          prob_in_f <- (prob_in_m + prob_in_t) / 2
        }
        th_f <- choose_threshold(y_tr, prob_in_f, pos_class, metric = threshold_metric)$threshold
      }
    }

    pred_m <- if (!is.null(Pm)) prob_to_class_threshold(Pm, pos_class, th_m, tune_threshold) else rep(NA_character_, length(te))
    pred_t <- if (!is.null(Pt)) prob_to_class_threshold(Pt, pos_class, th_t, tune_threshold) else rep(NA_character_, length(te))
    pred_f <- prob_to_class_threshold(Pf, pos_class, th_f, tune_threshold)

    prob_m <- if (!is.null(Pm)) extract_prob_pos(Pm, pos_class) else rep(NA_real_, length(te))
    prob_t <- if (!is.null(Pt)) extract_prob_pos(Pt, pos_class) else rep(NA_real_, length(te))
    prob_f <- extract_prob_pos(Pf, pos_class)

    fold_oof <- data.frame(
      repeat_id = r,
      fold = f,
      sample_id = sample_ids[te],
      y_true = as.character(y_te),
      pred_microbiome = pred_m,
      pred_metabolome = pred_t,
      pred_fused = pred_f,
      prob_micro = prob_m,
      prob_metab = prob_t,
      prob_fused = prob_f,
      stringsAsFactors = FALSE
    )
    colnames(fold_oof)[colnames(fold_oof) == "prob_micro"] <- prob_col_names$micro
    colnames(fold_oof)[colnames(fold_oof) == "prob_metab"] <- prob_col_names$metab
    colnames(fold_oof)[colnames(fold_oof) == "prob_fused"] <- prob_col_names$fused
    oof_rows[[length(oof_rows) + 1]] <- fold_oof

    rep_true <- c(rep_true, as.character(y_te))
    rep_pred_micro <- c(rep_pred_micro, pred_m)
    rep_pred_metab <- c(rep_pred_metab, pred_t)
    rep_pred_fused <- c(rep_pred_fused, pred_f)
    rep_prob_micro <- c(rep_prob_micro, prob_m)
    rep_prob_metab <- c(rep_prob_metab, prob_t)
    rep_prob_fused <- c(rep_prob_fused, prob_f)

    threshold_rows[[length(threshold_rows) + 1]] <- data.frame(
      repeat_id = r,
      fold = f,
      block = c("Microbiome", "Metabolome", "Fused"),
      threshold = c(th_m, th_t, th_f),
      stringsAsFactors = FALSE
    )
  }

  # metrics for this repeat
  m_micro <- compute_metrics(rep_true, rep_pred_micro, classes)
  m_metab <- compute_metrics(rep_true, rep_pred_metab, classes)
  m_fused <- compute_metrics(rep_true, rep_pred_fused, classes)

  auc_micro <- compute_auc(rep_true, rep_prob_micro, pos_class)
  auc_metab <- compute_auc(rep_true, rep_prob_metab, pos_class)
  auc_fused <- compute_auc(rep_true, rep_prob_fused, pos_class)

  pr_micro <- compute_pr_auc(rep_true, rep_prob_micro, pos_class)
  pr_metab <- compute_pr_auc(rep_true, rep_prob_metab, pos_class)
  pr_fused <- compute_pr_auc(rep_true, rep_prob_fused, pos_class)

  rep_rows[[length(rep_rows) + 1]] <- data.frame(
    repeat_id = r,
    overall_err_micro = m_micro$overall_err,
    ber_micro = m_micro$ber,
    overall_err_metab = m_metab$overall_err,
    ber_metab = m_metab$ber,
    overall_err_fused = m_fused$overall_err,
    ber_fused = m_fused$ber,
    auc_micro = auc_micro,
    auc_metab = auc_metab,
    auc_fused = auc_fused,
    pr_auc_micro = pr_micro,
    pr_auc_metab = pr_metab,
    pr_auc_fused = pr_fused,
    stringsAsFactors = FALSE
  )

  # pool predictions
  all_true <- c(all_true, rep_true)
  all_preds_micro <- c(all_preds_micro, rep_pred_micro)
  all_preds_metab <- c(all_preds_metab, rep_pred_metab)
  all_preds_fused <- c(all_preds_fused, rep_pred_fused)
  all_prob_micro <- c(all_prob_micro, rep_prob_micro)
  all_prob_metab <- c(all_prob_metab, rep_prob_metab)
  all_prob_fused <- c(all_prob_fused, rep_prob_fused)

  if (r %% 10 == 0) cat("Completed repeat", r, "of", outer_repeats, "\n")
}

rep_df <- do.call(rbind, rep_rows)
write.csv(rep_df, file.path(out_dir, "LF_perf_by_repeat.csv"), row.names = FALSE)

med_iqr <- function(x) {
  c(
    median = median(x, na.rm = TRUE),
    q25 = as.numeric(quantile(x, 0.25, na.rm = TRUE)),
    q75 = as.numeric(quantile(x, 0.75, na.rm = TRUE))
  )
}

make_perf_row <- function(block, ber, auc, pr_auc) {
  balacc <- 1 - ber
  b <- med_iqr(ber)
  ba <- med_iqr(balacc)
  a <- med_iqr(auc)
  p <- med_iqr(pr_auc)
  data.frame(
    block = block,
    ber_median = b["median"], ber_q25 = b["q25"], ber_q75 = b["q75"],
    balacc_median = ba["median"], balacc_q25 = ba["q25"], balacc_q75 = ba["q75"],
    auc_median = a["median"], auc_q25 = a["q25"], auc_q75 = a["q75"],
    pr_auc_median = p["median"], pr_auc_q25 = p["q25"], pr_auc_q75 = p["q75"],
    stringsAsFactors = FALSE
  )
}

perf_compact <- rbind(
  make_perf_row("Microbiome", rep_df$ber_micro, rep_df$auc_micro, rep_df$pr_auc_micro),
  make_perf_row("Metabolome", rep_df$ber_metab, rep_df$auc_metab, rep_df$pr_auc_metab),
  make_perf_row("Fused", rep_df$ber_fused, rep_df$auc_fused, rep_df$pr_auc_fused)
)
write.csv(perf_compact, file.path(out_dir, "LF_perf_compact_table.csv"), row.names = FALSE)

alpha_df <- do.call(rbind, alpha_rows)
write.csv(alpha_df, file.path(out_dir, "LF_alpha_choices.csv"), row.names = FALSE)

threshold_df <- do.call(rbind, threshold_rows)
write.csv(threshold_df, file.path(out_dir, "LF_threshold_choices.csv"), row.names = FALSE)

oof_df <- if (length(oof_rows) > 0) do.call(rbind, oof_rows) else data.frame()
oof_path <- file.path(out_dir, "LF_oof_predictions.csv")
if (nrow(oof_df) > 0) {
  write.csv(oof_df, oof_path, row.names = FALSE)
}

# pooled confusion matrix for fused model
m_pooled <- compute_metrics(all_true, all_preds_fused, classes)
cm <- as.matrix(m_pooled$cm)
write.csv(cm, file.path(out_dir, "LF_confusion_matrix_pooled.csv"))

# pooled confusion matrices for each block
cm_micro <- as.matrix(compute_metrics(all_true, all_preds_micro, classes)$cm)
cm_metab <- as.matrix(compute_metrics(all_true, all_preds_metab, classes)$cm)
write.csv(cm_micro, file.path(out_dir, "LF_confusion_matrix_pooled_microbiome.csv"))
write.csv(cm_metab, file.path(out_dir, "LF_confusion_matrix_pooled_metabolome.csv"))

conf_micro <- compute_binary_stats(all_true, all_preds_micro, pos_class)
conf_metab <- compute_binary_stats(all_true, all_preds_metab, pos_class)
conf_fused <- compute_binary_stats(all_true, all_preds_fused, pos_class)

conf_metrics <- rbind(
  data.frame(block = "Microbiome", conf_micro, stringsAsFactors = FALSE),
  data.frame(block = "Metabolome", conf_metab, stringsAsFactors = FALSE),
  data.frame(block = "Fused", conf_fused, stringsAsFactors = FALSE)
)
write.csv(conf_metrics, file.path(out_dir, "LF_confusion_metrics_pooled.csv"), row.names = FALSE)

roc_paths <- list(metabolome = NA_character_, all = NA_character_)
pr_paths <- list(microbiome = NA_character_, metabolome = NA_character_, fused = NA_character_)
oof_score_plot <- NA_character_
oof_score_plot_fused <- NA_character_
roc_ci_met <- NULL

if (nrow(oof_df) > 0 && length(classes) == 2) {
  class_levels <- classes
  oof_df$y_true_factor <- factor(oof_df$y_true, levels = class_levels)

  get_roc <- function(prob_col) {
    if (!(prob_col %in% colnames(oof_df))) return(NULL)
    vals <- oof_df[[prob_col]]
    if (all(is.na(vals))) return(NULL)
    tryCatch({
      roc(response = oof_df$y_true_factor,
          predictor = vals,
          levels = class_levels,
          direction = "<",
          quiet = TRUE)
    }, error = function(e) NULL)
  }

  roc_micro <- get_roc(prob_col_names$micro)
  roc_met <- get_roc(prob_col_names$metab)
  roc_fused <- get_roc(prob_col_names$fused)

  if (!is.null(roc_met)) {
    auc_met <- auc(roc_met)
    roc_ci_met <- tryCatch(ci.auc(roc_met), error = function(e) NULL)
    ci_txt <- if (!is.null(roc_ci_met)) {
      sprintf("DeLong 95%% CI: %.3f-%.3f", roc_ci_met[1], roc_ci_met[3])
    } else {
      NULL
    }
    p_met <- ggroc(roc_met) +
      labs(
        title = sprintf("ROC (OOF pooled) - Metabolome AUC = %.3f", auc_met),
        subtitle = ci_txt,
        x = "False positive rate (1 - specificity)",
        y = "True positive rate (sensitivity)"
      ) +
      theme_classic()
    roc_paths$metabolome <- file.path(out_dir, "ROC_metabolome_oof.svg")
    ggsave(roc_paths$metabolome, p_met,
           width = roc_plot_dims$width, height = roc_plot_dims$height,
           device = "svg")
  }

  roc_list <- list()
  if (!is.null(roc_micro)) roc_list$Microbiome <- roc_micro
  if (!is.null(roc_met)) roc_list$Metabolome <- roc_met
  if (!is.null(roc_fused)) roc_list$Fused <- roc_fused
  if (length(roc_list) > 0) {
    p_all <- ggroc(roc_list) +
      theme_classic() +
      labs(title = "ROC (OOF pooled)", x = "1 - specificity", y = "Sensitivity")
    roc_paths$all <- file.path(out_dir, "ROC_all_oof.svg")
    ggsave(roc_paths$all, p_all,
           width = roc_plot_dims$width, height = roc_plot_dims$height,
           device = "svg")
  }

  y01 <- ifelse(oof_df$y_true == pos_class, 1, 0)

  prob_micro_vals <- if (prob_col_names$micro %in% colnames(oof_df)) oof_df[[prob_col_names$micro]] else NULL
  if (!is.null(prob_micro_vals)) {
    scores_pos <- prob_micro_vals[y01 == 1]
    scores_neg <- prob_micro_vals[y01 == 0]
    if (length(scores_pos) > 0 && length(scores_neg) > 0) {
      pr_micro <- pr.curve(scores.class0 = scores_pos,
                           scores.class1 = scores_neg,
                           curve = TRUE)
      pr_paths$microbiome <- file.path(out_dir, "PR_microbiome_oof.svg")
      svg(pr_paths$microbiome, width = pr_plot_dims$width, height = pr_plot_dims$height)
      plot(pr_micro)
      dev.off()
    }
  }

  prob_met_vals <- if (prob_col_names$metab %in% colnames(oof_df)) oof_df[[prob_col_names$metab]] else NULL
  if (!is.null(prob_met_vals)) {
    scores_pos <- prob_met_vals[y01 == 1]
    scores_neg <- prob_met_vals[y01 == 0]
    if (length(scores_pos) > 0 && length(scores_neg) > 0) {
      pr_met <- pr.curve(scores.class0 = scores_pos,
                         scores.class1 = scores_neg,
                         curve = TRUE)
      pr_paths$metabolome <- file.path(out_dir, "PR_metabolome_oof.svg")
      svg(pr_paths$metabolome, width = pr_plot_dims$width, height = pr_plot_dims$height)
      plot(pr_met)
      dev.off()
    }
  }

  prob_fused_vals <- if (prob_col_names$fused %in% colnames(oof_df)) oof_df[[prob_col_names$fused]] else NULL
  if (!is.null(prob_fused_vals)) {
    scores_pos <- prob_fused_vals[y01 == 1]
    scores_neg <- prob_fused_vals[y01 == 0]
    if (length(scores_pos) > 0 && length(scores_neg) > 0) {
      pr_fused <- pr.curve(scores.class0 = scores_pos,
                           scores.class1 = scores_neg,
                           curve = TRUE)
      pr_paths$fused <- file.path(out_dir, "PR_fused_oof.svg")
      svg(pr_paths$fused, width = pr_plot_dims$width, height = pr_plot_dims$height)
      plot(pr_fused)
      dev.off()
    }
  }

  if (!is.null(prob_met_vals)) {
    plot_df <- data.frame(
      y_true = factor(oof_df$y_true, levels = class_levels),
      prob_pos = prob_met_vals,
      stringsAsFactors = FALSE
    )
    p_scores <- ggplot(plot_df, aes(x = y_true, y = prob_pos, fill = y_true)) +
      geom_violin(trim = FALSE, alpha = 0.6) +
      geom_boxplot(width = 0.15, outlier.size = 0.6, alpha = 0.8) +
      labs(
        title = sprintf("OOF P(%s) by class (Metabolome)", pos_class),
        x = "True class",
        y = sprintf("Predicted P(%s)", pos_class)
      ) +
      theme_classic() +
      theme(legend.position = "none")

    oof_score_plot <- file.path(out_dir, "OOF_prob_metabolome_violin.svg")
    ggsave(oof_score_plot, p_scores,
           width = pr_plot_dims$width, height = pr_plot_dims$height,
           device = "svg")
  }

  if (!is.null(prob_fused_vals)) {
    plot_df_f <- data.frame(
      y_true = factor(oof_df$y_true, levels = class_levels),
      prob_pos = prob_fused_vals,
      stringsAsFactors = FALSE
    )
    p_scores_f <- ggplot(plot_df_f, aes(x = y_true, y = prob_pos, fill = y_true)) +
      geom_violin(trim = FALSE, alpha = 0.6) +
      geom_boxplot(width = 0.15, outlier.size = 0.6, alpha = 0.8) +
      labs(
        title = sprintf("OOF P(%s) by class (Fused)", pos_class),
        x = "True class",
        y = sprintf("Predicted P(%s)", pos_class)
      ) +
      theme_classic() +
      theme(legend.position = "none")

    oof_score_plot_fused <- file.path(out_dir, "OOF_prob_fused_violin.svg")
    ggsave(oof_score_plot_fused, p_scores_f,
           width = pr_plot_dims$width, height = pr_plot_dims$height,
           device = "svg")
  }
}

perm_summary_rows <- list()
perm_null_store <- list()
perm_summary_path <- file.path(out_dir, "LF_permutation_summary.csv")
perm_results_path <- file.path(out_dir, "LF_permutation_results.rds")

if (n_permutations > 0 && nrow(oof_df) > 0 && length(unique(oof_df$y_true)) >= 2) {
  auc_blocks <- list(
    list(block = "Microbiome", col = prob_col_names$micro),
    list(block = "Metabolome", col = prob_col_names$metab),
    list(block = "Fused", col = prob_col_names$fused)
  )
  for (i in seq_along(auc_blocks)) {
    blk <- auc_blocks[[i]]
    if (!(blk$col %in% colnames(oof_df))) next
    res <- permutation_test_auc(oof_df$y_true, oof_df[[blk$col]], pos_class,
                                n_permutations, perm_seed + i)
    if (is.null(res)) next
    perm_summary_rows[[length(perm_summary_rows) + 1]] <- data.frame(
      block = blk$block,
      metric = "AUC",
      observed = res$observed,
      null_mean = mean(res$null, na.rm = TRUE),
      null_sd = sd(res$null, na.rm = TRUE),
      null_q025 = quantile(res$null, 0.025, na.rm = TRUE),
      null_q975 = quantile(res$null, 0.975, na.rm = TRUE),
      p_value = res$p_value,
      stringsAsFactors = FALSE
    )
    perm_null_store[[paste0(blk$block, "_AUC")]] <- res$null
  }

  pred_blocks <- list(
    list(block = "Microbiome", col = "pred_microbiome"),
    list(block = "Metabolome", col = "pred_metabolome"),
    list(block = "Fused", col = "pred_fused")
  )
  for (j in seq_along(pred_blocks)) {
    blk <- pred_blocks[[j]]
    if (!(blk$col %in% colnames(oof_df))) next
    res <- permutation_test_ber(oof_df$y_true, oof_df[[blk$col]], classes,
                                n_permutations, perm_seed + 100 + j)
    if (is.null(res)) next
    perm_summary_rows[[length(perm_summary_rows) + 1]] <- data.frame(
      block = blk$block,
      metric = "BER",
      observed = res$observed,
      null_mean = mean(res$null, na.rm = TRUE),
      null_sd = sd(res$null, na.rm = TRUE),
      null_q025 = quantile(res$null, 0.025, na.rm = TRUE),
      null_q975 = quantile(res$null, 0.975, na.rm = TRUE),
      p_value = res$p_value,
      stringsAsFactors = FALSE
    )
    perm_null_store[[paste0(blk$block, "_BER")]] <- res$null
  }
}

if (length(perm_summary_rows) > 0) {
  perm_summary_df <- do.call(rbind, perm_summary_rows)
  write.csv(perm_summary_df, perm_summary_path, row.names = FALSE)
  saveRDS(perm_null_store, perm_results_path)
} else {
  perm_summary_path <- NA_character_
  perm_results_path <- NA_character_
}

# feature stability
micro_stab <- data.frame(
  feature = names(micro_feat_counts),
  selected_count = as.integer(micro_feat_counts),
  prop_selected = as.numeric(micro_feat_counts) / max(1L, total_models_fit_micro),
  stringsAsFactors = FALSE
)
micro_beta_summary <- summarize_beta_matrix(micro_beta_mat)
micro_stab <- merge(micro_stab, micro_beta_summary, by = "feature", all.x = TRUE)
micro_stab$mean_beta[is.na(micro_stab$mean_beta)] <- 0
micro_stab$mean_abs_beta[is.na(micro_stab$mean_abs_beta)] <- 0
micro_stab$OR <- exp(micro_stab$mean_beta)
micro_stab <- micro_stab[order(-micro_stab$prop_selected, -micro_stab$mean_abs_beta), , drop = FALSE]
write.csv(micro_stab, file.path(out_dir, "LF_feature_stability_microbiome.csv"), row.names = FALSE)

micro_core <- micro_stab[micro_stab$prop_selected >= core_feature_prop_cutoff, , drop = FALSE]
write.csv(micro_core, file.path(out_dir, "LF_microbiome_core_features.csv"), row.names = FALSE)

metab_stab <- data.frame(
  feature = names(metab_feat_counts),
  selected_count = as.integer(metab_feat_counts),
  prop_selected = as.numeric(metab_feat_counts) / max(1L, total_models_fit_metab),
  stringsAsFactors = FALSE
)
metab_beta_summary <- summarize_beta_matrix(metab_beta_mat)
metab_stab <- merge(metab_stab, metab_beta_summary, by = "feature", all.x = TRUE)
metab_stab$mean_beta[is.na(metab_stab$mean_beta)] <- 0
metab_stab$mean_abs_beta[is.na(metab_stab$mean_abs_beta)] <- 0
metab_stab$OR <- exp(metab_stab$mean_beta)
metab_stab <- metab_stab[order(-metab_stab$prop_selected, -metab_stab$mean_abs_beta), , drop = FALSE]
write.csv(metab_stab, file.path(out_dir, "LF_feature_stability_metabolome.csv"), row.names = FALSE)

metab_core <- metab_stab[metab_stab$prop_selected >= core_feature_prop_cutoff, , drop = FALSE]
write.csv(metab_core, file.path(out_dir, "LF_metabolome_core_features.csv"), row.names = FALSE)

neg_class <- setdiff(classes, pos_class)
if (length(neg_class) != 1) neg_class <- "Other"

make_driver_table <- function(stab_df, block_name, top_n) {
  if (nrow(stab_df) == 0) return(data.frame())
  stab_df$direction <- ifelse(stab_df$mean_beta >= 0,
                              paste0(pos_class, "_up"),
                              paste0(neg_class, "_up"))
  stab_df$block <- block_name
  keep_cols <- c("block", "feature", "prop_selected", "mean_beta", "mean_abs_beta", "OR", "direction")
  out <- stab_df[, keep_cols, drop = FALSE]
  head(out, top_n)
}

drivers_micro <- make_driver_table(micro_stab, "Microbiome", driver_top_n)
drivers_metab <- make_driver_table(metab_stab, "Metabolome", driver_top_n)
drivers_top <- rbind(drivers_micro, drivers_metab)
if (nrow(drivers_top) == 0) {
  drivers_top <- data.frame(
    block = character(0), feature = character(0), prop_selected = numeric(0),
    mean_beta = numeric(0), mean_abs_beta = numeric(0), OR = numeric(0),
    direction = character(0), stringsAsFactors = FALSE
  )
}
write.csv(drivers_top, file.path(out_dir, "LF_driver_features_top.csv"), row.names = FALSE)

# summary text
summary_path <- file.path(out_dir, "LF_perf_summary.txt")
con <- file(summary_path, open = "wt")
writeLines("=== Late fusion (glmnet multinomial elastic net) performance ===", con)
writeLines(paste("outer_folds =", outer_folds,
                 "outer_repeats =", outer_repeats,
                 "inner_folds =", inner_folds,
                 "fusion_mode =", fusion_mode), con)
writeLines("", con)

summ_line <- function(x) {
  sprintf("mean=%.3f | median=%.3f | IQR=%.3f-%.3f", mean(x, na.rm=TRUE), median(x, na.rm=TRUE),
          quantile(x, 0.25, na.rm=TRUE), quantile(x, 0.75, na.rm=TRUE))
}

writeLines("BER (by repeat):", con)
writeLines(paste("  Microbiome:", summ_line(rep_df$ber_micro)), con)
writeLines(paste("  Metabolome:", summ_line(rep_df$ber_metab)), con)
writeLines(paste("  Fused:", summ_line(rep_df$ber_fused)), con)
writeLines("", con)

writeLines("Overall error (by repeat):", con)
writeLines(paste("  Microbiome:", summ_line(rep_df$overall_err_micro)), con)
writeLines(paste("  Metabolome:", summ_line(rep_df$overall_err_metab)), con)
writeLines(paste("  Fused:", summ_line(rep_df$overall_err_fused)), con)
writeLines("", con)

writeLines(paste("Positive class:", pos_class), con)
writeLines(paste("Threshold metric:", threshold_metric), con)
writeLines("AUC (by repeat):", con)
writeLines(paste("  Microbiome:", summ_line(rep_df$auc_micro)), con)
writeLines(paste("  Metabolome:", summ_line(rep_df$auc_metab)), con)
writeLines(paste("  Fused:", summ_line(rep_df$auc_fused)), con)
writeLines("", con)

writeLines("PR-AUC (by repeat):", con)
writeLines(paste("  Microbiome:", summ_line(rep_df$pr_auc_micro)), con)
writeLines(paste("  Metabolome:", summ_line(rep_df$pr_auc_metab)), con)
writeLines(paste("  Fused:", summ_line(rep_df$pr_auc_fused)), con)
writeLines("", con)

if (nrow(oof_df) > 0) {
  writeLines(paste("OOF predictions:", oof_path), con)
}
writeLines(paste("Compact performance table:", file.path(out_dir, "LF_perf_compact_table.csv")), con)
if (!is.na(roc_paths$metabolome)) {
  writeLines(paste("ROC (metabolome) SVG:", roc_paths$metabolome), con)
}
if (!is.na(roc_paths$all)) {
  writeLines(paste("ROC (all blocks) SVG:", roc_paths$all), con)
}
if (!is.na(pr_paths$microbiome)) {
  writeLines(paste("PR (microbiome) SVG:", pr_paths$microbiome), con)
}
if (!is.na(pr_paths$metabolome)) {
  writeLines(paste("PR (metabolome) SVG:", pr_paths$metabolome), con)
}
if (!is.na(pr_paths$fused)) {
  writeLines(paste("PR (fused) SVG:", pr_paths$fused), con)
}
if (!is.na(oof_score_plot)) {
  writeLines(paste("OOF score distribution (metabolome):", oof_score_plot), con)
}
if (!is.na(oof_score_plot_fused)) {
  writeLines(paste("OOF score distribution (fused):", oof_score_plot_fused), con)
}
if (!is.na(perm_summary_path)) {
  writeLines(paste("Permutation summary:", perm_summary_path), con)
}
writeLines(paste("Pooled confusion metrics:", file.path(out_dir, "LF_confusion_metrics_pooled.csv")), con)
writeLines(paste("Driver features (top", driver_top_n, "):", file.path(out_dir, "LF_driver_features_top.csv")), con)
writeLines(paste("Metabolome core features >=", core_feature_prop_cutoff,
                 "prop selected:", nrow(metab_core)), con)
writeLines(paste("Microbiome core features >=", core_feature_prop_cutoff,
                 "prop selected:", nrow(micro_core)), con)
writeLines("", con)

if (mean(rep_df$ber_metab, na.rm = TRUE) <= mean(rep_df$ber_fused, na.rm = TRUE)) {
  writeLines("Fusion did not improve BER vs metabolome; metabolome-only is the primary model.", con)
}

writeLines("Pooled confusion matrix (Fused) saved to LF_confusion_matrix_pooled.csv", con)
writeLines(paste("Total models fit (Microbiome):", total_models_fit_micro), con)
writeLines(paste("Total models fit (Metabolome):", total_models_fit_metab), con)
close(con)

cat("\nWrote outputs to:", out_dir, "\n")
cat("Summary:", summary_path, "\n")

sink()

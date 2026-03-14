#MintTea pipeline from grid search to final analsis

# PREPARE MINTTEA INPUT: MICROBIOME + METABOLOME
# Output table requirements implemented:
# - Rows are samples, columns are features
# - Includes: sample_id, group (healthy/disease)
# - Feature prefixes: T__ (microbiome), M__ (metabolome)

setwd("C:/PhD/Canterbury dataset/DIABLO")
options(stringsAsFactors = FALSE)

# Compatibility helper: allow accidental `minsource()` calls.
if (!exists("minsource", mode = "function")) {
  minsource <- function(file, ...) {
    base::source(file, ...)
  }
}

# USER CONFIG (for species level run use 'species')

microbiome_file <- "Cant_DIABLO_genus_Transposed_30.csv"
metabolome_file <- "metabolome_merged_gcms_primary.csv"
metadata_file   <- "Metadata_good.csv"

metadata_sample_id_col <- NULL
metadata_group_col <- "Group"

# Set TRUE when input matrices are already preprocessed (e.g., filtered, CLR/log transformed, centered/scaled).
input_preprocessed <- TRUE

# Optional strict pairwise analysis filter.
# Set to exactly two metadata group labels (in order: control/healthy, then case/disease)
# Set to NULL to disable pairwise filtering.
pairwise_groups <- NULL

# Optional explicit mapping. Leave NULL to auto-map.
healthy_levels <- NULL
# Example: healthy_levels <- c("Control", "Healthy")

disease_levels <- NULL
# Example: disease_levels <- c("ACeD", "TCeD", "Disease", "Case")

# Rare feature filtering
micro_min_presence_prop <- 0.10   # non-zero in at least 10% samples
micro_min_total_abundance <- 10
metab_max_missing_prop <- 0.20
metab_min_nonmissing <- 5

# Correlated feature clustering
corr_cutoff_micro <- 0.90
corr_cutoff_metab <- 0.90
corr_method <- "spearman"

# Other preprocessing
clr_pseudocount <- 1e-6
zscore_features <- TRUE

# Output
out_dir <- "mixomics_outputs_2block_species_final"
out_file <- file.path(out_dir, "minttea_input_microbiome_metabolome.csv")
out_rds  <- file.path(out_dir, "minttea_input_microbiome_metabolome.rds")

# MintTea grid search parameters
minttea_param_diablo_keepX <- c(10, 20, 30)
minttea_param_sgcca_design <- c(0.3, 0.5, 0.7)
minttea_param_n_repeats <- c(10)
minttea_param_n_folds <- c(5)
minttea_param_sgcca_ncomp <- c(3, 5)
minttea_param_edge_thresholds <- c(0.7, 0.8, 0.9)
minttea_n_evaluation_repeats <- 5
minttea_n_evaluation_folds <- 10
minttea_log_level <- "DEBUG"
minttea_return_main_results_only <- TRUE
minttea_seed <- 27

# Focused export target (single module chosen)
selected_minttea_setting <- "keep_20//des_0.3//nrep_10//nfol_5//ncom_5//edge_0.9"
selected_minttea_modules <- c("module1", "module2", "module3")
export_selected_module_svgs <- TRUE
use_selected_minttea_setting_only <- TRUE

# Repeated external CV for selected setting/modules
run_selected_setting_repeated_external_cv <- TRUE
selected_external_cv_repeats <- 100
selected_external_cv_folds <- 5
selected_external_cv_seed <- 20260306

# Optional: run MintTea and save all generated artifacts 
# Setting 'run_minttea' to TRUE will run the grid search and produce "minttea_settings_module_summary" which can be used to select optimal settings
run_minttea <- FALSE
minttea_dir <- file.path(out_dir, "minttea_outputs")
minttea_result_rds <- file.path(minttea_dir, "minttea_result.rds")
minttea_baseplots_pdf <- file.path(minttea_dir, "minttea_baseplots.pdf")
minttea_module_summary_csv <- file.path(minttea_dir, "minttea_settings_module_summary.csv")
selected_external_cv_results_csv <- file.path(minttea_dir, "selected_setting_repeated_external_cv_results.csv")
selected_external_cv_summary_csv <- file.path(minttea_dir, "selected_setting_repeated_external_cv_summary.csv")
selected_external_cv_predictions_csv <- file.path(minttea_dir, "selected_setting_repeated_external_cv_predictions.csv")
selected_external_cv_roc_points_csv <- file.path(minttea_dir, "selected_setting_repeated_external_cv_roc_points.csv")

# HELPERS

canon_id <- function(x) {
  toupper(gsub("[^A-Za-z0-9]+", "", trimws(as.character(x))))
}

safe_as_numeric <- function(v) {
  v <- trimws(as.character(v))
  v_lower <- tolower(v)

  miss <- is.na(v_lower) | v_lower %in% c(
    "", "na", "n/a", "nan", "null",
    "nd", "n.d.", "not detected", "below detection", "belowdetection",
    "bdl", "trace",
    "#num!", "#div/0!", "#value!", "#ref!"
  )
  v[miss] <- NA
  v <- gsub("^[<>]=?\\s*", "", v)
  v <- gsub(",", "", v)
  v <- gsub("[^0-9eE+\\-\\.]", "", v)

  num <- suppressWarnings(as.numeric(v))
  num[!is.finite(num)] <- NA
  num
}

auto_pick_col <- function(df, candidates) {
  nms <- colnames(df)
  low <- tolower(nms)
  for (c in candidates) {
    hit <- which(low == tolower(c))
    if (length(hit) == 1) return(nms[hit])
  }
  NULL
}

standardize_colnames <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

safe_filename <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if (!nzchar(x)) x <- "artifact"
  x
}

parse_minttea_setting <- function(setting_name) {
  tags <- strsplit(as.character(setting_name), "//", fixed = TRUE)[[1]]
  get_val <- function(prefix) {
    hit <- tags[startsWith(tags, paste0(prefix, "_"))]
    if (length(hit) == 0) return(NA_character_)
    sub(paste0("^", prefix, "_"), "", hit[1])
  }
  list(
    keepX = suppressWarnings(as.numeric(get_val("keep"))),
    design = suppressWarnings(as.numeric(get_val("des"))),
    nrep = suppressWarnings(as.integer(get_val("nrep"))),
    nfol = suppressWarnings(as.integer(get_val("nfol"))),
    ncom = suppressWarnings(as.integer(get_val("ncom"))),
    edge = suppressWarnings(as.numeric(get_val("edge")))
  )
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

compute_roc_points <- function(y_true_binary01, y_score) {
  y <- as.integer(y_true_binary01)
  s <- suppressWarnings(as.numeric(y_score))
  ok <- !is.na(y) & is.finite(s)
  y <- y[ok]
  s <- s[ok]

  n_pos <- sum(y == 1L)
  n_neg <- sum(y == 0L)
  if (length(y) < 2 || n_pos == 0L || n_neg == 0L) {
    return(data.frame())
  }

  thr <- sort(unique(s), decreasing = TRUE)
  thr <- c(Inf, thr, -Inf)

  out <- vector("list", length(thr))
  for (i in seq_along(thr)) {
    pred_pos <- s >= thr[i]
    tp <- sum(pred_pos & y == 1L)
    fp <- sum(pred_pos & y == 0L)
    tpr <- tp / n_pos
    fpr <- fp / n_neg
    out[[i]] <- data.frame(
      threshold = thr[i],
      fpr = fpr,
      tpr = tpr,
      stringsAsFactors = FALSE
    )
  }

  roc <- do.call(rbind, out)
  roc <- roc[order(roc$fpr, roc$tpr, decreasing = FALSE), , drop = FALSE]
  rownames(roc) <- NULL
  roc
}

make_stratified_fold_ids <- function(y, k, seed = NULL) {
  y <- as.character(y)
  n <- length(y)
  out <- integer(n)
  lev <- unique(y)
  if (!is.null(seed)) set.seed(seed)

  for (cl in lev) {
    idx <- which(y == cl)
    idx <- sample(idx, length(idx), replace = FALSE)
    out[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  out
}

run_repeated_external_cv_for_module <- function(proc_data,
                                                feature_names,
                                                repeats = 100,
                                                folds = 5,
                                                seed = 1) {
  if (!all(c("disease_state", feature_names) %in% colnames(proc_data))) {
    stop("proc_data is missing disease_state and/or module features.")
  }

  y_chr <- as.character(proc_data$disease_state)
  y <- ifelse(y_chr == "disease", 1L, 0L)
  if (length(unique(y)) < 2) stop("Need both classes for external CV.")

  min_class_n <- min(table(y))
  k_use <- min(as.integer(folds), as.integer(min_class_n))
  if (!is.finite(k_use) || k_use < 2L) {
    stop("Too few samples in one class for >=2-fold stratified CV.")
  }

  X <- as.matrix(proc_data[, feature_names, drop = FALSE])
  X <- apply(X, 2, function(v) suppressWarnings(as.numeric(v)))
  X <- as.matrix(X)

  rows <- list()
  pred_rows <- list()
  roc_rows <- list()
  rr <- 1L
  pp <- 1L
  qq <- 1L

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
      if (inherits(pca, "try-error")) next
      if (is.null(pca$x) || ncol(pca$x) < 1) next

      train_pc1 <- as.numeric(pca$x[, 1])
      test_pc <- try(stats::predict(pca, newdata = x_test), silent = TRUE)
      if (inherits(test_pc, "try-error") || is.null(test_pc) || ncol(test_pc) < 1) next
      test_pc1 <- as.numeric(test_pc[, 1])

      fit <- try(stats::glm(y_train ~ train_pc1, family = stats::binomial()), silent = TRUE)
      if (inherits(fit, "try-error")) next

      prob_test <- try(stats::predict(fit, newdata = data.frame(train_pc1 = test_pc1), type = "response"), silent = TRUE)
      if (inherits(prob_test, "try-error")) next
      prob_test <- suppressWarnings(as.numeric(prob_test))

      auc_fold <- compute_auc_rank(y_test, prob_test)
      rows[[rr]] <- data.frame(
        repeat_id = rep_i,
        fold_id = fold_i,
        n_train = length(train_idx),
        n_test = length(test_idx),
        auc = auc_fold,
        stringsAsFactors = FALSE
      )

      sample_id_test <- rep(NA_character_, length(test_idx))
      if ("sample_id" %in% colnames(proc_data)) {
        sample_id_test <- as.character(proc_data$sample_id[test_idx])
      }

      pred_rows[[pp]] <- data.frame(
        repeat_id = rep_i,
        fold_id = fold_i,
        sample_index = as.integer(test_idx),
        sample_id = sample_id_test,
        y_true = as.integer(y_test),
        y_score = prob_test,
        stringsAsFactors = FALSE
      )
      pp <- pp + 1L

      roc_df <- compute_roc_points(y_test, prob_test)
      if (nrow(roc_df) > 0) {
        roc_df$repeat_id <- rep_i
        roc_df$fold_id <- fold_i
        roc_rows[[qq]] <- roc_df[, c("repeat_id", "fold_id", "threshold", "fpr", "tpr")]
        qq <- qq + 1L
      }

      rr <- rr + 1L
    }
  }

  if (length(rows) == 0) {
    return(list(metrics = data.frame(), predictions = data.frame(), roc_points = data.frame()))
  }

  list(
    metrics = do.call(rbind, rows),
    predictions = if (length(pred_rows) > 0) do.call(rbind, pred_rows) else data.frame(),
    roc_points = if (length(roc_rows) > 0) do.call(rbind, roc_rows) else data.frame()
  )
}

if (isTRUE(use_selected_minttea_setting_only)) {
  sel <- parse_minttea_setting(selected_minttea_setting)
  if (!all(is.finite(c(sel$keepX, sel$design, sel$edge))) ||
      any(is.na(c(sel$nrep, sel$nfol, sel$ncom)))) {
    stop("selected_minttea_setting is not parseable into MintTea parameters: ", selected_minttea_setting)
  }
  minttea_param_diablo_keepX <- c(sel$keepX)
  minttea_param_sgcca_design <- c(sel$design)
  minttea_param_n_repeats <- c(sel$nrep)
  minttea_param_n_folds <- c(sel$nfol)
  minttea_param_sgcca_ncomp <- c(sel$ncom)
  minttea_param_edge_thresholds <- c(sel$edge)
}

build_minttea_module_summary <- function(minttea_obj) {
  if (is.null(minttea_obj) || !is.list(minttea_obj)) {
    return(data.frame())
  }

  rows <- list()
  r <- 1L

  settings_names <- names(minttea_obj)
  if (is.null(settings_names)) return(data.frame())

  for (settings in settings_names) {
    setting_obj <- minttea_obj[[settings]]
    if (!is.list(setting_obj)) next

    module_names <- names(setting_obj)
    if (is.null(module_names)) next

    module_names <- module_names[grepl("^module", module_names)]
    if (length(module_names) == 0) next

    parsed <- parse_minttea_setting(settings)
    n_modules_setting <- length(module_names)

    for (module in module_names) {
      mod_obj <- setting_obj[[module]]
      if (!is.list(mod_obj)) next

      features <- mod_obj$features
      module_size <- if (is.null(features)) NA_integer_ else length(features)

      auroc <- suppressWarnings(as.numeric(mod_obj$auroc)[1])
      shuf_auc <- suppressWarnings(as.numeric(mod_obj$shuffled_auroc))
      shuf_auc <- shuf_auc[is.finite(shuf_auc)]

      inter_view_corr <- suppressWarnings(as.numeric(mod_obj$inter_view_corr)[1])
      shuf_corr <- suppressWarnings(as.numeric(mod_obj$shuffled_inter_view_corr))
      shuf_corr <- shuf_corr[is.finite(shuf_corr)]

      auc_empirical_p <- if (is.finite(auroc) && length(shuf_auc) > 0) {
        (sum(shuf_auc >= auroc) + 1) / (length(shuf_auc) + 1)
      } else {
        NA_real_
      }

      corr_empirical_p <- if (is.finite(inter_view_corr) && length(shuf_corr) > 0) {
        (sum(shuf_corr >= inter_view_corr) + 1) / (length(shuf_corr) + 1)
      } else {
        NA_real_
      }

      rows[[r]] <- data.frame(
        settings = as.character(settings),
        keepX = parsed$keepX,
        design = parsed$design,
        nrep = parsed$nrep,
        nfol = parsed$nfol,
        ncom = parsed$ncom,
        edge = parsed$edge,
        n_modules_in_setting = n_modules_setting,
        module = as.character(module),
        module_size = module_size,
        auroc = auroc,
        shuffled_auroc_n = length(shuf_auc),
        shuffled_auroc_mean = if (length(shuf_auc) > 0) mean(shuf_auc) else NA_real_,
        shuffled_auroc_sd = if (length(shuf_auc) > 1) sd(shuf_auc) else NA_real_,
        shuffled_auroc_q95 = if (length(shuf_auc) > 0) as.numeric(stats::quantile(shuf_auc, 0.95, na.rm = TRUE)) else NA_real_,
        shuffled_auroc_q99 = if (length(shuf_auc) > 0) as.numeric(stats::quantile(shuf_auc, 0.99, na.rm = TRUE)) else NA_real_,
        auroc_minus_shuffled_mean = if (length(shuf_auc) > 0 && is.finite(auroc)) auroc - mean(shuf_auc) else NA_real_,
        auroc_empirical_p = auc_empirical_p,
        inter_view_corr = inter_view_corr,
        shuffled_inter_view_corr_n = length(shuf_corr),
        shuffled_inter_view_corr_mean = if (length(shuf_corr) > 0) mean(shuf_corr) else NA_real_,
        shuffled_inter_view_corr_sd = if (length(shuf_corr) > 1) sd(shuf_corr) else NA_real_,
        shuffled_inter_view_corr_q95 = if (length(shuf_corr) > 0) as.numeric(stats::quantile(shuf_corr, 0.95, na.rm = TRUE)) else NA_real_,
        shuffled_inter_view_corr_q99 = if (length(shuf_corr) > 0) as.numeric(stats::quantile(shuf_corr, 0.99, na.rm = TRUE)) else NA_real_,
        inter_view_corr_minus_shuffled_mean = if (length(shuf_corr) > 0 && is.finite(inter_view_corr)) inter_view_corr - mean(shuf_corr) else NA_real_,
        inter_view_corr_empirical_p = corr_empirical_p,
        shuffled_inter_view_corr = if (length(shuf_corr) > 0) paste(signif(shuf_corr, 6), collapse = ";") else NA_character_,
        stringsAsFactors = FALSE
      )
      r <- r + 1L
    }
  }

  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

save_minttea_artifacts <- function(obj, output_dir, root_name = "minttea") {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  walk_obj <- function(x, name_path, depth = 0L, max_depth = 3L) {
    if (depth > max_depth || is.null(x)) return(invisible(NULL))
    nm <- safe_filename(name_path)

    if (inherits(x, "ggplot")) {
      if (requireNamespace("ggplot2", quietly = TRUE)) {
        ggplot2::ggsave(
          filename = file.path(output_dir, paste0(nm, ".png")),
          plot = x,
          width = 10,
          height = 7,
          dpi = 300
        )
      }
      return(invisible(NULL))
    }

    if (is.data.frame(x)) {
      write.csv(x, file.path(output_dir, paste0(nm, ".csv")), row.names = FALSE)
      return(invisible(NULL))
    }

    if (is.matrix(x)) {
      write.csv(as.data.frame(x), file.path(output_dir, paste0(nm, ".csv")), row.names = TRUE)
      return(invisible(NULL))
    }

    if (is.list(x) && length(x) > 0) {
      nms <- names(x)
      if (is.null(nms)) nms <- paste0("item", seq_along(x))
      for (i in seq_along(x)) {
        child_name <- paste0(nm, "__", safe_filename(nms[i]))
        walk_obj(x[[i]], child_name, depth = depth + 1L, max_depth = max_depth)
      }
    }
  }

  walk_obj(obj, root_name)
  invisible(NULL)
}

read_metadata <- function(path, sample_col = NULL, group_col = "Group") {
  if (!file.exists(path)) stop("Metadata file not found: ", path)
  md <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)

  if (is.null(sample_col)) {
    sample_col <- auto_pick_col(md, c("sampleid", "sample_id", "sample", "id", "subject", "participant"))
    if (is.null(sample_col)) sample_col <- colnames(md)[1]
  }
  if (!(sample_col %in% colnames(md))) stop("Sample ID column not found in metadata: ", sample_col)

  if (!(group_col %in% colnames(md))) {
    g2 <- auto_pick_col(md, c("group", "condition", "status", "class"))
    if (!is.null(g2)) group_col <- g2
  }
  if (!(group_col %in% colnames(md))) stop("Group column not found in metadata: ", group_col)

  out <- data.frame(
    sample_id_raw = as.character(md[[sample_col]]),
    sample_can = canon_id(md[[sample_col]]),
    group_raw = as.character(md[[group_col]]),
    stringsAsFactors = FALSE
  )

  if (anyDuplicated(out$sample_can)) {
    stop("Duplicate canonical sample IDs in metadata.")
  }

  out
}

map_group_to_binary <- function(group_raw, healthy_levels = NULL, disease_levels = NULL) {
  g <- as.character(group_raw)
  g_trim <- trimws(g)
  g_low <- tolower(g_trim)

  out <- rep(NA_character_, length(g_low))

  if (!is.null(healthy_levels)) {
    out[g_low %in% tolower(healthy_levels)] <- "healthy"
  }
  if (!is.null(disease_levels)) {
    out[g_low %in% tolower(disease_levels)] <- "disease"
  }

  # Direct labels: only fill entries not already mapped by explicit levels.
  out[is.na(out) & g_low %in% c("healthy", "health", "control", "ctrl")] <- "healthy"
  out[is.na(out) & g_low %in% c("disease", "case", "patient", "affected", "aced", "tced", "coeliac", "celiac", "ce") ] <- "disease"

  # If we found at least one healthy level in multi-level data, map all remaining levels to disease
  if (any(out == "healthy", na.rm = TRUE) && any(is.na(out))) {
    out[is.na(out)] <- "disease"
  }

  if (all(is.na(out))) {
    lev <- unique(g_trim)
    if (length(lev) == 2) {
      lev_low <- tolower(lev)
      if (any(grepl("control|healthy|hc|well", lev_low)) && any(grepl("disease|case|patient|aced|tced|coeliac|celiac|ce", lev_low))) {
        healthy_lev <- lev[which(grepl("control|healthy|hc|well", lev_low))[1]]
        out[g_trim == healthy_lev] <- "healthy"
        out[is.na(out)] <- "disease"
      } else {
        out[g_trim == lev[1]] <- "healthy"
        out[g_trim == lev[2]] <- "disease"
        warning("Auto-mapped first metadata level to healthy and second to disease: ", lev[1], " -> healthy, ", lev[2], " -> disease")
      }
    }
  }

  if (any(is.na(out))) {
    bad <- unique(g_trim[is.na(out)])
    stop("Could not map all groups to healthy/disease. Unmapped values: ", paste(bad, collapse = ", "))
  }

  out
}

filter_pairwise_metadata <- function(md_df, pairwise_groups) {
  if (is.null(pairwise_groups)) return(md_df)

  pairwise_groups <- trimws(as.character(pairwise_groups))
  if (length(pairwise_groups) != 2 || any(!nzchar(pairwise_groups))) {
    stop("pairwise_groups must contain exactly two non-empty group labels.")
  }
  if (tolower(pairwise_groups[1]) == tolower(pairwise_groups[2])) {
    stop("pairwise_groups must contain two distinct group labels.")
  }

  md_group_low <- tolower(trimws(as.character(md_df$group_raw)))
  keep <- md_group_low %in% tolower(pairwise_groups)
  md_pair <- md_df[keep, , drop = FALSE]

  if (nrow(md_pair) == 0) {
    stop("No samples matched pairwise_groups: ", paste(pairwise_groups, collapse = " vs "))
  }

  present <- unique(tolower(trimws(as.character(md_pair$group_raw))))
  missing_levels <- pairwise_groups[!(tolower(pairwise_groups) %in% present)]
  if (length(missing_levels) > 0) {
    stop("Pairwise filtering found no samples for group(s): ", paste(missing_levels, collapse = ", "))
  }

  md_pair
}

read_block_oriented <- function(path, sample_ids_can, block_name = "Block") {
  if (!file.exists(path)) stop(block_name, " file not found: ", path)

  # Auto-detect delimiter so CSV/TSV inputs are both supported.
  ext <- tolower(tools::file_ext(path))
  sep_primary <- if (ext %in% c("tsv", "txt")) "\t" else ","
  sep_fallback <- if (sep_primary == "\t") "," else "\t"

  read_delim <- function(sep) {
    utils::read.table(
      path,
      header = TRUE,
      sep = sep,
      quote = "\"",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      fill = TRUE
    )
  }

  df <- tryCatch(read_delim(sep_primary), error = function(e) NULL)
  if (is.null(df) || ncol(df) < 2) {
    df2 <- tryCatch(read_delim(sep_fallback), error = function(e) NULL)
    if (!is.null(df2) && ncol(df2) >= 2) {
      df <- df2
    }
  }

  if (is.null(df) || ncol(df) < 2) {
    stop(block_name, " file should have at least 2 columns after delimiter detection (CSV/TSV).")
  }

  first_col <- as.character(df[[1]])
  first_col_can <- canon_id(first_col)
  col_can <- canon_id(colnames(df)[-1])

  # Option A: first column are sample IDs (rows = samples)
  match_rows <- sum(first_col_can %in% sample_ids_can)
  # Option B: column names (excluding first col) are sample IDs (rows = features)
  match_cols <- sum(col_can %in% sample_ids_can)

  if (match_rows >= match_cols) {
    X_raw <- df[-1]
    rownames(X_raw) <- first_col_can
    X <- as.data.frame(lapply(X_raw, safe_as_numeric), check.names = FALSE)
    rownames(X) <- first_col_can
  } else {
    feature_ids <- as.character(df[[1]])
    mat <- as.data.frame(df[-1], check.names = FALSE)
    mat_num <- as.data.frame(lapply(mat, safe_as_numeric), check.names = FALSE)
    m <- as.matrix(mat_num)
    rownames(m) <- feature_ids
    X <- as.data.frame(t(m), check.names = FALSE)
    rownames(X) <- canon_id(colnames(df)[-1])
  }

  if (anyDuplicated(rownames(X))) stop(block_name, " has duplicated canonical sample IDs.")
  X
}

median_impute <- function(X) {
  X2 <- X
  for (j in seq_len(ncol(X2))) {
    v <- X2[, j]
    if (anyNA(v)) {
      med <- median(v, na.rm = TRUE)
      if (!is.finite(med)) med <- 0
      v[is.na(v)] <- med
      X2[, j] <- v
    }
  }
  X2
}

zscore_matrix <- function(X) {
  X2 <- as.matrix(scale(X))
  if (is.null(dim(X2))) {
    X2 <- matrix(X2, ncol = 1)
    colnames(X2) <- colnames(X)
    rownames(X2) <- rownames(X)
  }
  bad_cols <- apply(X2, 2, function(v) all(is.na(v) | is.infinite(v)))
  if (any(bad_cols)) {
    X2[, bad_cols] <- 0
  }
  X2[!is.finite(X2)] <- 0
  X2
}

cluster_correlated_features <- function(X, prefix, cutoff = 0.90, method = "spearman") {
  if (ncol(X) <= 1) return(X)

  cor_mat <- suppressWarnings(cor(X, use = "pairwise.complete.obs", method = method))
  cor_mat[!is.finite(cor_mat)] <- 0
  diag(cor_mat) <- 1

  d <- as.dist(1 - abs(cor_mat))
  hc <- hclust(d, method = "average")
  cl <- cutree(hc, h = 1 - cutoff)

  grouped <- split(seq_len(ncol(X)), cl)
  X_new <- matrix(NA_real_, nrow = nrow(X), ncol = length(grouped))
  rownames(X_new) <- rownames(X)

  new_names <- character(length(grouped))
  k <- 1
  for (gid in names(grouped)) {
    idx <- grouped[[gid]]
    if (length(idx) == 1) {
      X_new[, k] <- X[, idx]
      nm <- colnames(X)[idx]
      if (!grepl(paste0("^", prefix, "__"), nm)) {
        nm <- paste0(prefix, "__", nm)
      }
      new_names[k] <- nm
    } else {
      X_new[, k] <- rowMeans(X[, idx, drop = FALSE], na.rm = TRUE)
      new_names[k] <- paste0(prefix, "__CL", sprintf("%04d", as.integer(gid)), "_n", length(idx))
    }
    k <- k + 1
  }

  colnames(X_new) <- make.unique(new_names)
  X_new
}


# MAIN

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (isTRUE(input_preprocessed)) {
  cat("Input preprocessed mode ON: disabling additional filtering/transforms/scaling/clustering.\n")
}

micro_min_presence_prop_use <- if (isTRUE(input_preprocessed)) 0 else micro_min_presence_prop
micro_min_total_abundance_use <- if (isTRUE(input_preprocessed)) 0 else micro_min_total_abundance
metab_max_missing_prop_use <- if (isTRUE(input_preprocessed)) 1 else metab_max_missing_prop
metab_min_nonmissing_use <- if (isTRUE(input_preprocessed)) 1 else metab_min_nonmissing
zscore_features_use <- if (isTRUE(input_preprocessed)) FALSE else zscore_features

apply_micro_relative_clr <- !isTRUE(input_preprocessed)
apply_metab_log1p <- !isTRUE(input_preprocessed)
apply_corr_clustering <- TRUE

md <- read_metadata(metadata_file, metadata_sample_id_col, metadata_group_col)
sample_ids_can <- md$sample_can

micro_raw <- read_block_oriented(microbiome_file, sample_ids_can, "Microbiome")
metab_raw <- read_block_oriented(metabolome_file, sample_ids_can, "Metabolome")

ids <- Reduce(base::intersect, list(sample_ids_can, rownames(micro_raw), rownames(metab_raw)))
if (length(ids) < 6) stop("Too few matched samples after alignment.")

md2 <- md[match(ids, md$sample_can), , drop = FALSE]
md2 <- filter_pairwise_metadata(md2, pairwise_groups)
ids <- md2$sample_can

if (length(ids) < 6) stop("Too few samples left after pairwise group filtering.")

healthy_levels_use <- healthy_levels
disease_levels_use <- disease_levels
if (!is.null(pairwise_groups)) {
  healthy_levels_use <- pairwise_groups[1]
  disease_levels_use <- pairwise_groups[2]
  cat("Pairwise mode ON:", pairwise_groups[1], "vs", pairwise_groups[2], "\n")
  cat("Pairwise metadata counts:\n")
  print(table(md2$group_raw))
}

group_bin <- map_group_to_binary(md2$group_raw, healthy_levels_use, disease_levels_use)

if (length(unique(group_bin)) < 2) {
  stop(
    "Group mapping collapsed to one class. Check pairwise_groups / healthy_levels / disease_levels. ",
    "Mapped counts: ",
    paste(names(table(group_bin)), as.integer(table(group_bin)), collapse = ", ")
  )
}

# Align
micro <- as.matrix(micro_raw[ids, , drop = FALSE])
metab <- as.matrix(metab_raw[ids, , drop = FALSE])


# PREPROCESS MICROBIOME (T__)
# 1) Rare feature filtering
present_prop <- colMeans(micro > 0, na.rm = TRUE)
total_abund <- colSums(micro, na.rm = TRUE)
keep_micro <- (present_prop >= micro_min_presence_prop_use) & (total_abund >= micro_min_total_abundance_use)
micro <- micro[, keep_micro, drop = FALSE]

if (ncol(micro) < 2) stop("Microbiome filtering left too few features.")

# 2) Normalize + CLR when non-negative (raw-input mode only)
if (isTRUE(apply_micro_relative_clr) && all(micro >= 0, na.rm = TRUE)) {
  rs <- rowSums(micro, na.rm = TRUE)
  rs[!is.finite(rs) | rs <= 0] <- 1
  micro <- micro / rs
  micro <- log(micro + clr_pseudocount)
  micro <- micro - rowMeans(micro, na.rm = TRUE)
} else if (isTRUE(apply_micro_relative_clr)) {
  warning("Microbiome matrix contains negatives; skipping relative-abundance CLR and using median-impute + scale.")
}

micro <- median_impute(micro)
if (isTRUE(zscore_features_use)) micro <- zscore_matrix(micro)

colnames(micro) <- paste0("T__", standardize_colnames(colnames(micro)))
if (isTRUE(apply_corr_clustering)) {
  micro <- cluster_correlated_features(micro, prefix = "T", cutoff = corr_cutoff_micro, method = corr_method)
}


# PREPROCESS METABOLOME (M__)
# 1) Rare feature filtering (missingness-based)
na_prop <- colMeans(is.na(metab))
nonmiss <- colSums(!is.na(metab))
keep_metab <- (na_prop <= metab_max_missing_prop_use) & (nonmiss >= metab_min_nonmissing_use)
metab <- metab[, keep_metab, drop = FALSE]

if (ncol(metab) < 2) stop("Metabolome filtering left too few features.")

# 2) Impute + optional log1p for non-negative data + scale
metab <- median_impute(metab)
if (isTRUE(apply_metab_log1p) && all(metab >= 0, na.rm = TRUE)) {
  metab <- log1p(metab)
}
if (isTRUE(zscore_features_use)) metab <- zscore_matrix(metab)

colnames(metab) <- paste0("M__", standardize_colnames(colnames(metab)))
if (isTRUE(apply_corr_clustering)) {
  metab <- cluster_correlated_features(metab, prefix = "M", cutoff = corr_cutoff_metab, method = corr_method)
}


# BUILD FINAL MINTTEA INPUT TABLE
out_df <- data.frame(
  sample_id = md2$sample_id_raw,
  micro,
  metab,
  disease_state = group_bin,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

write.csv(out_df, out_file, row.names = FALSE)
saveRDS(out_df, out_rds)

cat("=== MintTea input table created ===\n")
cat("Samples:", nrow(out_df), "\n")
cat("Microbiome features after preprocessing:", ncol(micro), "\n")
cat("Metabolome features after preprocessing:", ncol(metab), "\n")
cat("Output CSV:", out_file, "\n")
cat("Output RDS:", out_rds, "\n")
cat("Group counts:\n")
print(table(out_df$disease_state))

if (isTRUE(run_minttea)) {
  dir.create(minttea_dir, showWarnings = FALSE, recursive = TRUE)
  if (!requireNamespace("MintTea", quietly = TRUE)) {
    warning("MintTea package is not installed; skipping MintTea run.")
  } else {
    minttea_results <- MintTea::MintTea(
      proc_data = out_df,
      study_group_column = "disease_state",
      control_group_name = "healthy",
      case_group_name = "disease",
      sample_id_column = "sample_id",
      view_prefixes = c("T", "M"),
      param_diablo_keepX = minttea_param_diablo_keepX,
      param_sgcca_design = minttea_param_sgcca_design,
      param_n_repeats = minttea_param_n_repeats,
      param_n_folds = minttea_param_n_folds,
      param_sgcca_ncomp = minttea_param_sgcca_ncomp,
      param_edge_thresholds = minttea_param_edge_thresholds,
      n_evaluation_repeats = minttea_n_evaluation_repeats,
      n_evaluation_folds = minttea_n_evaluation_folds,
      log_level = minttea_log_level,
      return_main_results_only = minttea_return_main_results_only,
      seed = minttea_seed
    )

    # Keep both object names for downstream compatibility.
    minttea_result <- minttea_results
    saveRDS(minttea_results, minttea_result_rds)

    minttea_module_summary <- build_minttea_module_summary(minttea_results)
    if (nrow(minttea_module_summary) > 0) {
      utils::write.csv(minttea_module_summary, minttea_module_summary_csv, row.names = FALSE)
      cat("MintTea setting/module summary CSV:", minttea_module_summary_csv, "\n")
    } else {
      warning("Could not build MintTea module summary table (no module rows detected).")
    }

    save_minttea_artifacts(minttea_results, minttea_dir, root_name = "minttea")
    cat("MintTea grid search completed and saved to:", minttea_dir, "\n")
    cat("MintTea result RDS:", minttea_result_rds, "\n")
  }
}


# OPTIONAL: RUN MINTTEA + SAVE ARTIFACTS


plot_module_with_igraph <- function(minttea_module, prefix_colors = c('T' = 'darkred', 'M' = 'purple4', 'P' = 'darkgreen')) {
  require(igraph)
  nodes <- data.frame(node = minttea_module$features, stringsAsFactors = FALSE)
  nodes$prefix <- gsub('__.*$','', nodes$node)
  nodes$color <- unname(prefix_colors[nodes$prefix])
  nodes$color[is.na(nodes$color)] <- 'grey40'

  edges <- minttea_module$module_edges
  if (is.data.frame(edges) && ('edge_weight' %in% colnames(edges))) {
    edges$width <- edges$edge_weight
  }
  g <- graph_from_data_frame(edges, directed=F, vertices = nodes)
  plot(g, vertex.label.cex = 0.8, vertex.label.dist = -2, vertex.label.color = 'black', edge.color = 'grey80')
}

plot_module_stats <- function(minttea_results, settings, module) {
  require(dplyr)
  require(ggplot2)

  tmp <- data.frame(AUC = minttea_results[[settings]][[module]]$shuffled_auroc, Mode = 'Null')
  tmp$Mode <- factor(tmp$Mode, levels = c('Null', 'True module'))
  true_auc <- minttea_results[[settings]][[module]]$auroc
  p1 <- ggplot(tmp, aes(fill = Mode, x = AUC)) +
    geom_histogram(color = 'black', bins = 20) +
    scale_fill_manual(values = c('Null' = 'grey80', 'True module' = 'orangered3'), drop = FALSE) +
    scale_y_continuous(expand = c(0,0)) +
    geom_vline(color = 'orangered3', xintercept = true_auc, linewidth = 2) +
    ylab('Count') +
    xlab("AUC of module's 1st PC") +
    theme_classic() +
    theme(legend.title = element_blank())

  tmp <- data.frame(Cross_Omic_Correlation = minttea_results[[settings]][[module]]$shuffled_inter_view_corr, Mode = 'Null')
  tmp$Mode <- factor(tmp$Mode, levels = c('Null', 'True module'))
  true_cor <- minttea_results[[settings]][[module]]$inter_view_corr
  p2 <- ggplot(tmp, aes(fill = Mode, x = Cross_Omic_Correlation)) +
    geom_histogram(color = 'black', bins = 20) +
    scale_fill_manual(values = c('Null' = 'grey80', 'True module' = 'orangered3'), drop = FALSE) +
    scale_y_continuous(expand = c(0,0)) +
    geom_vline(color = 'orangered3', xintercept = true_cor, linewidth = 2) +
    ylab('Count') +
    xlab('Avg. correlation between\nfeatures cross-omic') +
    theme_classic() +
    theme(legend.title = element_blank())

  return(list(p1 = p1, p2 = p2))
}

extract_overview_tables <- function(setting_obj) {
  if (!is.list(setting_obj)) return(NULL)

  if (all(c("sens_analysis_modules", "modules_overview") %in% names(setting_obj))) {
    return(list(
      sens_analysis_modules = setting_obj$sens_analysis_modules,
      modules_overview = setting_obj$modules_overview
    ))
  }

  for (nm in names(setting_obj)) {
    x <- setting_obj[[nm]]
    if (is.list(x) && all(c("sens_analysis_modules", "modules_overview") %in% names(x))) {
      return(list(
        sens_analysis_modules = x$sens_analysis_modules,
        modules_overview = x$modules_overview
      ))
    }
  }

  NULL
}

plot_modules_overview <- function(sens_analysis_modules,
                                  modules_overview,
                                  feature_type_color_map = c(T = "darkred", P = "darkgreen", M = "purple4"),
                                  hide_y_axis_text = FALSE,
                                  show_rf = FALSE,
                                  mean_auc_rf = NA_real_) {
  if (!requireNamespace("dplyr", quietly = TRUE) || !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_modules_overview requires dplyr and ggplot2.")
  }

  dataset_order <- unique(as.character(modules_overview$dataset))

  tmp1 <- sens_analysis_modules |>
    dplyr::mutate(feature_type = substr(feature, 1, 1)) |>
    dplyr::mutate(module_num = as.numeric(gsub("module", "", module))) |>
    dplyr::group_by(dataset, module_num, feature_type) |>
    dplyr::summarise(N = dplyr::n(), .groups = "drop")

  modules_to_plot <- modules_overview |>
    dplyr::mutate(module_num = as.numeric(gsub("module", "", module))) |>
    dplyr::filter(.data$multi_view) |>
    dplyr::select(dataset, module_num, is_interesting)

  tmp1 <- dplyr::inner_join(tmp1, modules_to_plot, by = c("dataset", "module_num"))
  if (nrow(tmp1) == 0) stop("No multi-view modules available for overview plot.")

  tmp1$module2 <- factor(paste("Module", tmp1$module_num), levels = paste("Module", sort(unique(tmp1$module_num), decreasing = TRUE)))

  p1 <- ggplot2::ggplot(
    tmp1 |>
      dplyr::filter(.data$dataset %in% dataset_order) |>
      dplyr::mutate(dataset = factor(.data$dataset, levels = dataset_order)),
    ggplot2::aes(x = module2, y = N, fill = feature_type)
  ) +
    ggplot2::geom_bar(ggplot2::aes(alpha = is_interesting, color = is_interesting), width = 0.8, stat = "identity") +
    ggplot2::scale_alpha_manual(values = c("FALSE" = 0.3, "TRUE" = 0.9), guide = "none") +
    ggplot2::scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "black"), guide = "none") +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) +
    ggplot2::scale_x_discrete(expand = c(0, 0.5)) +
    ggplot2::coord_flip() +
    ggplot2::theme_classic() +
    ggplot2::facet_grid(rows = ggplot2::vars(dataset), space = "free_y", scales = "free_y", switch = "y") +
    ggplot2::xlab(NULL) +
    ggplot2::ylab("No. of features\nof each type") +
    ggplot2::scale_fill_manual(name = "View", values = feature_type_color_map) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_line(linewidth = 0.5, color = "grey93")) +
    ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank()) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_line(linewidth = 0.5, color = "grey93")) +
    ggplot2::theme(strip.background = ggplot2::element_blank()) +
    ggplot2::theme(strip.placement = "outside")

  if (isTRUE(hide_y_axis_text)) p1 <- p1 + ggplot2::theme(axis.text.y = ggplot2::element_blank())

  tmp2 <- modules_overview |>
    dplyr::mutate(module_num = as.numeric(gsub("module", "", module))) |>
    dplyr::inner_join(modules_to_plot, by = c("dataset", "module_num"), suffix = c("", ".sel")) |>
    dplyr::mutate(is_interesting = ifelse(!is.na(.data$is_interesting.sel), .data$is_interesting.sel, .data$is_interesting))

  tmp2$module2 <- factor(paste("Module", tmp2$module_num), levels = levels(tmp1$module2))
  points_size <- ifelse(isTRUE(hide_y_axis_text), 3, 3.7)

  p2 <- ggplot2::ggplot(
    tmp2 |>
      dplyr::filter(.data$dataset %in% dataset_order) |>
      dplyr::mutate(dataset = factor(.data$dataset, levels = dataset_order)),
    ggplot2::aes(x = module2)
  ) +
    ggplot2::geom_hline(yintercept = 0.5, color = "darkred", linetype = "dashed", linewidth = 1) +
    ggplot2::geom_linerange(ggplot2::aes(ymax = mean_module_auc_shuffled + sd_module_auc_shuffled, ymin = mean_module_auc_shuffled - sd_module_auc_shuffled), alpha = 0.4, linewidth = 2, color = "grey70") +
    ggplot2::geom_point(ggplot2::aes(y = mean_module_auc_shuffled), shape = 16, size = points_size - 0.5, color = "grey60", alpha = 0.8) +
    ggplot2::geom_point(ggplot2::aes(y = mean_module_auc, fill = is_interesting, color = is_interesting), shape = 23, size = points_size, alpha = 0.9) +
    ggplot2::scale_fill_manual(values = c("FALSE" = "#D7E7ED", "TRUE" = "skyblue4"), guide = "none") +
    ggplot2::scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "black"), guide = "none") +
    ggplot2::scale_y_continuous(breaks = seq(0.5, 1, 0.1), expand = ggplot2::expansion(mult = c(0.05, 0.1))) +
    ggplot2::scale_x_discrete(expand = c(0, 0.5)) +
    ggplot2::coord_flip() +
    ggplot2::theme_classic() +
    ggplot2::xlab(NULL) +
    ggplot2::ylab("Module AUC") +
    ggplot2::facet_grid(rows = ggplot2::vars(dataset), space = "free_y", scales = "free_y", switch = "y") +
    ggplot2::theme(axis.text.y = ggplot2::element_blank()) +
    ggplot2::theme(strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_blank())

  if (isTRUE(show_rf) && is.finite(mean_auc_rf)) {
    p2 <- p2 + ggplot2::geom_hline(yintercept = mean_auc_rf, color = "goldenrod2", linewidth = 2, alpha = 0.7)
  }

  p3 <- ggplot2::ggplot(
    tmp2 |>
      dplyr::filter(.data$dataset %in% dataset_order) |>
      dplyr::mutate(dataset = factor(.data$dataset, levels = dataset_order)),
    ggplot2::aes(x = module2)
  ) +
    ggplot2::geom_linerange(ggplot2::aes(ymax = avg_spear_corr_shuffled + sd_spear_corr_shuffled, ymin = avg_spear_corr_shuffled - sd_spear_corr_shuffled), alpha = 0.4, linewidth = 2, color = "grey70") +
    ggplot2::geom_point(ggplot2::aes(y = avg_spear_corr_shuffled), shape = 16, size = points_size - 0.5, color = "grey60", alpha = 0.8) +
    ggplot2::geom_point(ggplot2::aes(y = avg_spear_corr, fill = is_interesting, color = is_interesting), shape = 23, size = points_size, alpha = 0.9) +
    ggplot2::scale_fill_manual(values = c("FALSE" = "#D9B9AB", "TRUE" = "sienna4"), guide = "none") +
    ggplot2::scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "black"), guide = "none") +
    ggplot2::scale_x_discrete(expand = c(0, 0.5)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.1))) +
    ggplot2::coord_flip() +
    ggplot2::theme_classic() +
    ggplot2::xlab(NULL) +
    ggplot2::ylab("Cross-view avg.\ncorrelation") +
    ggplot2::facet_grid(rows = ggplot2::vars(dataset), space = "free_y", scales = "free_y", switch = "y") +
    ggplot2::theme(axis.text.y = ggplot2::element_blank()) +
    ggplot2::theme(strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_blank())

  list(p1 = p1, p2 = p2, p3 = p3)
}


# SAVE CUSTOM MODULE PLOTS 

minttea_obj <- NULL
if (exists("minttea_result") && !is.null(minttea_result)) minttea_obj <- minttea_result
if (is.null(minttea_obj) && exists("minttea_results") && !is.null(minttea_results)) minttea_obj <- minttea_results

if (!is.null(minttea_obj)) {
  custom_plot_dir <- file.path(minttea_dir, "custom_module_plots")
  dir.create(custom_plot_dir, showWarnings = FALSE, recursive = TRUE)

  settings_names <- names(minttea_obj)
  if (length(settings_names) > 0) {
    for (settings in settings_names) {
      setting_obj <- minttea_obj[[settings]]
      if (!is.list(setting_obj)) next

      module_names <- names(setting_obj)
      if (is.null(module_names) || length(module_names) == 0) next

      for (module in module_names) {
        mod_obj <- setting_obj[[module]]
        if (!is.list(mod_obj)) next
        settings_tag <- safe_filename(settings)
        module_tag <- safe_filename(module)

        # Save igraph module network if structure is present
        if (all(c("features", "module_edges") %in% names(mod_obj))) {
          network_png <- file.path(custom_plot_dir, paste0("network_", settings_tag, "_", module_tag, ".png"))
          try({
            grDevices::png(network_png, width = 2400, height = 1800, res = 300)
            plot_module_with_igraph(mod_obj)
            grDevices::dev.off()
          }, silent = TRUE)
        }

        # Save module stats plots if AUROC/correlation shuffle outputs exist
        if (all(c("shuffled_auroc", "auroc", "shuffled_inter_view_corr", "inter_view_corr") %in% names(mod_obj))) {
          stats_plots <- try(plot_module_stats(minttea_obj, settings, module), silent = TRUE)
          if (!inherits(stats_plots, "try-error") && requireNamespace("ggplot2", quietly = TRUE)) {
            p1_file <- file.path(custom_plot_dir, paste0("stats_auc_", settings_tag, "_", module_tag, ".png"))
            p2_file <- file.path(custom_plot_dir, paste0("stats_cor_", settings_tag, "_", module_tag, ".png"))
            try(ggplot2::ggsave(filename = p1_file, plot = stats_plots$p1, width = 8, height = 6, dpi = 300), silent = TRUE)
            try(ggplot2::ggsave(filename = p2_file, plot = stats_plots$p2, width = 8, height = 6, dpi = 300), silent = TRUE)
          }
        }
      }

      # Save overview strips for this setting if overview tables are available
      overview_tables <- extract_overview_tables(setting_obj)
      if (!is.null(overview_tables)) {
        settings_tag <- safe_filename(settings)
        overview_plots <- try(
          plot_modules_overview(
            sens_analysis_modules = overview_tables$sens_analysis_modules,
            modules_overview = overview_tables$modules_overview
          ),
          silent = TRUE
        )

        if (!inherits(overview_plots, "try-error") && requireNamespace("ggplot2", quietly = TRUE)) {
          try(ggplot2::ggsave(
            filename = file.path(custom_plot_dir, paste0("overview_strip1_features_", settings_tag, ".png")),
            plot = overview_plots$p1,
            width = 10,
            height = 8,
            dpi = 300
          ), silent = TRUE)

          try(ggplot2::ggsave(
            filename = file.path(custom_plot_dir, paste0("overview_strip2_cor_", settings_tag, ".png")),
            plot = overview_plots$p3,
            width = 8,
            height = 8,
            dpi = 300
          ), silent = TRUE)

          try(ggplot2::ggsave(
            filename = file.path(custom_plot_dir, paste0("overview_strip3_auc_", settings_tag, ".png")),
            plot = overview_plots$p2,
            width = 8,
            height = 8,
            dpi = 300
          ), silent = TRUE)

          if (requireNamespace("cowplot", quietly = TRUE)) {
            combined <- cowplot::plot_grid(overview_plots$p1, overview_plots$p3, overview_plots$p2, nrow = 1, rel_widths = c(8, 3, 3), align = "h", axis = "tb")
            try(ggplot2::ggsave(
              filename = file.path(custom_plot_dir, paste0("overview_combined_", settings_tag, ".png")),
              plot = combined,
              width = 20,
              height = 8,
              dpi = 300
            ), silent = TRUE)
          }
        }
      }
    }
  }

  cat("Custom module plots saved to:", custom_plot_dir, "\n")
}

resolve_existing_minttea_dir <- function(primary_dir, study_group_col = "disease_state") {
  candidates <- unique(c(
    primary_dir,
    paste0(primary_dir, "_", study_group_col),
    paste0(primary_dir, "_disease_status"),
    file.path(dirname(primary_dir), "minttea_outputs_disease_status"),
    file.path(dirname(primary_dir), "minttea_outputs_disease_state")
  ))

  hit <- candidates[file.exists(file.path(candidates, "minttea_result.rds"))]
  if (length(hit) > 0) return(hit[1])
  primary_dir
}

load_minttea_obj_if_needed <- function(existing_obj, minttea_output_dir) {
  if (!is.null(existing_obj)) return(existing_obj)
  rds_path <- file.path(minttea_output_dir, "minttea_result.rds")
  if (file.exists(rds_path)) {
    return(readRDS(rds_path))
  }
  NULL
}

minttea_result_rds_exists <- function(minttea_output_dir) {
  file.exists(file.path(minttea_output_dir, "minttea_result.rds"))
}

export_cytoscape_graphml <- function(node_table, edge_table, out_graphml) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    warning("igraph is required for GraphML export.")
    return(invisible(FALSE))
  }
  if (is.null(node_table) || !is.data.frame(node_table) || nrow(node_table) == 0) {
    warning("No node table available for GraphML export.")
    return(invisible(FALSE))
  }
  if (is.null(edge_table) || !is.data.frame(edge_table) || nrow(edge_table) == 0) {
    warning("No edge table available for GraphML export.")
    return(invisible(FALSE))
  }

  vdf <- node_table
  if (!("id" %in% colnames(vdf))) {
    if ("node" %in% colnames(vdf)) {
      vdf$id <- as.character(vdf$node)
    } else {
      warning("Node table needs 'id' or 'node' column for GraphML export.")
      return(invisible(FALSE))
    }
  }
  vdf$id <- as.character(vdf$id)

  edf <- edge_table
  if (!all(c("source", "target") %in% colnames(edf))) {
    if (all(c("feature.x", "feature.y") %in% colnames(edf))) {
      edf$source <- as.character(edf$feature.x)
      edf$target <- as.character(edf$feature.y)
    } else {
      warning("Edge table needs source/target (or feature.x/feature.y) for GraphML export.")
      return(invisible(FALSE))
    }
  }
  edf$source <- as.character(edf$source)
  edf$target <- as.character(edf$target)

 
  edge_nodes <- unique(c(edf$source, edf$target))
  missing_nodes <- setdiff(edge_nodes, vdf$id)
  if (length(missing_nodes) > 0) {
    pad <- data.frame(
      id = missing_nodes,
      feature_type = NA_character_,
      degree = NA_real_,
      stringsAsFactors = FALSE
    )
    vdf <- rbind(vdf, pad)
  }

  # igraph uses first two columns as from/to; force correct ordering.
  edge_attr <- edf[, setdiff(colnames(edf), c("source", "target", "feature.x", "feature.y")), drop = FALSE]
  edge_core <- data.frame(from = edf$source, to = edf$target, stringsAsFactors = FALSE)
  d_use <- cbind(edge_core, edge_attr)

  vertex_attr <- vdf[, setdiff(colnames(vdf), "id"), drop = FALSE]
  g <- igraph::graph_from_data_frame(
    d = d_use,
    directed = TRUE,
    vertices = data.frame(name = vdf$id, vertex_attr, check.names = FALSE)
  )

  dir.create(dirname(out_graphml), recursive = TRUE, showWarnings = FALSE)
  igraph::write_graph(g, file = out_graphml, format = "graphml")
  invisible(TRUE)
}

plot_module_with_degree_and_signed_edges_svg <- function(
    minttea_module,
    out_svg,
    proc_data = NULL,
    corr_method = "spearman",
    edge_color_mode = c("weight", "signed_corr"),
    prefix_colors = c("T" = "darkred", "M" = "purple4", "P" = "darkgreen"),
    weight_edge_color = "#3a7ca5",
    pos_edge_color = "#1f77b4",
    neg_edge_color = "#d62728",
    na_corr_color = "#9e9e9e") {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph is required for network plotting.")
  }
  edge_color_mode <- match.arg(edge_color_mode)

  edges <- minttea_module$module_edges
  if (!is.data.frame(edges) || !all(c("feature.x", "feature.y") %in% colnames(edges))) {
    stop("Module does not contain a valid module_edges table.")
  }

  edges$feature.x <- as.character(edges$feature.x)
  edges$feature.y <- as.character(edges$feature.y)
  edges$edge_weight <- suppressWarnings(as.numeric(edges$edge_weight))
  edges$edge_weight_sign <- ifelse(is.finite(edges$edge_weight) & edges$edge_weight < 0, "negative", "positive")

  edges$signed_corr <- NA_real_
  if (!is.null(proc_data) && is.data.frame(proc_data) && edge_color_mode == "signed_corr") {
    for (i in seq_len(nrow(edges))) {
      fx <- edges$feature.x[i]
      fy <- edges$feature.y[i]
      if (fx %in% colnames(proc_data) && fy %in% colnames(proc_data)) {
        vx <- suppressWarnings(as.numeric(proc_data[[fx]]))
        vy <- suppressWarnings(as.numeric(proc_data[[fy]]))
        ok <- is.finite(vx) & is.finite(vy)
        if (sum(ok) >= 3) {
          edges$signed_corr[i] <- suppressWarnings(stats::cor(vx[ok], vy[ok], method = corr_method))
        }
      }
    }
  }
  edges$corr_sign <- ifelse(
    is.finite(edges$signed_corr) & edges$signed_corr < 0,
    "negative",
    ifelse(is.finite(edges$signed_corr), "positive", "unknown")
  )

  nodes <- sort(unique(c(edges$feature.x, edges$feature.y)))
  node_df <- data.frame(node = nodes, stringsAsFactors = FALSE)
  node_df$prefix <- sub("__.*$", "", node_df$node)
  node_df$color <- unname(prefix_colors[node_df$prefix])
  node_df$color[is.na(node_df$color)] <- "grey40"

  g <- igraph::graph_from_data_frame(
    d = edges[, c("feature.x", "feature.y")],
    directed = FALSE,
    vertices = node_df
  )

  deg <- igraph::degree(g)
  igraph::V(g)$degree <- as.numeric(deg)
  igraph::V(g)$size <- 6 + 2 * sqrt(pmax(igraph::V(g)$degree, 0))
  igraph::V(g)$color <- node_df$color[match(igraph::V(g)$name, node_df$node)]
  igraph::V(g)$frame.color <- "grey30"
  igraph::V(g)$label.color <- "black"
  igraph::V(g)$label.cex <- 0.75

  if (edge_color_mode == "signed_corr") {
    igraph::E(g)$color <- ifelse(
      edges$corr_sign == "negative",
      neg_edge_color,
      ifelse(edges$corr_sign == "positive", pos_edge_color, na_corr_color)
    )
  } else {
    igraph::E(g)$color <- weight_edge_color
  }
  abs_w <- abs(edges$edge_weight)
  abs_w[!is.finite(abs_w)] <- 1
  max_w <- max(abs_w, na.rm = TRUE)
  if (!is.finite(max_w) || max_w <= 0) max_w <- 1
  igraph::E(g)$width <- 1 + 4 * (abs_w / max_w)

  dir.create(dirname(out_svg), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(filename = out_svg, width = 12, height = 10)
  } else {
    grDevices::svg(filename = out_svg, width = 12, height = 10)
  }

  lay <- igraph::layout_with_fr(g)
  plot(
    g,
    layout = lay,
    vertex.label.dist = 0.4,
    edge.curved = 0.08,
    main = ifelse(edge_color_mode == "signed_corr", "Selected Module Network (signed correlation)", "Selected Module Network (edge strength)")
  )
  if (edge_color_mode == "signed_corr") {
    legend(
      "topleft",
      legend = c("Positive correlation", "Negative correlation", "Unknown"),
      col = c(pos_edge_color, neg_edge_color, na_corr_color),
      lwd = 3,
      bty = "n"
    )
  } else {
    legend(
      "topleft",
      legend = c("Edge strength (weight)"),
      col = c(weight_edge_color),
      lwd = 3,
      bty = "n"
    )
  }
  grDevices::dev.off()

  list(
    node_degree = data.frame(
      node = igraph::V(g)$name,
      degree = as.numeric(igraph::V(g)$degree),
      stringsAsFactors = FALSE
    ),
    edge_table = edges
  )
}


# FOCUSED SVG EXPORT FOR SELECTED MODULE

if (isTRUE(export_selected_module_svgs)) {
  minttea_dir_effective <- resolve_existing_minttea_dir(minttea_dir, study_group_col = "disease_state")
  if (!identical(minttea_dir_effective, minttea_dir)) {
    cat("Using detected MintTea output directory:", minttea_dir_effective, "\n")
  }

  selected_svg_dir <- file.path(minttea_dir_effective, "selected_module_svgs")
  dir.create(selected_svg_dir, recursive = TRUE, showWarnings = FALSE)

  minttea_obj_selected <- load_minttea_obj_if_needed(minttea_obj, minttea_dir_effective)
  has_saved_minttea <- minttea_result_rds_exists(minttea_dir_effective)
  if (is.null(minttea_obj_selected)) {
    if (isTRUE(run_minttea) || has_saved_minttea) {
      warning("Could not load MintTea object for selected module SVG export.")
    } else {
      message(
        "Skipping selected module SVG export: no MintTea object in memory and no saved 'minttea_result.rds'. ",
        "Set run_minttea <- TRUE (or provide an existing MintTea output) to export selected modules."
      )
    }
  } else if (!(selected_minttea_setting %in% names(minttea_obj_selected))) {
    warning("Selected setting not found in MintTea object: ", selected_minttea_setting)
  } else {
    available_modules <- names(minttea_obj_selected[[selected_minttea_setting]])
    modules_to_export <- intersect(selected_minttea_modules, available_modules)
    missing_modules <- setdiff(selected_minttea_modules, modules_to_export)

    if (length(missing_modules) > 0) {
      warning("Requested module(s) not found for selected setting: ", paste(missing_modules, collapse = ", "))
    }

    if (length(modules_to_export) == 0) {
      warning("No requested modules were found for selected setting; skipping selected module SVG export.")
    }

    setting_tag <- safe_filename(selected_minttea_setting)

    for (selected_minttea_module in modules_to_export) {
      selected_mod <- minttea_obj_selected[[selected_minttea_setting]][[selected_minttea_module]]
      module_tag <- safe_filename(selected_minttea_module)

      network_weight_svg <- file.path(selected_svg_dir, paste0("network_weight_", setting_tag, "_", module_tag, ".svg"))
      network_signed_corr_svg <- file.path(selected_svg_dir, paste0("network_signedcorr_", setting_tag, "_", module_tag, ".svg"))

      weight_plot_obj <- try(
        plot_module_with_degree_and_signed_edges_svg(
          minttea_module = selected_mod,
          out_svg = network_weight_svg,
          proc_data = out_df,
          corr_method = corr_method,
          edge_color_mode = "weight"
        ),
        silent = TRUE
      )
      signed_plot_obj <- try(
        plot_module_with_degree_and_signed_edges_svg(
          minttea_module = selected_mod,
          out_svg = network_signed_corr_svg,
          proc_data = out_df,
          corr_method = corr_method,
          edge_color_mode = "signed_corr"
        ),
        silent = TRUE
      )

      if (!inherits(weight_plot_obj, "try-error") && is.list(weight_plot_obj) && !is.null(weight_plot_obj$node_degree)) {
        utils::write.csv(
          weight_plot_obj$node_degree,
          file.path(selected_svg_dir, paste0("node_degree_", setting_tag, "_", module_tag, ".csv")),
          row.names = FALSE
        )
      } else {
        warning("Failed to build weight-based selected module network SVG: ", selected_minttea_module)
      }

      if (inherits(signed_plot_obj, "try-error")) {
        warning("Failed to build signed-correlation selected module network SVG: ", selected_minttea_module)
      }

      if (!inherits(signed_plot_obj, "try-error") && is.list(signed_plot_obj) && is.data.frame(signed_plot_obj$edge_table)) {
        edges_out <- signed_plot_obj$edge_table
        utils::write.csv(
          edges_out,
          file.path(selected_svg_dir, paste0("edges_with_weight_and_corr_", setting_tag, "_", module_tag, ".csv")),
          row.names = FALSE
        )
      } else if (is.data.frame(selected_mod$module_edges)) {
        edges_out <- selected_mod$module_edges
        edges_out$edge_weight <- suppressWarnings(as.numeric(edges_out$edge_weight))
        utils::write.csv(
          edges_out,
          file.path(selected_svg_dir, paste0("edges_weight_only_", setting_tag, "_", module_tag, ".csv")),
          row.names = FALSE
        )
      }

      # Cytoscape exports: nodes + edges for both network variants
      node_table <- NULL
      if (!inherits(weight_plot_obj, "try-error") && is.list(weight_plot_obj) && is.data.frame(weight_plot_obj$node_degree)) {
        node_table <- weight_plot_obj$node_degree
      } else if (!inherits(signed_plot_obj, "try-error") && is.list(signed_plot_obj) && is.data.frame(signed_plot_obj$node_degree)) {
        node_table <- signed_plot_obj$node_degree
      }

      if (!is.null(node_table)) {
        node_table$id <- as.character(node_table$node)
        node_table$feature_type <- sub("__.*$", "", node_table$id)
        node_table <- node_table[, c("id", "feature_type", "degree")]
        utils::write.csv(
          node_table,
          file.path(selected_svg_dir, paste0("cytoscape_nodes_", setting_tag, "_", module_tag, ".csv")),
          row.names = FALSE
        )
      }

      cyt_weight_edges <- NULL
      if (is.data.frame(selected_mod$module_edges) && all(c("feature.x", "feature.y") %in% colnames(selected_mod$module_edges))) {
        cyt_weight_edges <- selected_mod$module_edges
        cyt_weight_edges$source <- as.character(cyt_weight_edges$feature.x)
        cyt_weight_edges$target <- as.character(cyt_weight_edges$feature.y)
        cyt_weight_edges$edge_weight <- suppressWarnings(as.numeric(cyt_weight_edges$edge_weight))
        cyt_weight_edges$interaction <- "weight"
        cyt_weight_edges <- cyt_weight_edges[, c("source", "interaction", "target", "edge_weight")]
        utils::write.csv(
          cyt_weight_edges,
          file.path(selected_svg_dir, paste0("cytoscape_edges_weight_", setting_tag, "_", module_tag, ".csv")),
          row.names = FALSE
        )
      }

      cyt_signed_edges <- NULL
      if (!inherits(signed_plot_obj, "try-error") && is.list(signed_plot_obj) && is.data.frame(signed_plot_obj$edge_table)) {
        cyt_signed_edges <- signed_plot_obj$edge_table
        cyt_signed_edges$source <- as.character(cyt_signed_edges$feature.x)
        cyt_signed_edges$target <- as.character(cyt_signed_edges$feature.y)
        cyt_signed_edges$interaction <- ifelse(
          cyt_signed_edges$corr_sign %in% c("positive", "negative"),
          as.character(cyt_signed_edges$corr_sign),
          "unknown"
        )
        cyt_signed_edges <- cyt_signed_edges[, c("source", "interaction", "target", "edge_weight", "signed_corr", "corr_sign")]
        utils::write.csv(
          cyt_signed_edges,
          file.path(selected_svg_dir, paste0("cytoscape_edges_signedcorr_", setting_tag, "_", module_tag, ".csv")),
          row.names = FALSE
        )
      }

      # Cytoscape-ready network objects (GraphML) with directed edges
      if (!is.null(node_table) && is.data.frame(cyt_weight_edges) && nrow(cyt_weight_edges) > 0) {
        export_cytoscape_graphml(
          node_table = node_table,
          edge_table = cyt_weight_edges,
          out_graphml = file.path(selected_svg_dir, paste0("cytoscape_weight_", setting_tag, "_", module_tag, ".graphml"))
        )
      }
      if (!is.null(node_table) && is.data.frame(cyt_signed_edges) && nrow(cyt_signed_edges) > 0) {
        export_cytoscape_graphml(
          node_table = node_table,
          edge_table = cyt_signed_edges,
          out_graphml = file.path(selected_svg_dir, paste0("cytoscape_signedcorr_", setting_tag, "_", module_tag, ".graphml"))
        )
      }

      if (all(c("shuffled_auroc", "auroc", "shuffled_inter_view_corr", "inter_view_corr") %in% names(selected_mod)) &&
          requireNamespace("ggplot2", quietly = TRUE)) {
        stats_plots <- try(plot_module_stats(minttea_obj_selected, selected_minttea_setting, selected_minttea_module), silent = TRUE)
        if (!inherits(stats_plots, "try-error")) {
          ggplot2::ggsave(
            filename = file.path(selected_svg_dir, paste0("stats_auc_", setting_tag, "_", module_tag, ".svg")),
            plot = stats_plots$p1,
            width = 8,
            height = 6,
            dpi = 300,
            device = "svg"
          )
          ggplot2::ggsave(
            filename = file.path(selected_svg_dir, paste0("stats_cor_", setting_tag, "_", module_tag, ".svg")),
            plot = stats_plots$p2,
            width = 8,
            height = 6,
            dpi = 300,
            device = "svg"
          )
        }
      }
    }

    cat("Selected module SVG exports saved to:", selected_svg_dir, "\n")
  }
}


# REPEATED EXTERNAL CV FOR SELECTED SETTING

if (isTRUE(run_selected_setting_repeated_external_cv)) {
  minttea_dir_effective <- resolve_existing_minttea_dir(minttea_dir, study_group_col = "disease_state")
  minttea_obj_selected <- load_minttea_obj_if_needed(minttea_obj, minttea_dir_effective)

  if (is.null(minttea_obj_selected)) {
    warning("Skipping repeated external CV: MintTea object not available.")
  } else if (!(selected_minttea_setting %in% names(minttea_obj_selected))) {
    warning("Skipping repeated external CV: selected setting not found: ", selected_minttea_setting)
  } else {
    available_modules <- names(minttea_obj_selected[[selected_minttea_setting]])
    module_names <- intersect(selected_minttea_modules, available_modules)
    module_names <- module_names[grepl("^module", module_names)]

    if (length(module_names) == 0) {
      warning("Skipping repeated external CV: no selected modules found in selected setting.")
    } else {
      cv_rows <- list()
      cv_summ <- list()
      cv_pred_rows <- list()
      cv_roc_rows <- list()
      kk <- 1L
      ss <- 1L
      pp <- 1L
      qq <- 1L

      for (m in module_names) {
        mod_obj <- minttea_obj_selected[[selected_minttea_setting]][[m]]
        feats <- as.character(mod_obj$features)
        feats <- feats[feats %in% colnames(out_df)]

        if (length(feats) < 2) {
          warning("Skipping module for repeated external CV (too few matched features): ", m)
          next
        }

        cv_obj <- run_repeated_external_cv_for_module(
          proc_data = out_df,
          feature_names = feats,
          repeats = selected_external_cv_repeats,
          folds = selected_external_cv_folds,
          seed = selected_external_cv_seed
        )

        cv_df <- cv_obj$metrics

        if (nrow(cv_df) == 0) {
          warning("No valid repeated external CV folds produced for module: ", m)
          next
        }

        cv_df$settings <- selected_minttea_setting
        cv_df$module <- m
        cv_df$module_size <- length(feats)
        cv_rows[[kk]] <- cv_df
        kk <- kk + 1L

        if (is.data.frame(cv_obj$predictions) && nrow(cv_obj$predictions) > 0) {
          pred_df <- cv_obj$predictions
          pred_df$settings <- selected_minttea_setting
          pred_df$module <- m
          pred_df$module_size <- length(feats)
          cv_pred_rows[[pp]] <- pred_df
          pp <- pp + 1L
        }

        if (is.data.frame(cv_obj$roc_points) && nrow(cv_obj$roc_points) > 0) {
          roc_df <- cv_obj$roc_points
          roc_df$settings <- selected_minttea_setting
          roc_df$module <- m
          roc_df$module_size <- length(feats)
          cv_roc_rows[[qq]] <- roc_df
          qq <- qq + 1L
        }

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

      if (length(cv_rows) > 0) {
        cv_out <- do.call(rbind, cv_rows)
        cv_sum <- do.call(rbind, cv_summ)

        dir.create(dirname(selected_external_cv_results_csv), recursive = TRUE, showWarnings = FALSE)
        utils::write.csv(cv_out, selected_external_cv_results_csv, row.names = FALSE)
        utils::write.csv(cv_sum, selected_external_cv_summary_csv, row.names = FALSE)

        if (length(cv_pred_rows) > 0) {
          cv_pred_out <- do.call(rbind, cv_pred_rows)
          utils::write.csv(cv_pred_out, selected_external_cv_predictions_csv, row.names = FALSE)
          cat("Repeated external CV predictions CSV:", selected_external_cv_predictions_csv, "\n")
        }

        if (length(cv_roc_rows) > 0) {
          cv_roc_out <- do.call(rbind, cv_roc_rows)
          utils::write.csv(cv_roc_out, selected_external_cv_roc_points_csv, row.names = FALSE)
          cat("Repeated external CV ROC points CSV:", selected_external_cv_roc_points_csv, "\n")
        }

        best_idx <- order(cv_sum$auc_mean, decreasing = TRUE, na.last = TRUE)
        cv_sum <- cv_sum[best_idx, , drop = FALSE]

        cat("Repeated external CV results CSV:", selected_external_cv_results_csv, "\n")
        cat("Repeated external CV summary CSV:", selected_external_cv_summary_csv, "\n")
        cat("Top module by repeated external CV AUC:\n")
        print(utils::head(cv_sum, 1))
      } else {
        warning("Repeated external CV requested but no module produced valid CV outputs.")
      }
    }
  }
}

# Plots and runs bootstrapped external validation of MintTea models 

# Requires both genus and species model results
# Applies bootstrapping for accurate validation

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
})

base_dir <- "C:/PhD/Canterbury Dataset/DIABLO"

genus_results <- file.path(
  base_dir,
  "mixomics_outputs_2block_genus_final",
  "minttea_outputs",
  "selected_setting_repeated_external_cv_results.csv"
)

genus_summary <- file.path(
  base_dir,
  "mixomics_outputs_2block_genus_final",
  "minttea_outputs",
  "selected_setting_repeated_external_cv_summary.csv"
)

genus_roc_points <- file.path(
  base_dir,
  "mixomics_outputs_2block_genus_final",
  "minttea_outputs",
  "selected_setting_repeated_external_cv_roc_points.csv"
)

genus_predictions <- file.path(
  base_dir,
  "mixomics_outputs_2block_genus_final",
  "minttea_outputs",
  "selected_setting_repeated_external_cv_predictions.csv"
)

species_results <- file.path(
  base_dir,
  "mixomics_outputs_2block_species_final",
  "minttea_outputs",
  "selected_setting_repeated_external_cv_results.csv"
)

species_summary <- file.path(
  base_dir,
  "mixomics_outputs_2block_species_final",
  "minttea_outputs",
  "selected_setting_repeated_external_cv_summary.csv"
)

species_roc_points <- file.path(
  base_dir,
  "mixomics_outputs_2block_species_final",
  "minttea_outputs",
  "selected_setting_repeated_external_cv_roc_points.csv"
)

species_predictions <- file.path(
  base_dir,
  "mixomics_outputs_2block_species_final",
  "minttea_outputs",
  "selected_setting_repeated_external_cv_predictions.csv"
)

out_svg <- file.path(
  base_dir,
  "mixomics_outputs_2block_consensus",
  "minttea_outputs",
  "consensus_cv_roc_curve_genus_vs_species.svg"
)

out_svg_pooled <- file.path(
  base_dir,
  "mixomics_outputs_2block_consensus",
  "minttea_outputs",
  "consensus_cv_roc_curve_genus_vs_species_pooled_oof.svg"
)

out_svg_pooled_boot_ci <- file.path(
  base_dir,
  "mixomics_outputs_2block_consensus",
  "minttea_outputs",
  "consensus_cv_roc_curve_genus_vs_species_pooled_oof_bootstrap_ci.svg"
)

read_csv_checked <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

read_csv_if_exists <- function(path) {
  if (!file.exists(path)) return(NULL)
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

pick_best_module <- function(summary_df) {
  req <- c("module", "auc_mean")
  if (!all(req %in% colnames(summary_df))) {
    stop("Summary CSV must include columns: ", paste(req, collapse = ", "))
  }
  x <- summary_df[is.finite(summary_df$auc_mean), , drop = FALSE]
  if (nrow(x) == 0) stop("No finite auc_mean rows found in summary CSV.")
  x <- x[order(x$auc_mean, decreasing = TRUE), , drop = FALSE]
  as.character(x$module[1])
}

# Convert AUC into a smooth ROC under equal-variance binormal assumption.
roc_from_auc_binormal <- function(auc_vec, fpr_grid = seq(0, 1, by = 0.0025)) {
  auc_vec <- suppressWarnings(as.numeric(auc_vec))
  auc_vec <- auc_vec[is.finite(auc_vec) & auc_vec > 0 & auc_vec < 1]
  if (length(auc_vec) == 0) {
    stop("No valid AUC values in (0,1) were available for ROC reconstruction.")
  }

  d_vals <- sqrt(2) * stats::qnorm(auc_vec)

  fold_curves <- vapply(
    d_vals,
    function(d) stats::pnorm(stats::qnorm(fpr_grid) + d),
    numeric(length(fpr_grid))
  )

  if (is.null(dim(fold_curves))) {
    fold_curves <- matrix(fold_curves, ncol = 1)
  }

  tpr_mean <- rowMeans(fold_curves, na.rm = TRUE)
  tpr_lo <- apply(fold_curves, 1, stats::quantile, probs = 0.025, na.rm = TRUE)
  tpr_hi <- apply(fold_curves, 1, stats::quantile, probs = 0.975, na.rm = TRUE)

  data.frame(
    fpr = fpr_grid,
    tpr_mean = pmin(pmax(tpr_mean, 0), 1),
    tpr_lo = pmin(pmax(tpr_lo, 0), 1),
    tpr_hi = pmin(pmax(tpr_hi, 0), 1),
    stringsAsFactors = FALSE
  )
}

summarize_fold_roc_curves <- function(roc_df, fpr_grid = seq(0, 1, by = 0.0025)) {
  req <- c("repeat_id", "fold_id", "fpr", "tpr")
  if (!all(req %in% colnames(roc_df))) {
    stop("ROC points CSV must include columns: ", paste(req, collapse = ", "))
  }

  roc_df$fpr <- suppressWarnings(as.numeric(roc_df$fpr))
  roc_df$tpr <- suppressWarnings(as.numeric(roc_df$tpr))
  roc_df <- roc_df[is.finite(roc_df$fpr) & is.finite(roc_df$tpr), , drop = FALSE]
  if (nrow(roc_df) == 0) stop("ROC points table has no finite fpr/tpr rows.")

  fold_key <- paste(roc_df$repeat_id, roc_df$fold_id, sep = "__")
  split_df <- split(roc_df, fold_key)
  if (length(split_df) == 0) stop("No fold-level ROC groups were found.")

  fold_mat <- lapply(split_df, function(df) {
    df <- df[order(df$fpr, df$tpr), c("fpr", "tpr")]
    agg <- aggregate(tpr ~ fpr, data = df, FUN = max)
    agg <- agg[order(agg$fpr), , drop = FALSE]
    agg$tpr <- cummax(pmin(pmax(agg$tpr, 0), 1))

    x <- c(0, agg$fpr, 1)
    y <- c(0, agg$tpr, 1)
    ord <- order(x, y)
    x <- x[ord]
    y <- y[ord]

    xy <- aggregate(y ~ x, FUN = max)
    stats::approx(x = xy$x, y = xy$y, xout = fpr_grid, rule = 2, ties = "ordered")$y
  })

  curve_mat <- do.call(cbind, fold_mat)
  if (is.null(dim(curve_mat))) curve_mat <- matrix(curve_mat, ncol = 1)

  data.frame(
    fpr = fpr_grid,
    tpr_mean = rowMeans(curve_mat, na.rm = TRUE),
    tpr_lo = apply(curve_mat, 1, stats::quantile, probs = 0.025, na.rm = TRUE),
    tpr_hi = apply(curve_mat, 1, stats::quantile, probs = 0.975, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

compute_auc_rank <- function(y_true_binary01, y_score) {
  y <- as.integer(y_true_binary01)
  s <- suppressWarnings(as.numeric(y_score))
  ok <- !is.na(y) & is.finite(s)
  y <- y[ok]
  s <- s[ok]

  n_pos <- sum(y == 1L)
  n_neg <- sum(y == 0L)
  if (length(y) < 2 || n_pos == 0L || n_neg == 0L) return(NA_real_)

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
  if (length(y) < 2 || n_pos == 0L || n_neg == 0L) return(data.frame())

  # Fast ROC construction: sort once, then use cumulative class counts at score change points.
  ord <- order(s, decreasing = TRUE)
  s_ord <- s[ord]
  y_ord <- y[ord]

  tp_cum <- cumsum(y_ord == 1L)
  fp_cum <- cumsum(y_ord == 0L)
  idx_last <- c(which(diff(s_ord) != 0), length(s_ord))

  thr <- s_ord[idx_last]
  tp <- tp_cum[idx_last]
  fp <- fp_cum[idx_last]

  roc <- data.frame(
    threshold = c(Inf, thr, -Inf),
    fpr = c(0, fp / n_neg, 1),
    tpr = c(0, tp / n_pos, 1),
    stringsAsFactors = FALSE
  )
  roc <- roc[order(roc$fpr, roc$tpr), , drop = FALSE]
  rownames(roc) <- NULL
  roc
}

build_pooled_curve_from_predictions <- function(pred_df, best_module, model_label) {
  req <- c("module", "y_true", "y_score")
  if (!all(req %in% colnames(pred_df))) return(NULL)

  x <- pred_df[as.character(pred_df$module) == best_module, , drop = FALSE]
  if (nrow(x) == 0) return(NULL)

  roc <- compute_roc_points(x$y_true, x$y_score)
  if (nrow(roc) == 0) return(NULL)

  roc$model <- model_label
  roc$module <- best_module
  roc$auc_pooled <- compute_auc_rank(x$y_true, x$y_score)
  roc
}

summarize_pooled_bootstrap_roc <- function(pred_df,
                                           best_module,
                                           model_label,
                                           fpr_grid = seq(0, 1, by = 0.01),
                                           n_boot = 500,
                                           seed = 20260309) {
  req <- c("module", "y_true", "y_score")
  if (!all(req %in% colnames(pred_df))) return(NULL)

  x <- pred_df[as.character(pred_df$module) == best_module, c("y_true", "y_score"), drop = FALSE]
  x$y_true <- suppressWarnings(as.integer(x$y_true))
  x$y_score <- suppressWarnings(as.numeric(x$y_score))
  x <- x[!is.na(x$y_true) & is.finite(x$y_score), , drop = FALSE]
  if (nrow(x) < 2) return(NULL)

  pos_idx <- which(x$y_true == 1L)
  neg_idx <- which(x$y_true == 0L)
  if (length(pos_idx) < 2 || length(neg_idx) < 2) return(NULL)

  interpolate_roc_to_grid <- function(roc_in, grid) {
    if (nrow(roc_in) == 0) return(rep(NA_real_, length(grid)))
    tmp <- roc_in[, c("fpr", "tpr"), drop = FALSE]
    tmp <- tmp[is.finite(tmp$fpr) & is.finite(tmp$tpr), , drop = FALSE]
    if (nrow(tmp) == 0) return(rep(NA_real_, length(grid)))
    tmp <- aggregate(tpr ~ fpr, data = tmp, FUN = max)
    tmp <- tmp[order(tmp$fpr), , drop = FALSE]
    tmp$tpr <- cummax(pmin(pmax(tmp$tpr, 0), 1))

    x_all <- c(0, tmp$fpr, 1)
    y_all <- c(0, tmp$tpr, 1)
    ord <- order(x_all, y_all)
    xy <- data.frame(x = x_all[ord], y = y_all[ord], stringsAsFactors = FALSE)
    xy <- aggregate(y ~ x, data = xy, FUN = max)
    stats::approx(x = xy$x, y = xy$y, xout = grid, rule = 2, ties = "ordered")$y
  }

  # Mean pooled ROC from the full pooled OOF predictions.
  roc_mean <- compute_roc_points(x$y_true, x$y_score)
  if (nrow(roc_mean) == 0) return(NULL)
  tpr_mean <- interpolate_roc_to_grid(roc_mean, fpr_grid)

  set.seed(seed)
  boot_curves <- matrix(NA_real_, nrow = length(fpr_grid), ncol = as.integer(n_boot))

  for (b in seq_len(as.integer(n_boot))) {
    bpos <- sample(pos_idx, size = length(pos_idx), replace = TRUE)
    bneg <- sample(neg_idx, size = length(neg_idx), replace = TRUE)
    bb <- x[c(bpos, bneg), , drop = FALSE]

    roc_b <- compute_roc_points(bb$y_true, bb$y_score)
    if (nrow(roc_b) == 0) next

    boot_curves[, b] <- interpolate_roc_to_grid(roc_b, fpr_grid)
  }

  data.frame(
    fpr = fpr_grid,
    tpr_mean = pmin(pmax(tpr_mean, 0), 1),
    tpr_lo = pmin(pmax(apply(boot_curves, 1, stats::quantile, probs = 0.025, na.rm = TRUE), 0), 1),
    tpr_hi = pmin(pmax(apply(boot_curves, 1, stats::quantile, probs = 0.975, na.rm = TRUE), 0), 1),
    model = model_label,
    module = best_module,
    auc_pooled = compute_auc_rank(x$y_true, x$y_score),
    stringsAsFactors = FALSE
  )
}

build_model_curve <- function(results_df, summary_df, model_label, roc_points_df = NULL) {
  req <- c("module", "auc")
  if (!all(req %in% colnames(results_df))) {
    stop("Results CSV must include columns: ", paste(req, collapse = ", "))
  }

  best_module <- pick_best_module(summary_df)
  auc_vals <- results_df$auc[results_df$module == best_module]

  curve_df <- NULL
  if (!is.null(roc_points_df)) {
    if ("module" %in% colnames(roc_points_df)) {
      roc_mod <- roc_points_df[as.character(roc_points_df$module) == best_module, , drop = FALSE]
      if (nrow(roc_mod) > 0) {
        curve_df <- tryCatch(
          summarize_fold_roc_curves(roc_mod),
          error = function(e) NULL
        )
      }
    }
  }

  if (is.null(curve_df)) {
    curve_df <- roc_from_auc_binormal(auc_vals)
  }

  curve_df$model <- model_label
  curve_df$module <- best_module
  curve_df$auc_mean <- mean(as.numeric(auc_vals), na.rm = TRUE)
  curve_df
}

genus_res <- read_csv_checked(genus_results)
genus_sum <- read_csv_checked(genus_summary)
species_res <- read_csv_checked(species_results)
species_sum <- read_csv_checked(species_summary)
genus_roc <- read_csv_if_exists(genus_roc_points)
species_roc <- read_csv_if_exists(species_roc_points)
genus_pred <- read_csv_if_exists(genus_predictions)
species_pred <- read_csv_if_exists(species_predictions)

curve_genus <- build_model_curve(genus_res, genus_sum, "Genus", roc_points_df = genus_roc)
curve_species <- build_model_curve(species_res, species_sum, "Species", roc_points_df = species_roc)
plot_df <- rbind(curve_genus, curve_species)

# Keep explicit model-specific data frames so each CI is rendered independently.
genus_df <- plot_df[plot_df$model == "Genus", , drop = FALSE]
species_df <- plot_df[plot_df$model == "Species", , drop = FALSE]

legend_df <- unique(plot_df[, c("model", "module", "auc_mean")])
legend_df$label <- sprintf("%s (%s, mean AUC = %.3f)", legend_df$model, legend_df$module, legend_df$auc_mean)
label_map <- stats::setNames(legend_df$label, legend_df$model)

palette <- c("Genus" = "#1f77b4", "Species" = "#d62728")

p <- ggplot2::ggplot() +
  ggplot2::geom_ribbon(
    data = genus_df,
    ggplot2::aes(x = fpr, ymin = tpr_lo, ymax = tpr_hi),
    inherit.aes = FALSE,
    fill = palette[["Genus"]],
    alpha = 0.18,
    color = NA
  ) +
  ggplot2::geom_ribbon(
    data = species_df,
    ggplot2::aes(x = fpr, ymin = tpr_lo, ymax = tpr_hi),
    inherit.aes = FALSE,
    fill = palette[["Species"]],
    alpha = 0.18,
    color = NA
  ) +
  ggplot2::geom_line(
    data = genus_df,
    ggplot2::aes(x = fpr, y = tpr_mean, color = model),
    linewidth = 1.15
  ) +
  ggplot2::geom_line(
    data = species_df,
    ggplot2::aes(x = fpr, y = tpr_mean, color = model),
    linewidth = 1.15
  ) +
  ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey55", linewidth = 0.6) +
  ggplot2::scale_color_manual(values = palette, labels = label_map) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "Consensus CV ROC Curve: Genus vs Species",
    subtitle = if (!is.null(genus_roc) && !is.null(species_roc)) {
      "ROC from exported out-of-fold raw ROC points"
    } else {
      "ROC reconstructed from fold-level AUC values (equal-variance binormal model)"
    },
    x = "False Positive Rate",
    y = "True Positive Rate",
    color = NULL
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.box = "vertical"
  )

dir.create(dirname(out_svg), recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggplot2::ggsave(out_svg, plot = p, width = 8.5, height = 7, device = svglite::svglite)
} else {
  ggplot2::ggsave(out_svg, plot = p, width = 8.5, height = 7, device = "svg")
}

# Option 3: pooled out-of-fold ROC (single smoother curve per model, no CI ribbon).
genus_best <- pick_best_module(genus_sum)
species_best <- pick_best_module(species_sum)
pool_genus <- if (!is.null(genus_pred)) build_pooled_curve_from_predictions(genus_pred, genus_best, "Genus") else NULL
pool_species <- if (!is.null(species_pred)) build_pooled_curve_from_predictions(species_pred, species_best, "Species") else NULL

if (!is.null(pool_genus) && !is.null(pool_species)) {
  pool_df <- rbind(pool_genus, pool_species)
  leg_pool <- unique(pool_df[, c("model", "module", "auc_pooled")])
  leg_pool$label <- sprintf("%s (%s, pooled OOF AUC = %.3f)", leg_pool$model, leg_pool$module, leg_pool$auc_pooled)
  label_pool <- stats::setNames(leg_pool$label, leg_pool$model)

  p_pool <- ggplot2::ggplot(pool_df, ggplot2::aes(x = fpr, y = tpr, color = model)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey55", linewidth = 0.6) +
    ggplot2::scale_color_manual(values = palette, labels = label_pool) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "Consensus CV ROC Curve: Genus vs Species (Pooled OOF)",
      subtitle = "Pooled out-of-fold predictions across repeated external CV",
      x = "False Positive Rate",
      y = "True Positive Rate",
      color = NULL
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical"
    )

  if (requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(out_svg_pooled, plot = p_pool, width = 8.5, height = 7, device = svglite::svglite)
  } else {
    ggplot2::ggsave(out_svg_pooled, plot = p_pool, width = 8.5, height = 7, device = "svg")
  }
  cat("Consensus pooled OOF ROC SVG written to:\n", out_svg_pooled, "\n", sep = "")

  boot_genus <- tryCatch(
    summarize_pooled_bootstrap_roc(genus_pred, genus_best, "Genus"),
    error = function(e) {
      message("Bootstrap failed for Genus: ", conditionMessage(e))
      NULL
    }
  )
  boot_species <- tryCatch(
    summarize_pooled_bootstrap_roc(species_pred, species_best, "Species"),
    error = function(e) {
      message("Bootstrap failed for Species: ", conditionMessage(e))
      NULL
    }
  )

  if (!is.null(boot_genus) && !is.null(boot_species)) {
    boot_df <- rbind(boot_genus, boot_species)
    bg <- boot_df[boot_df$model == "Genus", , drop = FALSE]
    bs <- boot_df[boot_df$model == "Species", , drop = FALSE]

    leg_boot <- unique(boot_df[, c("model", "module", "auc_pooled")])
    leg_boot$label <- sprintf("%s (%s, pooled OOF AUC = %.3f)", leg_boot$model, leg_boot$module, leg_boot$auc_pooled)
    label_boot <- stats::setNames(leg_boot$label, leg_boot$model)

    p_pool_boot <- ggplot2::ggplot() +
      ggplot2::geom_ribbon(
        data = bg,
        ggplot2::aes(x = fpr, ymin = tpr_lo, ymax = tpr_hi),
        inherit.aes = FALSE,
        fill = palette[["Genus"]],
        alpha = 0.16,
        color = NA
      ) +
      ggplot2::geom_ribbon(
        data = bs,
        ggplot2::aes(x = fpr, ymin = tpr_lo, ymax = tpr_hi),
        inherit.aes = FALSE,
        fill = palette[["Species"]],
        alpha = 0.16,
        color = NA
      ) +
      ggplot2::geom_line(
        data = bg,
        ggplot2::aes(x = fpr, y = tpr_mean, color = model),
        linewidth = 1.15
      ) +
      ggplot2::geom_line(
        data = bs,
        ggplot2::aes(x = fpr, y = tpr_mean, color = model),
        linewidth = 1.15
      ) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey55", linewidth = 0.6) +
      ggplot2::scale_color_manual(values = palette, labels = label_boot) +
      ggplot2::coord_equal() +
      ggplot2::labs(
        title = "Consensus CV ROC: Pooled OOF with Bootstrap CI",
        subtitle = "Stratified bootstrap (95% CI) on pooled out-of-fold predictions",
        x = "False Positive Rate",
        y = "True Positive Rate",
        color = NULL
      ) +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(
        legend.position = "bottom",
        legend.box = "vertical"
      )

    if (requireNamespace("svglite", quietly = TRUE)) {
      ggplot2::ggsave(out_svg_pooled_boot_ci, plot = p_pool_boot, width = 8.5, height = 7, device = svglite::svglite)
    } else {
      ggplot2::ggsave(out_svg_pooled_boot_ci, plot = p_pool_boot, width = 8.5, height = 7, device = "svg")
    }
    cat("Consensus pooled OOF ROC (bootstrap CI) SVG written to:\n", out_svg_pooled_boot_ci, "\n", sep = "")
  } else {
    cat("Pooled OOF ROC bootstrap CI SVG skipped (bootstrap summary failed for genus and/or species).\n")
  }
} else {
  cat("Pooled OOF ROC SVG skipped (missing predictions for genus and/or species).\n")
}

cat("Consensus ROC SVG written to:\n", out_svg, "\n", sep = "")

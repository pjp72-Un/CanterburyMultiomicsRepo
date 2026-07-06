# =====================================================================
# CROSS-BLOCK DIABLO (block sPLS-DA)  —  multivariate alternative to the
# per-edge correlation network (Cross block network.R)
# ---------------------------------------------------------------------
# Rationale:
#   At n = 27 the per-edge Spearman network is underpowered: glmnet-driver
#   grids pass only a handful of edges, and the full feature grid leaves no
#   single edge bootstrap-stable. DIABLO (mixOmics block sPLS-DA) instead
#   borrows strength ACROSS features — it finds sparse latent components that
#   are simultaneously (a) correlated between the microbiome and metabolome
#   blocks and (b) discriminative for CeD vs Control. The deliverable is a
#   compact, CV-stable set of co-selected features plus a relevance network
#   built from the latent structure, rather than thousands of fragile edges.
#
# Inputs (aligned_blocks/):
#   aligned_blocks_<level>_<prev>_status.rds  ->  list(Microbiome, Metabolome, Group)
#     Microbiome : n x p  (CLR-transformed, prevalence-filtered)
#     Metabolome : n x m  (log2-transformed)
#     Group      : factor, binary CeD / Control  (the *_status files)
#
# Outputs (cross_block_diablo_<level>/):
#   diablo_model_<level>.rds                 fitted final model + tuning
#   diablo_selected_features_<level>.csv     selected features per comp + stability
#   diablo_network_edges_<level>.csv         cross-block relevance-network edges
#   diablo_performance_<level>.csv           CV error rate / AUC per component
#   diablo_cv_keepX_<level>.csv              tuned keepX grid result
#   plots: sample (plotIndiv), plotDiablo, loadings, circos, network, cross-block scatter  (svg+png)
# =====================================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(mixOmics)
  library(ggplot2)
  library(dplyr)
  library(igraph)
})

# ---------------------------------------------------------------------
# USER CONFIG
# ---------------------------------------------------------------------
base_dir   <- "C:/PhD/Canterbury Dataset/crossomic"
block_dir  <- file.path(base_dir, "aligned_blocks")

level      <- "species"     # "genus" | "species"  (genus won the feature scout)
prev       <- "30"        # prevalence-filter tag in the filename
# Outcome / design:
#   "status" = binary CeD vs Control            (pooled disease; main analysis)
#   "group"  = three-way ACeD / Control / TCeD   (untreated vs control vs treated;
#              adds a treatment axis the binary model collapses, but n=27 split
#              three ways is underpowered — exploratory)
outcome    <- "status"    # "status" | "group"
aligned_rds <- file.path(block_dir, sprintf("aligned_blocks_%s_%s_%s.rds", level, prev, outcome))

out_dir    <- file.path(base_dir, sprintf("cross_block_diablo_%s_%s", level, outcome))

# DIABLO design: off-diagonal weight links the two X blocks.
#   ~0.1  -> prioritise CLASS DISCRIMINATION  (DIABLO behaves like sgccda for Y)
#   ~1.0  -> prioritise CROSS-BLOCK CORRELATION (recover shared latent structure)
#   0.5   -> balanced compromise (used here: we want a correlated, interpretable
#            module that still separates CeD from Control)
design_weight <- 0.5

# Components retained. The binary disease contrast is a single axis, so comp1 is
# the deliverable (comp2 added no CV discrimination and was unstable). The
# three-way outcome needs two axes (disease + treatment), so it keeps comp2.
ncomp_final <- if (identical(outcome, "group")) 2L else 1L
ncomp_max   <- ncomp_final   # components evaluated during tuning / perf

# keepX grid searched per block during tuning (features kept per component).
test_keepX_micro <- c(5, 8, 10, 15, 20, 30)
test_keepX_metab <- c(5, 8, 10, 15, 20, 30)

do_tune    <- TRUE        # FALSE -> skip CV tuning, use default_keepX below
default_keepX <- list(Microbiome = c(5, 5), Metabolome = c(20, 8))

cv_dist    <- "centroids.dist"   # distance for class assignment in CV
validation <- "loo"       # "loo" (deterministic, ideal at n=27) | "Mfold"
folds      <- 5           # used only when validation = "Mfold"
nrepeat    <- 1           # >1 only meaningful for Mfold
tune_seed  <- 27L

net_cutoff <- 0.50        # |similarity| cutoff for relevance-network edges/circos

# Colorblind-safe palette (matches the rest of the chapter)
col_ced     <- "#D55E00"
col_control <- "#0072B2"
col_treated <- "#009E73"
group_cols  <- if (identical(outcome, "group"))
  c(Control = col_control, ACeD = col_ced, TCeD = col_treated) else
  c(Control = col_control, CeD = col_ced)

# ---------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ab <- readRDS(aligned_rds)
stopifnot(all(c("Microbiome", "Metabolome", "Group") %in% names(ab)))

X <- list(
  Microbiome = as.matrix(ab$Microbiome),
  Metabolome = as.matrix(ab$Metabolome)
)
stopifnot(identical(rownames(X$Microbiome), rownames(X$Metabolome)))

# Outcome: binary CeD vs Control (status blocks) or three-way ACeD/Control/TCeD
# (group blocks). Direction labels downstream assume the last level is the
# "disease/treated" reference for sign interpretation.
if (identical(outcome, "group")) {
  Y <- factor(as.character(ab$Group), levels = c("Control", "ACeD", "TCeD"))
} else {
  Y <- factor(ifelse(as.character(ab$Group) == "Control", "Control", "CeD"),
              levels = c("Control", "CeD"))
}
Y <- droplevels(Y)
stopifnot(length(Y) == nrow(X$Microbiome))

cat(sprintf("DIABLO | level=%s prev=%s | n=%d | taxa=%d | metabolites=%d\n",
            level, prev, nrow(X$Microbiome), ncol(X$Microbiome), ncol(X$Metabolome)))
cat("Class balance:", paste(names(table(Y)), table(Y), sep = "="), "\n")

# ---------------------------------------------------------------------
# DESIGN MATRIX  (blocks x blocks; Y link is automatic)
# ---------------------------------------------------------------------
design <- matrix(design_weight, nrow = length(X), ncol = length(X),
                 dimnames = list(names(X), names(X)))
diag(design) <- 0

# ---------------------------------------------------------------------
# TUNE keepX  (CV over the feature grid)
# ---------------------------------------------------------------------
if (do_tune) {
  cat("\nTuning keepX via", validation, "CV (this can take a few minutes)...\n")
  set.seed(tune_seed)
  test_keepX <- list(Microbiome = test_keepX_micro, Metabolome = test_keepX_metab)
  tune_res <- tune.block.splsda(
    X = X, Y = Y, ncomp = ncomp_max,
    test.keepX = test_keepX, design = design,
    validation = validation, folds = folds, nrepeat = nrepeat,
    dist = cv_dist, BPPARAM = BiocParallel::SerialParam()
  )
  keepX <- tune_res$choice.keepX
  # ensure each block has ncomp_final entries (pad with last value if needed)
  keepX <- lapply(keepX, function(v) {
    if (length(v) >= ncomp_final) v[seq_len(ncomp_final)]
    else c(v, rep(tail(v, 1), ncomp_final - length(v)))
  })
  cv_tab <- data.frame(
    block = rep(names(keepX), each = ncomp_final),
    comp  = rep(seq_len(ncomp_final), times = length(keepX)),
    keepX = unlist(keepX)
  )
  write.csv(cv_tab, file.path(out_dir, sprintf("diablo_cv_keepX_%s.csv", level)),
            row.names = FALSE)
  cat("Tuned keepX:\n"); print(cv_tab, row.names = FALSE)
} else {
  keepX <- default_keepX
  cat("\nUsing default keepX (tuning skipped):\n"); print(keepX)
}

# ---------------------------------------------------------------------
# FIT FINAL MODEL
# ---------------------------------------------------------------------
model <- block.splsda(X = X, Y = Y, ncomp = ncomp_final,
                      keepX = keepX, design = design)

# Realised cross-block correlation of the latent variates (DIABLO analogue of
# MintTea's cross-view correlation). comp-1 is the discriminative/shared axis.
xcor <- sapply(seq_len(ncomp_final), function(k)
  stats::cor(model$variates$Microbiome[, k], model$variates$Metabolome[, k]))
cat("\nCross-block latent correlation (Microbiome vs Metabolome) per comp:\n")
print(round(setNames(xcor, paste0("comp", seq_len(ncomp_final))), 3))

# ---------------------------------------------------------------------
# PERFORMANCE  (CV error rate per component + block AUC)
# ---------------------------------------------------------------------
set.seed(tune_seed)
perf_res <- perf(model, validation = validation, folds = folds,
                 nrepeat = nrepeat, dist = cv_dist,
                 BPPARAM = BiocParallel::SerialParam())

perf_rows <- list()
# Weighted-vote (consensus across blocks) error rate per component.
vote <- tryCatch(perf_res$WeightedVote.error.rate[[cv_dist]],
                 error = function(e) NULL)
if (!is.null(vote)) {
  for (rn in c("Overall.ER", "Overall.BER")) {
    if (rn %in% rownames(vote)) {
      perf_rows[[rn]] <- data.frame(
        metric = paste0("WeightedVote_", rn),
        comp   = seq_len(ncol(vote)),
        value  = as.numeric(vote[rn, ])
      )
    }
  }
}
# Block-wise AUC per component via auroc (resubstitution AUC; optimistic, so the
# CV WeightedVote BER above is the honest discrimination metric).
auc_res <- tryCatch(
  auroc(model, print = FALSE, plot = FALSE),
  error = function(e) NULL)
if (!is.null(auc_res)) {
  for (blk in names(auc_res)) {
    comps <- auc_res[[blk]]
    if (is.null(comps)) next
    vals <- vapply(comps, function(cc) as.numeric(cc[1, 1]), numeric(1))
    perf_rows[[paste0("AUC_", blk)]] <- data.frame(
      metric = paste0("AUCresub_", blk),
      comp   = seq_along(vals),
      value  = as.numeric(vals)
    )
  }
}
perf_tab <- dplyr::bind_rows(perf_rows)
if (nrow(perf_tab)) {
  write.csv(perf_tab, file.path(out_dir, sprintf("diablo_performance_%s.csv", level)),
            row.names = FALSE)
  cat("\nPerformance (CV):\n"); print(perf_tab, row.names = FALSE)
}

# ---------------------------------------------------------------------
# SELECTED FEATURES + CV STABILITY
#   perf()$features$stable gives the frequency with which each feature is
#   selected across CV folds — the multivariate analogue of edge stability.
# ---------------------------------------------------------------------
# perf()$features$stable is nested as stable[[nrep]][[block]][[comp]], where each
# leaf is a named frequency table (selection rate of that feature across CV folds).
stab <- tryCatch(perf_res$features$stable[[1]], error = function(e) NULL)

get_stability <- function(block, comp, feat) {
  if (is.null(stab)) return(NA_real_)
  s <- stab[[block]][[paste0("comp", comp)]]
  if (is.null(s)) return(NA_real_)
  val <- s[feat]
  if (length(val) == 0 || is.na(val)) 0 else as.numeric(val)
}

sel_rows <- list()
for (k in seq_len(ncomp_final)) {
  sv <- selectVar(model, comp = k)
  for (blk in names(X)) {
    if (is.null(sv[[blk]])) next
    feats <- sv[[blk]]$name
    loads <- sv[[blk]]$value$value.var
    if (length(feats) == 0) next
    sel_rows[[paste(blk, k)]] <- data.frame(
      block     = blk,
      comp      = k,
      feature   = feats,
      loading   = as.numeric(loads),
      direction = ifelse(as.numeric(loads) >= 0, "CeD_up", "Control_up"),
      cv_stability = vapply(feats, function(f) get_stability(blk, k, f), numeric(1)),
      stringsAsFactors = FALSE
    )
  }
}
sel_tab <- dplyr::bind_rows(sel_rows) %>%
  dplyr::arrange(.data$comp, .data$block, dplyr::desc(abs(.data$loading)))
write.csv(sel_tab, file.path(out_dir, sprintf("diablo_selected_features_%s.csv", level)),
          row.names = FALSE)
cat(sprintf("\nSelected features: %d total across %d components\n",
            nrow(sel_tab), ncomp_final))
print(sel_tab %>% dplyr::group_by(.data$comp, .data$block) %>%
        dplyr::summarise(n = dplyr::n(), .groups = "drop"), n = Inf)

# ---------------------------------------------------------------------
# RELEVANCE NETWORK  (cross-block similarity from the latent structure)
#   network() returns, per block pair, a similarity matrix M_<b1>_<b2> based on
#   the selected features' latent associations. We export it as an edge table so
#   the microbe<->metabolite links are auditable, then draw the graph.
# ---------------------------------------------------------------------
edge_tab <- NULL
net <- tryCatch({
  grDevices::pdf(NULL)              # network() opens a device; capture & discard
  on.exit(grDevices::dev.off(), add = TRUE)
  network(model, blocks = c(1, 2), cutoff = net_cutoff,
          color.node = c(col_control, col_ced))
}, error = function(e) { message("network() failed: ", conditionMessage(e)); NULL })

if (!is.null(net)) {
  simM <- net[[grep("^M", names(net), value = TRUE)[1]]]   # micro x metab similarity
  if (!is.null(simM)) {
    edge_tab <- as.data.frame(as.table(as.matrix(simM)),
                              stringsAsFactors = FALSE)
    names(edge_tab) <- c("microbe", "metabolite", "similarity")
    edge_tab <- edge_tab %>%
      dplyr::mutate(abs_sim = abs(.data$similarity),
                    passes_cutoff = .data$abs_sim >= net_cutoff,
                    sign = ifelse(.data$similarity >= 0, "positive", "negative")) %>%
      dplyr::arrange(dplyr::desc(.data$abs_sim))
    write.csv(edge_tab, file.path(out_dir, sprintf("diablo_network_edges_%s.csv", level)),
              row.names = FALSE)
    cat(sprintf("\nRelevance-network edges: %d total, %d with |similarity| >= %.2f\n",
                nrow(edge_tab), sum(edge_tab$passes_cutoff), net_cutoff))
  }
}

# ---------------------------------------------------------------------
# PLOTS  (svg + png)
# ---------------------------------------------------------------------
save_plot <- function(name, expr, w = 18 / 2.54, h = 18 / 2.54) {
  for (dev in c("svg", "png")) {
    f <- file.path(out_dir, sprintf("%s_%s.%s", name, level, dev))
    if (dev == "svg") grDevices::svg(f, width = w, height = h)
    else grDevices::png(f, width = w, height = h, units = "in", res = 300)
    ok <- tryCatch({ force(expr); TRUE },
                   error = function(e) { message(name, " failed: ", conditionMessage(e)); FALSE })
    grDevices::dev.off()
    if (!ok && file.exists(f)) unlink(f)
  }
}

# Sample scatter needs >= 2 components; only drawn for the three-way model.
if (ncomp_final >= 2) {
  save_plot("diablo_sample", plotIndiv(
    model, ind.names = FALSE, legend = TRUE, ellipse = TRUE,
    col.per.group = group_cols, title = sprintf("DIABLO (%s) — samples", level)))
}

# Comp-1 latent scores by group (the honest 1-D view; works for any ncomp).
score_df <- rbind(
  data.frame(group = Y, block = "Microbiome",
             score = model$variates$Microbiome[, 1]),
  data.frame(group = Y, block = "Metabolome",
             score = model$variates$Metabolome[, 1])
)
p_score <- ggplot2::ggplot(score_df,
    ggplot2::aes(x = .data$group, y = .data$score, fill = .data$group)) +
  ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  ggplot2::geom_jitter(width = 0.15, height = 0, size = 1.6, alpha = 0.8) +
  ggplot2::facet_wrap(~ block) +
  ggplot2::scale_fill_manual(values = group_cols) +
  ggplot2::labs(x = NULL, y = "Component 1 latent score",
                title = sprintf("DIABLO (%s) — comp 1 scores", level)) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(legend.position = "none")
ggplot2::ggsave(file.path(out_dir, sprintf("diablo_comp1_scores_%s.svg", level)),
                p_score, width = 18 / 2.54, height = 10 / 2.54)
ggplot2::ggsave(file.path(out_dir, sprintf("diablo_comp1_scores_%s.png", level)),
                p_score, width = 18 / 2.54, height = 10 / 2.54, dpi = 300)

# Cross-block scatter: microbiome vs metabolome comp-1 latent score, one point
# per sample. This is the direct visual of the cross-block correlation (xcor[1])
# reported above — it shows that the two blocks co-vary per sample without any
# edge-thresholding, the honest one-component analogue of the relevance network.
xcb_df <- data.frame(
  group       = Y,
  micro_score = model$variates$Microbiome[, 1],
  metab_score = model$variates$Metabolome[, 1]
)
# Within-group correlations: how much cross-block coupling survives once the
# CeD/Control separation along the shared axis is removed (pooled r is inflated
# by the between-group gap because DIABLO builds the axis to discriminate).
xcb_within <- sapply(split(xcb_df, xcb_df$group), function(d)
  if (nrow(d) > 2) stats::cor(d$micro_score, d$metab_score) else NA_real_)
xcb_sub <- paste0(
  sprintf("comp 1 cross-block r = %.2f", xcor[1]),
  paste0(c("", sprintf("  |  %s r = %.2f", names(xcb_within), xcb_within)), collapse = "")
)
p_xcb <- ggplot2::ggplot(xcb_df,
    ggplot2::aes(x = .data$micro_score, y = .data$metab_score, colour = .data$group)) +
  ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                       colour = "grey40", linewidth = 0.6,
                       inherit.aes = FALSE,
                       ggplot2::aes(x = .data$micro_score, y = .data$metab_score)) +
  ggplot2::geom_point(size = 2.2, alpha = 0.85) +
  ggplot2::scale_colour_manual(values = group_cols) +
  ggplot2::labs(
    x = "Microbiome component 1 score",
    y = "Metabolome component 1 score",
    colour = NULL,
    title = sprintf("DIABLO (%s) \u2014 cross-block scores", level),
    subtitle = xcb_sub
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(out_dir, sprintf("diablo_crossblock_scatter_%s.svg", level)),
                p_xcb, width = 12 / 2.54, height = 12 / 2.54)
ggplot2::ggsave(file.path(out_dir, sprintf("diablo_crossblock_scatter_%s.png", level)),
                p_xcb, width = 12 / 2.54, height = 12 / 2.54, dpi = 300)

save_plot("diablo_blocks", plotDiablo(model, ncomp = 1))

# Component-1 loadings as a robust ggplot barplot (plotLoadings chokes on long
# metabolite names). Built from the selected-feature table.
load_df <- sel_tab %>%
  dplyr::filter(.data$comp == 1) %>%
  dplyr::mutate(
    flab = ifelse(nchar(.data$feature) > 34,
                  paste0(substr(.data$feature, 1, 32), "\u2026"), .data$feature)
  )
load_df$flab <- factor(load_df$flab, levels = load_df$flab[order(load_df$loading)])
p_load <- ggplot2::ggplot(load_df,
    ggplot2::aes(x = .data$loading, y = .data$flab, fill = .data$direction)) +
  ggplot2::geom_col() +
  ggplot2::facet_grid(block ~ ., scales = "free_y", space = "free_y") +
  ggplot2::scale_fill_manual(values = c(CeD_up = col_ced, Control_up = col_control)) +
  ggplot2::labs(x = "Component 1 loading", y = NULL, fill = NULL,
                title = sprintf("DIABLO (%s) \u2014 comp 1 loadings", level)) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(legend.position = "top",
                 panel.grid.major.y = ggplot2::element_blank())
ggplot2::ggsave(file.path(out_dir, sprintf("diablo_loadings_%s.svg", level)),
                p_load, width = 18 / 2.54, height = 18 / 2.54)
ggplot2::ggsave(file.path(out_dir, sprintf("diablo_loadings_%s.png", level)),
                p_load, width = 18 / 2.54, height = 18 / 2.54, dpi = 300)

if (!is.null(net)) {
  save_plot("diablo_network", network(
    model, blocks = c(1, 2), cutoff = net_cutoff,
    color.node = c(col_control, col_ced)))
  save_plot("diablo_circos", circosPlot(
    model, cutoff = net_cutoff, line = TRUE,
    color.blocks = c(col_control, col_ced), size.variables = 0.6))
}

# ---------------------------------------------------------------------
# SAVE MODEL
# ---------------------------------------------------------------------
saveRDS(list(model = model,
             keepX = keepX,
             design = design,
             cross_block_cor = xcor,
             performance = if (exists("perf_tab")) perf_tab else NULL,
             selected = sel_tab,
             edges = edge_tab),
        file.path(out_dir, sprintf("diablo_model_%s.rds", level)))

cat("\nDone. Outputs written to:", out_dir, "\n")

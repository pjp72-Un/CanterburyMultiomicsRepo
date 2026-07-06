suppressPackageStartupMessages({
  library(phyloseq)
  library(ANCOMBC)
  library(dplyr)
})

options(stringsAsFactors = FALSE)


# USER CONFIG

seq_object_dir <- "C:/PhD/Sequencing Results/Files"
out_dir <- "C:/PhD/Canterbury Dataset/16S/robustness/species"

n_boot <- 500
n_perm <- 0
seed <- 27

tax_rank <- "Species"
prv_cut <- 0.15
lib_cut <- 1000
alpha_cutoff <- 0.05
q_cutoff <- 0.05
lfc_thresh <- 0.50
n_cl <- 4

# Criteria for flagging a taxon as bootstrap-robust (all three must hold):
#   1. selection_frequency >= robust_selfreq_min  (re-selected across resamples)
#   2. 95% bootstrap LFC interval excludes zero    (stable effect direction)
#   3. rank_median within the top robust_rank_top_frac of testable features
#      (relative cut, since the number of ranked features differs per contrast).
robust_selfreq_min  <- 0.55
robust_rank_top_frac <- 0.30

analysis_configs <- list(
  list(
    analysis = "untreated_vs_control",
    rds_file = file.path(seq_object_dir, "Cant_5%_untreated.rds"),
    group_var = "CD",
    source_group_var = "Diagnosed.CD",
    subject_var = "Subject",
    paired = FALSE,
    fix_formula = "CD + Age + Sex",
    rand_formula = NULL,
    reference_level = "FALSE",
    alternate_level = "TRUE"
  ),
  list(
    analysis = "treated_vs_control",
    rds_file = file.path(seq_object_dir, "Cant_5%_treated.rds"),
    group_var = "CD",
    source_group_var = "Diagnosed.CD",
    subject_var = "Subject",
    paired = FALSE,
    fix_formula = "CD + Age + Sex",
    rand_formula = NULL,
    reference_level = "FALSE",
    alternate_level = "TRUE"
  ),
  list(
    analysis = "paired_treated_vs_untreated",
    rds_file = file.path(seq_object_dir, "Cant_5%_untreated_treated.rds"),
    group_var = "GFD",
    source_group_var = "GFD",
    subject_var = "Subject",
    paired = TRUE,
    fix_formula = "GFD",
    rand_formula = "(1 | Subject)",
    reference_level = "FALSE",
    alternate_level = "TRUE"
  )
)

get_cli_option <- function(args, name) {
  option_prefix <- paste0("--", name, "=")
  matches <- args[startsWith(args, option_prefix)]

  if (length(matches) > 1L) {
    stop("Command-line option supplied more than once: --", name)
  }
  if (length(matches) == 0L) {
    return(NULL)
  }

  sub(option_prefix, "", matches[[1]], fixed = TRUE)
}

parse_cli_nonnegative_integer <- function(args, name, default_value) {
  raw_value <- get_cli_option(args, name)
  if (is.null(raw_value)) {
    return(default_value)
  }

  parsed_value <- suppressWarnings(as.integer(raw_value))
  if (length(parsed_value) != 1L || is.na(parsed_value) || parsed_value < 0L) {
    stop("Command-line option must be a single non-negative integer: --", name)
  }

  parsed_value
}

parse_cli_string <- function(args, name, default_value) {
  raw_value <- get_cli_option(args, name)
  if (is.null(raw_value)) {
    return(default_value)
  }
  if (!nzchar(raw_value)) {
    stop("Command-line option must not be empty: --", name)
  }

  raw_value
}

cli_args <- commandArgs(trailingOnly = TRUE)
n_boot <- parse_cli_nonnegative_integer(cli_args, "n-boot", n_boot)
n_perm <- parse_cli_nonnegative_integer(cli_args, "n-perm", n_perm)
out_dir <- parse_cli_string(cli_args, "out-dir", out_dir)


# HELPERS

safe_factor <- function(x, levels = NULL) {
  x <- as.character(x)
  if (is.null(levels)) {
    factor(x)
  } else {
    factor(x, levels = levels)
  }
}

validate_prepared_sample_data <- function(sd, cfg, require_subject = FALSE) {
  required <- cfg$group_var
  if (require_subject) {
    required <- c(required, cfg$subject_var)
  }

  missing <- required[!(required %in% colnames(sd))]
  if (length(missing) > 0L) {
    stop(
      "Prepared sample_data is missing required columns: ",
      paste(missing, collapse = ", "),
      ". Run prepare_phyloseq() before resampling or permutation."
    )
  }
}

build_phyloseq_with_duplicates <- function(ps, sample_ids, new_sample_ids, subject_var = NULL, new_subject_ids = NULL) {
  otu <- as(otu_table(ps), "matrix")
  taxa_rows <- taxa_are_rows(ps)

  if (taxa_rows) {
    otu_new <- otu[, sample_ids, drop = FALSE]
    colnames(otu_new) <- new_sample_ids
  } else {
    otu_new <- otu[sample_ids, , drop = FALSE]
    rownames(otu_new) <- new_sample_ids
  }

  sd <- as.data.frame(sample_data(ps), check.names = FALSE)[sample_ids, , drop = FALSE]
  rownames(sd) <- new_sample_ids
  if (!is.null(subject_var) && !is.null(new_subject_ids)) {
    sd[[subject_var]] <- new_subject_ids
  }

  parts <- list(
    otu_table(otu_new, taxa_are_rows = taxa_rows),
    sample_data(sd)
  )

  tax <- tryCatch(tax_table(ps), error = function(e) NULL)
  if (!is.null(tax)) parts[[length(parts) + 1L]] <- tax

  tree <- tryCatch(phy_tree(ps), error = function(e) NULL)
  if (!is.null(tree)) parts[[length(parts) + 1L]] <- tree

  ps_new <- do.call(phyloseq, parts)

  ref <- tryCatch(refseq(ps), error = function(e) NULL)
  if (!is.null(ref)) {
    ps_new <- merge_phyloseq(ps_new, ref)
  }

  ps_new
}

bootstrap_unpaired_phyloseq <- function(ps, cfg) {
  sd <- as.data.frame(sample_data(ps), check.names = FALSE)
  validate_prepared_sample_data(sd, cfg, require_subject = FALSE)
  sd$sample_id <- rownames(sd)

  groups <- split(sd$sample_id, sd[[cfg$group_var]])
  draw_ids <- unlist(lapply(groups, function(ids) sample(ids, length(ids), replace = TRUE)), use.names = FALSE)
  new_ids <- paste0("boot_sample_", seq_along(draw_ids))

  build_phyloseq_with_duplicates(ps, draw_ids, new_ids)
}

bootstrap_paired_phyloseq <- function(ps, cfg) {
  sd <- as.data.frame(sample_data(ps), check.names = FALSE)
  validate_prepared_sample_data(sd, cfg, require_subject = TRUE)
  sd$sample_id <- rownames(sd)

  subject_ids <- unique(as.character(sd[[cfg$subject_var]]))
  draw_subjects <- sample(subject_ids, length(subject_ids), replace = TRUE)

  sample_ids <- character(0)
  new_ids <- character(0)
  new_subject_ids <- character(0)
  cursor <- 1L

  for (i in seq_along(draw_subjects)) {
    sid <- draw_subjects[i]
    ids <- sd$sample_id[as.character(sd[[cfg$subject_var]]) == sid]
    n_ids <- length(ids)
    sample_ids <- c(sample_ids, ids)
    new_ids <- c(new_ids, paste0("boot_sample_", cursor:(cursor + n_ids - 1L)))
    new_subject_ids <- c(new_subject_ids, rep(paste0("boot_subject_", i), n_ids))
    cursor <- cursor + n_ids
  }

  build_phyloseq_with_duplicates(ps, sample_ids, new_ids, cfg$subject_var, new_subject_ids)
}

permute_unpaired_phyloseq <- function(ps, cfg) {
  sd <- as.data.frame(sample_data(ps), check.names = FALSE)
  validate_prepared_sample_data(sd, cfg, require_subject = FALSE)
  levs <- levels(sd[[cfg$group_var]])
  sd[[cfg$group_var]] <- factor(sample(as.character(sd[[cfg$group_var]])), levels = levs)
  sample_data(ps) <- sample_data(sd)
  ps
}

permute_paired_phyloseq <- function(ps, cfg) {
  sd <- as.data.frame(sample_data(ps), check.names = FALSE)
  validate_prepared_sample_data(sd, cfg, require_subject = TRUE)
  levs <- levels(sd[[cfg$group_var]])
  group_vals <- as.character(sd[[cfg$group_var]])

  for (sid in unique(as.character(sd[[cfg$subject_var]]))) {
    idx <- which(as.character(sd[[cfg$subject_var]]) == sid)
    if (length(unique(group_vals[idx])) < 2L) next
    if (runif(1) < 0.5) {
      group_vals[idx] <- rev(group_vals[idx])
    }
  }

  sd[[cfg$group_var]] <- factor(group_vals, levels = levs)
  sample_data(ps) <- sample_data(sd)
  ps
}

make_taxon_labels <- function(ps, taxrank = "Genus") {
  tax_df <- as.data.frame(tax_table(ps), check.names = FALSE)
  if (!(taxrank %in% colnames(tax_df))) {
    return(taxa_names(ps))
  }

  labels <- trimws(as.character(tax_df[[taxrank]]))

  # For species-level labels, prepend the genus so the (often bare) epithet is
  # unambiguous (e.g. "Bacteroides ovatus" rather than just "ovatus").
  if (identical(taxrank, "Species") && "Genus" %in% colnames(tax_df)) {
    genus <- trimws(as.character(tax_df[["Genus"]]))
    have_label <- !(is.na(labels) | labels == "")
    have_genus <- !(is.na(genus) | genus == "")
    combine <- have_label & have_genus
    labels[combine] <- paste(genus[combine], labels[combine])
  }

  missing <- is.na(labels) | labels == ""
  labels[missing] <- taxa_names(ps)[missing]
  make.unique(labels)
}

prepare_phyloseq <- function(ps, cfg) {
  if (identical(tax_rank, "Species")) {
    # Build a genus-aware species key so tax_glom does not merge identical
    # epithets across different genera, then drop ASVs lacking a species call.
    tax_df <- as.data.frame(tax_table(ps), check.names = FALSE)
    genus <- trimws(as.character(tax_df[["Genus"]]))
    species <- trimws(as.character(tax_df[["Species"]]))
    have_sp <- !(is.na(species) | species == "")
    have_ge <- !(is.na(genus) | genus == "")
    tax_df[["GenusSpecies"]] <- ifelse(
      have_sp & have_ge, paste(genus, species),
      ifelse(have_sp, species, NA_character_)
    )
    tax_table(ps) <- phyloseq::tax_table(as.matrix(tax_df))
    ps <- tax_glom(ps, taxrank = "GenusSpecies", NArm = TRUE)
    taxa_names(ps) <- make_taxon_labels(ps, taxrank = "Species")
  } else {
    ps <- tax_glom(ps, taxrank = tax_rank)
    taxa_names(ps) <- make_taxon_labels(ps, taxrank = tax_rank)
  }

  sd <- as.data.frame(sample_data(ps), check.names = FALSE)
  required_cols <- cfg$source_group_var
  if (!is.null(cfg$subject_var)) {
    required_cols <- c(required_cols, cfg$subject_var)
  }
  missing_cols <- required_cols[!(required_cols %in% colnames(sd))]
  if (length(missing_cols) > 0L) {
    stop("sample_data is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  sd[[cfg$group_var]] <- safe_factor(sd[[cfg$source_group_var]], levels = c(cfg$reference_level, cfg$alternate_level))
  if (!is.null(cfg$subject_var) && cfg$subject_var %in% colnames(sd)) {
    sd[[cfg$subject_var]] <- safe_factor(sd[[cfg$subject_var]])
  }
  sample_data(ps) <- sample_data(sd)

  ps
}

as_ancombc_input <- function(ps) {
  if (
    requireNamespace("mia", quietly = TRUE) &&
      exists("convertFromPhyloseq", where = asNamespace("mia"), inherits = FALSE)
  ) {
    return(mia::convertFromPhyloseq(ps))
  }

  ps
}

extract_ancom_columns <- function(res_df, group_var) {
  lfc_col <- grep(paste0("^lfc_.*", group_var), names(res_df), value = TRUE)[1]
  se_col <- grep(paste0("^se_.*", group_var), names(res_df), value = TRUE)[1]
  p_col <- if (!is.na(lfc_col)) sub("^lfc_", "p_", lfc_col) else NA_character_
  q_col <- if (!is.na(lfc_col)) sub("^lfc_", "q_", lfc_col) else NA_character_

  cols <- c(lfc_col, se_col, p_col, q_col)
  if (any(is.na(cols)) || !all(cols %in% names(res_df))) {
    stop("Could not locate expected lfc/se/p/q columns for ", group_var)
  }

  list(lfc = lfc_col, se = se_col, p = p_col, q = q_col)
}

run_ancom_on_prepared <- function(ps_prepped, cfg) {
  ancom_input <- as_ancombc_input(ps_prepped)

  fit <- tryCatch(
    ancombc2(
      data = ancom_input,
      fix_formula = cfg$fix_formula,
      rand_formula = cfg$rand_formula,
      p_adj_method = "fdr",
      group = cfg$group_var,
      global = FALSE,
      pairwise = FALSE,
      dunnet = FALSE,
      trend = FALSE,
      struc_zero = TRUE,
      pseudo_sens = TRUE,
      prv_cut = prv_cut,
      lib_cut = lib_cut,
      alpha = alpha_cutoff,
      n_cl = n_cl,
      verbose = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(fit) || is.null(fit$res)) {
    return(NULL)
  }

  cols <- extract_ancom_columns(fit$res, cfg$group_var)

  out <- fit$res %>%
    transmute(
      taxon = .data$taxon,
      estimate = .data[[cols$lfc]],
      se = .data[[cols$se]],
      p_value = .data[[cols$p]],
      q_value = .data[[cols$q]],
      selected = .data[[cols$q]] < q_cutoff & abs(.data[[cols$lfc]]) >= lfc_thresh
    )

  out$taxon <- as.character(out$taxon)

  # Within-fit feature ranking by ANCOM-BC2 test statistic |W| = |lfc / se|.
  # Rank 1 = strongest evidence of differential abundance. Features that could
  # not be fit (non-finite statistic) receive NA and are excluded from ranking.
  out$W <- out$estimate / out$se
  rankable <- is.finite(out$W)
  out$rank <- NA_real_
  out$rank[rankable] <- rank(-abs(out$W[rankable]), ties.method = "average")
  out$n_ranked <- sum(rankable)

  out
}

run_ancom_once <- function(ps, cfg) {
  ps_prepped <- prepare_phyloseq(ps, cfg)
  run_ancom_on_prepared(ps_prepped, cfg)
}

summarize_bootstrap <- function(observed_df, boot_tables) {
  if (length(boot_tables) == 0) {
    observed_df$bootstrap_n <- 0L
    observed_df$selection_frequency <- NA_real_
    observed_df$bootstrap_median <- NA_real_
    observed_df$bootstrap_ci_low <- NA_real_
    observed_df$bootstrap_ci_high <- NA_real_
    observed_df$rank_n <- 0L
    observed_df$rank_median <- NA_real_
    observed_df$rank_ci_low <- NA_real_
    observed_df$rank_ci_high <- NA_real_
    return(observed_df)
  }

  boot_long <- bind_rows(boot_tables, .id = "bootstrap_id")

  boot_summary <- boot_long %>%
    group_by(.data$taxon) %>%
    summarise(
      bootstrap_n = sum(is.finite(.data$estimate)),
      selection_frequency = mean(.data$selected, na.rm = TRUE),
      bootstrap_median = median(.data$estimate, na.rm = TRUE),
      bootstrap_ci_low = as.numeric(quantile(.data$estimate, 0.025, na.rm = TRUE)),
      bootstrap_ci_high = as.numeric(quantile(.data$estimate, 0.975, na.rm = TRUE)),
      rank_n = sum(is.finite(.data$rank)),
      rank_median = median(.data$rank, na.rm = TRUE),
      rank_ci_low = as.numeric(quantile(.data$rank, 0.025, na.rm = TRUE)),
      rank_ci_high = as.numeric(quantile(.data$rank, 0.975, na.rm = TRUE)),
      .groups = "drop"
    )

  left_join(observed_df, boot_summary, by = "taxon")
}

summarize_permutations <- function(observed_df, perm_tables) {
  if (length(perm_tables) == 0) {
    observed_df$perm_n <- 0L
    observed_df$perm_p_value <- NA_real_
    observed_df$perm_q_value <- NA_real_
    return(list(
      observed = observed_df,
      discovery_summary = data.frame(
        observed_n_selected = sum(observed_df$selected, na.rm = TRUE),
        perm_mean_n_selected = NA_real_,
        perm_sd_n_selected = NA_real_,
        perm_empirical_p = NA_real_,
        stringsAsFactors = FALSE
      )
    ))
  }

  perm_long <- bind_rows(perm_tables, .id = "perm_id")
  perm_joined <- left_join(
    perm_long,
    observed_df %>% select(.data$taxon, observed_estimate = .data$estimate),
    by = "taxon"
  )

  perm_summary <- perm_joined %>%
    group_by(.data$taxon) %>%
    summarise(
      perm_n = sum(is.finite(.data$estimate)),
      perm_p_value = (1 + sum(abs(.data$estimate) >= abs(.data$observed_estimate), na.rm = TRUE)) /
        (1 + sum(is.finite(.data$estimate))),
      .groups = "drop"
    )

  observed_out <- left_join(observed_df, perm_summary, by = "taxon")
  observed_out$perm_q_value <- p.adjust(observed_out$perm_p_value, method = "BH")

  perm_disc <- perm_long %>%
    group_by(.data$perm_id) %>%
    summarise(n_selected = sum(.data$selected, na.rm = TRUE), .groups = "drop")

  observed_n_selected <- sum(observed_df$selected, na.rm = TRUE)
  discovery_summary <- data.frame(
    observed_n_selected = observed_n_selected,
    perm_mean_n_selected = mean(perm_disc$n_selected, na.rm = TRUE),
    perm_sd_n_selected = stats::sd(perm_disc$n_selected, na.rm = TRUE),
    perm_empirical_p = (1 + sum(perm_disc$n_selected >= observed_n_selected, na.rm = TRUE)) / (1 + nrow(perm_disc)),
    stringsAsFactors = FALSE
  )

  list(observed = observed_out, discovery_summary = discovery_summary)
}

run_one_analysis <- function(cfg) {
  ps_raw <- readRDS(cfg$rds_file)
  ps <- prepare_phyloseq(ps_raw, cfg)
  observed_df <- run_ancom_on_prepared(ps, cfg)
  if (is.null(observed_df)) {
    warning("ANCOM-BC2 failed for observed analysis: ", cfg$analysis)
    return(NULL)
  }

  boot_tables <- vector("list", n_boot)
  for (i in seq_len(n_boot)) {
    boot_ps <- if (cfg$paired) bootstrap_paired_phyloseq(ps, cfg) else bootstrap_unpaired_phyloseq(ps, cfg)
    boot_tables[[i]] <- run_ancom_on_prepared(boot_ps, cfg)
  }
  boot_tables <- boot_tables[!vapply(boot_tables, is.null, logical(1))]

  observed_df <- summarize_bootstrap(observed_df, boot_tables)
  discovery_summary <- NULL

  if (n_perm > 0L) {
    perm_tables <- vector("list", n_perm)
    for (i in seq_len(n_perm)) {
      perm_ps <- if (cfg$paired) permute_paired_phyloseq(ps, cfg) else permute_unpaired_phyloseq(ps, cfg)
      perm_tables[[i]] <- run_ancom_on_prepared(perm_ps, cfg)
    }
    perm_tables <- perm_tables[!vapply(perm_tables, is.null, logical(1))]

    perm_out <- summarize_permutations(observed_df, perm_tables)
    observed_df <- perm_out$observed
    discovery_summary <- perm_out$discovery_summary
    discovery_summary$analysis <- cfg$analysis
  }

  observed_df$analysis <- cfg$analysis
  observed_df <- observed_df %>%
    rename(observed_W = .data$W, observed_rank = .data$rank)

  # Explicit bootstrap-robustness classification. Rank is a relative cut: a
  # taxon must sit within the top fraction of testable (rankable) features for
  # this contrast, because the feature count (n_ranked) varies substantially
  # across contrasts/levels. n_ranked is the number of rankable features in the
  # observed fit, NOT the number of resamples.
  rank_robust_cut <- robust_rank_top_frac * observed_df$n_ranked
  observed_df <- observed_df %>%
    mutate(
      ci_excludes_zero = is.finite(.data$bootstrap_ci_low) &
        is.finite(.data$bootstrap_ci_high) &
        (.data$bootstrap_ci_low > 0 | .data$bootstrap_ci_high < 0),
      rank_robust = is.finite(.data$rank_median) &
        .data$rank_median <= rank_robust_cut,
      robust = .data$selected &
        is.finite(.data$selection_frequency) &
        .data$selection_frequency >= robust_selfreq_min &
        .data$ci_excludes_zero &
        .data$rank_robust
    )

  result_cols <- c(
    "analysis",
    "taxon",
    "estimate",
    "se",
    "p_value",
    "q_value",
    "selected",
    "observed_W",
    "observed_rank",
    "bootstrap_n",
    "selection_frequency",
    "bootstrap_median",
    "bootstrap_ci_low",
    "bootstrap_ci_high",
    "rank_n",
    "rank_median",
    "rank_ci_low",
    "rank_ci_high",
    "n_ranked",
    "ci_excludes_zero",
    "rank_robust",
    "robust"
  )
  if (n_perm > 0L) {
    result_cols <- c(result_cols, "perm_n", "perm_p_value", "perm_q_value")
  }
  observed_df <- observed_df %>% select(dplyr::all_of(result_cols))

  list(results = observed_df, discovery_summary = discovery_summary)
}


# MAIN

set.seed(seed)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat(
  "Running differential abundance robustness with",
  "n_boot =", n_boot,
  "n_perm =", n_perm,
  "out_dir =", out_dir,
  "\n"
)

for (cfg in analysis_configs) {
  cat("Running ANCOM-BC2 robustness for", cfg$analysis, "...\n")
  analysis_out <- run_one_analysis(cfg)
  if (is.null(analysis_out)) next

  utils::write.csv(
    analysis_out$results,
    file.path(out_dir, paste0(cfg$analysis, "_ancombc2_robustness.csv")),
    row.names = FALSE
  )

  if (!is.null(analysis_out$discovery_summary)) {
    utils::write.csv(
      analysis_out$discovery_summary,
      file.path(out_dir, paste0(cfg$analysis, "_ancombc2_permutation_discovery_summary.csv")),
      row.names = FALSE
    )
  }
}

cat("Differential abundance robustness outputs written to:", out_dir, "\n")
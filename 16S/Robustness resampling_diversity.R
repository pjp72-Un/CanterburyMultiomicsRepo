suppressPackageStartupMessages({
  library(phyloseq)
  library(microbiome)
  library(picante)
  library(dplyr)
  library(tidyr)
  library(lmerTest)
  library(vegan)
  library(permute)
})

options(stringsAsFactors = FALSE)


# USER CONFIG

seq_object_dir <- "C:/PhD/Sequencing Results/Files"
faith_tree_file <- "C:/PhD/Canterbury Dataset/objects/Cant_UF_tree.RDS"

ps1_untreated_file <- file.path(seq_object_dir, "Cant_untreated.rds")
ps1_treated_file <- file.path(seq_object_dir, "Cant_treated.rds")
ps1_paired_file <- file.path(seq_object_dir, "Cant_untreated_treated.rds")

ps2_all_file <- file.path(seq_object_dir, "Cant_5%.rds")
ps2_untreated_file <- file.path(seq_object_dir, "Cant_5%_untreated.rds")
ps2_treated_file <- file.path(seq_object_dir, "Cant_5%_treated.rds")
ps2_paired_file <- file.path(seq_object_dir, "Cant_5%_untreated_treated.rds")

out_dir <- "C:/PhD/Canterbury Dataset/16S/robustness"

n_boot <- 500
n_perm <- 4999
seed <- 27

alpha_pseudocount <- 1
clr_pseudocount <- 1e-6

selected_phyla <- c(
  "Verrucomicrobia",
  "Actinobacteria",
  "Proteobacteria",
  "Firmicutes",
  "Bacteroidetes"
)
top_n_phyla <- 10


# HELPERS

safe_factor <- function(x, levels = NULL) {
  x <- as.character(x)
  if (is.null(levels)) {
    factor(x)
  } else {
    factor(x, levels = levels)
  }
}

as_flag <- function(x) {
  if (is.logical(x)) return(x)

  if (is.numeric(x)) {
    out <- rep(NA, length(x))
    out[!is.na(x)] <- x[!is.na(x)] != 0
    return(out)
  }

  x_chr <- trimws(tolower(as.character(x)))
  out <- rep(NA, length(x_chr))
  out[x_chr %in% c("true", "t", "1", "yes", "y")] <- TRUE
  out[x_chr %in% c("false", "f", "0", "no", "n")] <- FALSE
  out
}

extract_lm_term <- function(fit, term) {
  coefs <- summary(fit)$coefficients
  if (!(term %in% rownames(coefs))) {
    return(c(estimate = NA_real_, p_value = NA_real_))
  }

  p_col <- grep("^Pr\\(", colnames(coefs), value = TRUE)[1]
  if (is.na(p_col)) p_col <- "Pr(>|t|)"

  c(
    estimate = unname(coefs[term, "Estimate"]),
    p_value = if (p_col %in% colnames(coefs)) unname(coefs[term, p_col]) else NA_real_
  )
}

extract_lmer_term <- function(fit, term) {
  coefs <- summary(fit)$coefficients
  if (!(term %in% rownames(coefs))) {
    return(c(estimate = NA_real_, p_value = NA_real_))
  }

  if ("Pr(>|t|)" %in% colnames(coefs)) {
    p_value <- unname(coefs[term, "Pr(>|t|)"])
  } else if ("t value" %in% colnames(coefs)) {
    p_value <- 2 * pnorm(-abs(unname(coefs[term, "t value"])))
  } else {
    p_value <- NA_real_
  }

  c(
    estimate = unname(coefs[term, "Estimate"]),
    p_value = p_value
  )
}

fit_lm_effect <- function(df, formula, term) {
  fit <- try(stats::lm(formula, data = df, na.action = stats::na.omit), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(c(estimate = NA_real_, p_value = NA_real_))
  }
  extract_lm_term(fit, term)
}

fit_lmer_effect <- function(df, formula, term, reml = FALSE) {
  fit <- try(lmerTest::lmer(formula, data = df, REML = reml, na.action = stats::na.omit), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(c(estimate = NA_real_, p_value = NA_real_))
  }
  extract_lmer_term(fit, term)
}

resample_clusters <- function(df, cluster_var, strata_var = NULL) {
  out <- vector("list", 0)
  next_id <- 1L

  if (is.null(strata_var)) {
    cluster_ids <- unique(as.character(df[[cluster_var]]))
    draws <- sample(cluster_ids, length(cluster_ids), replace = TRUE)

    for (draw_id in draws) {
      part <- df[as.character(df[[cluster_var]]) == draw_id, , drop = FALSE]
      part[[cluster_var]] <- paste0(draw_id, "__boot", next_id)
      out[[length(out) + 1L]] <- part
      next_id <- next_id + 1L
    }
  } else {
    strata_levels <- unique(as.character(df[[strata_var]]))

    for (stratum in strata_levels) {
      sub_df <- df[as.character(df[[strata_var]]) == stratum, , drop = FALSE]
      cluster_ids <- unique(as.character(sub_df[[cluster_var]]))
      draws <- sample(cluster_ids, length(cluster_ids), replace = TRUE)

      for (draw_id in draws) {
        part <- sub_df[as.character(sub_df[[cluster_var]]) == draw_id, , drop = FALSE]
        part[[cluster_var]] <- paste0(draw_id, "__boot", next_id)
        out[[length(out) + 1L]] <- part
        next_id <- next_id + 1L
      }
    }
  }

  out_df <- do.call(rbind, out)
  rownames(out_df) <- NULL
  out_df
}

permute_group_by_cluster <- function(df, group_var, cluster_var) {
  out <- df
  cluster_df <- unique(out[, c(cluster_var, group_var), drop = FALSE])
  cluster_ids <- as.character(cluster_df[[cluster_var]])
  permuted_groups <- sample(as.character(cluster_df[[group_var]]), length(cluster_ids), replace = FALSE)
  names(permuted_groups) <- cluster_ids

  levs <- levels(out[[group_var]])
  out[[group_var]] <- factor(permuted_groups[as.character(out[[cluster_var]])], levels = levs)
  out
}

flip_binary_within_cluster <- function(df, group_var, cluster_var) {
  out <- df
  levs <- levels(out[[group_var]])
  if (length(levs) != 2L) {
    stop("flip_binary_within_cluster() requires exactly two levels.")
  }

  cluster_ids <- unique(as.character(out[[cluster_var]]))
  flip_ids <- sample(c(FALSE, TRUE), length(cluster_ids), replace = TRUE)
  names(flip_ids) <- cluster_ids

  current <- as.character(out[[group_var]])
  to_flip <- flip_ids[as.character(out[[cluster_var]])]
  current[to_flip & current == levs[1]] <- levs[2]
  current[to_flip & current == levs[2]] <- levs[1]
  out[[group_var]] <- factor(current, levels = levs)
  out
}

collect_bootstrap_estimates <- function(df, fit_fun, resample_fun, n_iter) {
  estimates <- rep(NA_real_, n_iter)
  for (i in seq_len(n_iter)) {
    boot_df <- resample_fun(df)
    estimates[i] <- unname(fit_fun(boot_df)["estimate"])
  }
  estimates[is.finite(estimates)]
}

collect_permutation_estimates <- function(df, fit_fun, permute_fun, n_iter) {
  estimates <- rep(NA_real_, n_iter)
  for (i in seq_len(n_iter)) {
    perm_df <- permute_fun(df)
    estimates[i] <- unname(fit_fun(perm_df)["estimate"])
  }
  estimates[is.finite(estimates)]
}

summarize_resampling <- function(analysis, comparison, feature, observed, boot_est, perm_est) {
  observed_est <- unname(observed["estimate"])
  observed_p <- unname(observed["p_value"])

  if (length(boot_est) == 0) {
    boot_ci_low <- NA_real_
    boot_ci_high <- NA_real_
    boot_median <- NA_real_
    boot_sign_consistency <- NA_real_
  } else {
    boot_ci_low <- as.numeric(stats::quantile(boot_est, 0.025, na.rm = TRUE))
    boot_ci_high <- as.numeric(stats::quantile(boot_est, 0.975, na.rm = TRUE))
    boot_median <- stats::median(boot_est, na.rm = TRUE)
    if (is.finite(observed_est) && observed_est != 0) {
      boot_sign_consistency <- mean(sign(boot_est) == sign(observed_est), na.rm = TRUE)
    } else {
      boot_sign_consistency <- NA_real_
    }
  }

  if (!is.finite(observed_est) || length(perm_est) == 0) {
    permutation_p <- NA_real_
  } else {
    permutation_p <- (1 + sum(abs(perm_est) >= abs(observed_est), na.rm = TRUE)) / (length(perm_est) + 1)
  }

  data.frame(
    analysis = analysis,
    comparison = comparison,
    feature = feature,
    observed_estimate = observed_est,
    observed_p = observed_p,
    observed_q = NA_real_,
    bootstrap_n = length(boot_est),
    bootstrap_median = boot_median,
    bootstrap_ci_low = boot_ci_low,
    bootstrap_ci_high = boot_ci_high,
    bootstrap_sign_consistency = boot_sign_consistency,
    permutation_n = length(perm_est),
    permutation_p = permutation_p,
    permutation_q = NA_real_,
    stringsAsFactors = FALSE
  )
}

add_q_values <- function(df) {
  if (nrow(df) == 0) return(df)

  df %>%
    group_by(.data$analysis, .data$comparison) %>%
    mutate(
      observed_q = p.adjust(.data$observed_p, method = "BH"),
      permutation_q = p.adjust(.data$permutation_p, method = "BH")
    ) %>%
    ungroup()
}

metric_mapping <- c(
  observed = "observed",
  shannon = "diversity_shannon",
  fisher = "diversity_fisher",
  bulla = "evenness_bulla",
  simpson = "evenness_simpson",
  dominance = "dominance_simpson",
  rare_abundance = "rarity_low_abundance"
)

compute_alpha_table <- function(ps) {
  alpha_df <- microbiome::alpha(
    ps,
    index = c(
      "observed", "shannon", "fisher", "bulla",
      "simpson", "dominance", "low_abundance", "rare_abundance"
    )
  )
  alpha_df <- data.frame(
    SampleID = rownames(alpha_df),
    alpha_df,
    row.names = NULL,
    check.names = FALSE
  )

  meta_df <- as(sample_data(ps), "data.frame")
  meta_df <- data.frame(
    SampleID = rownames(meta_df),
    meta_df,
    row.names = NULL,
    check.names = FALSE
  )
  meta_df$read_depth <- sample_sums(ps)

  merged <- dplyr::left_join(alpha_df, meta_df, by = "SampleID")
  merged$age <- suppressWarnings(as.numeric(merged$Age))
  merged$sex <- safe_factor(merged$Sex)
  merged$subject <- safe_factor(merged$Subject)
  merged
}

run_alpha_robustness <- function(ps, comparison, paired = FALSE) {
  merged <- compute_alpha_table(ps)
  rows <- list()
  next_row <- 1L

  if (paired) {
    merged$state <- safe_factor(ifelse(as_flag(merged$GFD), "Treated", "untreated"), levels = c("untreated", "Treated"))
    fit_formula <- stats::as.formula("outcome ~ state + read_depth + (1 | subject)")
    term_name <- "stateTreated"
    fit_fun <- function(df) fit_lmer_effect(df, fit_formula, term_name)
    boot_fun <- function(df) resample_clusters(df, cluster_var = "subject")
  } else {
    merged$group <- safe_factor(ifelse(as_flag(merged$Diagnosed.CD), "CeD", "Control"), levels = c("Control", "CeD"))
    fit_formula <- stats::as.formula("outcome ~ group + read_depth + age + sex")
    term_name <- "groupCeD"
    fit_fun <- function(df) fit_lm_effect(df, fit_formula, term_name)
    boot_fun <- function(df) resample_clusters(df, cluster_var = "subject", strata_var = "group")
  }

  for (metric_name in names(metric_mapping)) {
    metric_col <- metric_mapping[[metric_name]]
    if (!(metric_col %in% colnames(merged))) next

    work_df <- merged
    work_df$outcome <- log2(work_df[[metric_col]] + alpha_pseudocount)

    observed <- fit_fun(work_df)
    boot_est <- collect_bootstrap_estimates(work_df, fit_fun, boot_fun, n_boot)
    perm_est <- numeric(0)

    rows[[next_row]] <- summarize_resampling(
      analysis = "alpha_diversity",
      comparison = comparison,
      feature = metric_name,
      observed = observed,
      boot_est = boot_est,
      perm_est = perm_est
    )
    next_row <- next_row + 1L
  }

  add_q_values(bind_rows(rows))
}

compute_faith_pd_df <- function(ps) {
  comm <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) comm <- t(comm)

  pd_res <- picante::pd(comm, phy_tree(ps), include.root = FALSE)
  meta_df <- as(sample_data(ps), "data.frame")
  meta_df <- data.frame(
    SampleID = rownames(meta_df),
    meta_df,
    row.names = NULL,
    check.names = FALSE
  )
  meta_df$FaithPD <- pd_res$PD[match(meta_df$SampleID, rownames(pd_res))]
  meta_df$age <- suppressWarnings(as.numeric(meta_df$Age))
  meta_df$sex <- safe_factor(meta_df$Sex)
  meta_df$subject <- safe_factor(meta_df$Subject)
  meta_df
}

build_faith_objects <- function(ps) {
  meta_df <- as.data.frame(sample_data(ps))

  ps_untreated <- prune_samples(
    rownames(meta_df)[meta_df$GFD == FALSE & !is.na(meta_df$Diagnosed.CD)],
    ps
  )
  ps_untreated <- prune_taxa(taxa_sums(ps_untreated) > 0, ps_untreated)

  ps_treated <- prune_samples(
    rownames(meta_df)[
      (meta_df$GFD == TRUE & meta_df$Diagnosed.CD == TRUE) |
        (meta_df$GFD == FALSE & meta_df$Diagnosed.CD == FALSE)
    ],
    ps
  )
  ps_treated <- prune_taxa(taxa_sums(ps_treated) > 0, ps_treated)

  is_untreated_cd <- with(meta_df, Diagnosed.CD == TRUE & GFD == FALSE & !is.na(Subject))
  is_treated_cd <- with(meta_df, Diagnosed.CD == TRUE & GFD == TRUE & !is.na(Subject))
  paired_subjects <- intersect(unique(meta_df$Subject[is_untreated_cd]), unique(meta_df$Subject[is_treated_cd]))

  paired_samples <- rownames(meta_df)[
    meta_df$Subject %in% paired_subjects &
      meta_df$Diagnosed.CD == TRUE &
      meta_df$GFD %in% c(FALSE, TRUE)
  ]

  ps_paired <- prune_samples(paired_samples, ps)
  ps_paired <- prune_taxa(taxa_sums(ps_paired) > 0, ps_paired)

  list(
    untreated = ps_untreated,
    treated = ps_treated,
    paired = ps_paired
  )
}

run_faith_robustness <- function(ps_tree) {
  faith_objects <- build_faith_objects(ps_tree)
  rows <- list()
  next_row <- 1L

  for (comparison in c("untreated_vs_control", "treated_vs_control")) {
    faith_df <- compute_faith_pd_df(faith_objects[[sub("_vs_control", "", comparison)]])
    faith_df$group <- safe_factor(ifelse(as_flag(faith_df$Diagnosed.CD), "CeD", "Control"), levels = c("Control", "CeD"))

    fit_formula <- stats::as.formula("FaithPD ~ group + age + sex")
    fit_fun <- function(df) fit_lm_effect(df, fit_formula, "groupCeD")
    boot_fun <- function(df) resample_clusters(df, cluster_var = "subject", strata_var = "group")

    observed <- fit_fun(faith_df)
    boot_est <- collect_bootstrap_estimates(faith_df, fit_fun, boot_fun, n_boot)
    perm_est <- numeric(0)

    rows[[next_row]] <- summarize_resampling(
      analysis = "faith_pd",
      comparison = comparison,
      feature = "FaithPD",
      observed = observed,
      boot_est = boot_est,
      perm_est = perm_est
    )
    next_row <- next_row + 1L
  }

  faith_pair_df <- compute_faith_pd_df(faith_objects$paired)
  faith_pair_df$state <- safe_factor(ifelse(as_flag(faith_pair_df$GFD), "Treated", "untreated"), levels = c("untreated", "Treated"))

  fit_formula <- stats::as.formula("FaithPD ~ state + age + (1 | subject)")
  fit_fun <- function(df) fit_lmer_effect(df, fit_formula, "stateTreated")
  boot_fun <- function(df) resample_clusters(df, cluster_var = "subject")

  observed <- fit_fun(faith_pair_df)
  boot_est <- collect_bootstrap_estimates(faith_pair_df, fit_fun, boot_fun, n_boot)
  perm_est <- numeric(0)

  rows[[next_row]] <- summarize_resampling(
    analysis = "faith_pd",
    comparison = "paired_treated_vs_untreated",
    feature = "FaithPD",
    observed = observed,
    boot_est = boot_est,
    perm_est = perm_est
  )

  add_q_values(bind_rows(rows))
}

get_phylum_rel_abund_df <- function(ps) {
  ps_ph <- tax_glom(ps, taxrank = "Phylum")
  ps_ph <- prune_taxa(taxa_sums(ps_ph) > 0, ps_ph)
  ps_ph_rel <- transform_sample_counts(ps_ph, function(x) x / sum(x))

  comm <- as.data.frame(otu_table(ps_ph_rel))
  if (taxa_are_rows(ps_ph_rel)) comm <- t(comm)
  comm <- as.data.frame(comm, check.names = FALSE)

  tax_df <- as.data.frame(tax_table(ps_ph_rel))
  phylum_names <- as.character(tax_df$Phylum)
  phylum_names[is.na(phylum_names) | phylum_names == ""] <- "Unidentified"
  taxa_names(ps_ph_rel) <- make.unique(phylum_names)

  comm <- as.data.frame(otu_table(ps_ph_rel))
  if (taxa_are_rows(ps_ph_rel)) comm <- t(comm)
  comm <- as.data.frame(comm, check.names = FALSE)

  comm$SampleID <- rownames(comm)
  long_df <- tidyr::pivot_longer(comm, cols = -all_of("SampleID"), names_to = "Phylum", values_to = "Abundance")

  meta_df <- as.data.frame(sample_data(ps_ph_rel))
  meta_df$SampleID <- rownames(meta_df)

  merged <- dplyr::left_join(long_df, meta_df, by = "SampleID")
  merged$age <- suppressWarnings(as.numeric(merged$Age))
  merged$sex <- safe_factor(merged$Sex)
  merged$subject <- safe_factor(merged$Subject)
  merged
}

pick_phyla <- function(paired_df) {
  if (!is.null(selected_phyla) && length(selected_phyla) > 0) {
    keep <- intersect(selected_phyla, unique(paired_df$Phylum))
    if (length(keep) > 0) return(keep)
  }

  paired_df %>%
    group_by(.data$Phylum) %>%
    summarise(mean_abund = mean(.data$Abundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(.data$mean_abund)) %>%
    slice_head(n = top_n_phyla) %>%
    pull(.data$Phylum)
}

run_phylum_robustness <- function(ps_untreated, ps_treated, ps_paired) {
  # Match 'Phylum abundance analysis.R': restrict each object to the selected
  # phyla BEFORE computing relative abundance, so abundances are normalized
  # within the selected phyla rather than the whole community.
  subset_to_selected_phyla <- function(ps) {
    tt <- as.data.frame(tax_table(ps))
    if (!("Phylum" %in% colnames(tt)) || is.null(selected_phyla) || length(selected_phyla) == 0) {
      return(ps)
    }
    keep_taxa <- rownames(tt)[as.character(tt$Phylum) %in% selected_phyla]
    ps_sub <- prune_taxa(keep_taxa, ps)
    prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  }

  ps_untreated <- subset_to_selected_phyla(ps_untreated)
  ps_treated <- subset_to_selected_phyla(ps_treated)
  ps_paired <- subset_to_selected_phyla(ps_paired)

  untreated_df <- get_phylum_rel_abund_df(ps_untreated)
  treated_df <- get_phylum_rel_abund_df(ps_treated)
  paired_df <- get_phylum_rel_abund_df(ps_paired)
  phyla_to_test <- pick_phyla(paired_df)

  rows <- list()
  next_row <- 1L

  run_one <- function(df, comparison, paired = FALSE) {
    local_rows <- list()
    local_next <- 1L

    if (paired) {
      df$state <- safe_factor(ifelse(as_flag(df$GFD), "Treated", "untreated"), levels = c("untreated", "Treated"))
      fit_formula <- stats::as.formula("outcome ~ state + age + (1 | subject)")
      term_name <- "stateTreated"
      # REML = TRUE to match the default lmer() fit in 'Phylum abundance analysis.R'
      fit_fun <- function(x) fit_lmer_effect(x, fit_formula, term_name, reml = TRUE)
      boot_fun <- function(x) resample_clusters(x, cluster_var = "subject")
      perm_fun <- function(x) flip_binary_within_cluster(x, group_var = "state", cluster_var = "subject")
    } else {
      df$group <- safe_factor(ifelse(as_flag(df$Diagnosed.CD), "CeD", "Control"), levels = c("Control", "CeD"))
      fit_formula <- stats::as.formula("outcome ~ group + age + sex")
      term_name <- "groupCeD"
      fit_fun <- function(x) fit_lm_effect(x, fit_formula, term_name)
      boot_fun <- function(x) resample_clusters(x, cluster_var = "subject", strata_var = "group")
      perm_fun <- function(x) permute_group_by_cluster(x, group_var = "group", cluster_var = "subject")
    }

    for (phylum_name in phyla_to_test) {
      work_df <- df[df$Phylum == phylum_name, , drop = FALSE]
      if (nrow(work_df) == 0) next

      work_df$outcome <- asin(sqrt(work_df$Abundance))
      observed <- fit_fun(work_df)
      boot_est <- collect_bootstrap_estimates(work_df, fit_fun, boot_fun, n_boot)
      perm_est <- collect_permutation_estimates(work_df, fit_fun, perm_fun, n_perm)

      local_rows[[local_next]] <- summarize_resampling(
        analysis = "phylum_abundance",
        comparison = comparison,
        feature = phylum_name,
        observed = observed,
        boot_est = boot_est,
        perm_est = perm_est
      )
      local_next <- local_next + 1L
    }

    bind_rows(local_rows)
  }

  rows[[next_row]] <- run_one(untreated_df, "untreated_vs_control", paired = FALSE)
  next_row <- next_row + 1L
  rows[[next_row]] <- run_one(treated_df, "treated_vs_control", paired = FALSE)
  next_row <- next_row + 1L
  rows[[next_row]] <- run_one(paired_df, "paired_treated_vs_untreated", paired = TRUE)

  add_q_values(bind_rows(rows))
}

make_group_meta <- function(ps) {
  meta_df <- as(sample_data(ps), "data.frame")
  meta_df$SampleID <- rownames(meta_df)
  meta_df$Group <- dplyr::case_when(
    meta_df$Diagnosed.CD == FALSE ~ "Control",
    meta_df$Diagnosed.CD == TRUE & meta_df$GFD == FALSE ~ "untreated",
    meta_df$Diagnosed.CD == TRUE & meta_df$GFD == TRUE ~ "Treated",
    TRUE ~ NA_character_
  )

  meta_df$Subject <- safe_factor(meta_df$Subject)
  meta_df$Group <- safe_factor(meta_df$Group, levels = c("Control", "untreated", "Treated"))
  meta_df$Sex <- safe_factor(meta_df$Sex)
  meta_df$Age <- suppressWarnings(as.numeric(meta_df$Age))
  meta_df
}

compute_aitchison_distance <- function(ps) {
  otu_mat <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) otu_mat <- t(otu_mat)
  otu_pseudo <- otu_mat + clr_pseudocount
  clr_mat <- t(apply(otu_pseudo, 1, function(x) log(x) - mean(log(x))))
  stats::dist(clr_mat, method = "euclidean")
}

summarize_adonis <- function(adonis_obj, comparison, test_name) {
  out <- as.data.frame(adonis_obj)
  out$term <- rownames(out)
  rownames(out) <- NULL
  out <- out[!(out$term %in% c("Residual", "Total")), , drop = FALSE]

  data.frame(
    comparison = comparison,
    test = test_name,
    term = out$term,
    df = out$Df,
    sum_of_squares = out$SumOfSqs,
    mean_square = out$SumOfSqs / out$Df,
    r2 = out$R2,
    f_statistic = out$F,
    p_value = out[["Pr(>F)"]],
    n_perm = rep(n_perm, nrow(out)),
    stringsAsFactors = FALSE
  )
}

run_betadisper_summary <- function(dist_obj, groups, comparison) {
  if (length(unique(groups)) < 2L) {
    return(data.frame(
      comparison = comparison,
      test = "betadisper",
      term = "Group",
      df = NA_real_,
      sum_of_squares = NA_real_,
      mean_square = NA_real_,
      r2 = NA_real_,
      f_statistic = NA_real_,
      p_value = NA_real_,
      n_perm = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  bd <- vegan::betadisper(dist_obj, groups)
  pt <- vegan::permutest(bd, permutations = n_perm)
  tab <- as.data.frame(pt$tab)
  tab$term <- rownames(tab)
  rownames(tab) <- NULL
  tab <- tab[tab$term == "Groups", , drop = FALSE]

  data.frame(
    comparison = comparison,
    test = "betadisper",
    term = "Group",
    df = tab$Df,
    sum_of_squares = tab[["Sum Sq"]],
    mean_square = tab[["Mean Sq"]],
    r2 = NA_real_,
    f_statistic = tab$F,
    p_value = tab[["Pr(>F)"]],
    n_perm = tab$N.Perm,
    stringsAsFactors = FALSE
  )
}

build_paired_subject_permutations <- function(subjects, n_perm) {
  subject_ids <- as.character(subjects)
  block_index <- split(seq_along(subject_ids), subject_ids)
  block_sizes <- vapply(block_index, length, integer(1))

  if (length(block_index) == 0L || any(block_sizes != 2L)) {
    stop("Paired beta-diversity permutations require exactly two samples per subject.")
  }

  n_blocks <- length(block_index)
  total_perm <- (2^n_blocks) - 1L
  if (total_perm <= 0L) {
    return(matrix(integer(0), nrow = 0L, ncol = length(subject_ids)))
  }

  make_perm_matrix <- function(combo_mat) {
    perm_mat <- matrix(
      rep(seq_along(subject_ids), each = nrow(combo_mat)),
      nrow = nrow(combo_mat)
    )

    for (j in seq_len(n_blocks)) {
      idx <- block_index[[j]]
      swap_rows <- combo_mat[, j]
      perm_mat[swap_rows, idx[1]] <- idx[2]
      perm_mat[swap_rows, idx[2]] <- idx[1]
    }

    perm_mat
  }

  if (total_perm <= n_perm) {
    combo_mat <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), n_blocks), KEEP.OUT.ATTRS = FALSE))
    combo_mat <- combo_mat[rowSums(combo_mat) > 0, , drop = FALSE]
    return(make_perm_matrix(combo_mat))
  }

  target_n <- min(n_perm, total_perm)
  combo_keys <- character(0)
  combo_rows <- vector("list", 0)

  while (length(combo_rows) < target_n) {
    combo <- sample(c(FALSE, TRUE), n_blocks, replace = TRUE)
    if (!any(combo)) next

    combo_key <- paste(as.integer(combo), collapse = "")
    if (combo_key %in% combo_keys) next

    combo_keys <- c(combo_keys, combo_key)
    combo_rows[[length(combo_rows) + 1L]] <- combo
  }

  combo_mat <- do.call(rbind, combo_rows)
  make_perm_matrix(combo_mat)
}

run_beta_robustness <- function(ps) {
  dist_all <- compute_aitchison_distance(ps)
  meta_df <- make_group_meta(ps)

  rows <- list()
  next_row <- 1L

  adonis_global <- vegan::adonis2(
    dist_all ~ Group + Sex + Age,
    data = meta_df,
    permutations = n_perm
  )
  rows[[next_row]] <- summarize_adonis(adonis_global, "all_groups", "adonis2")
  next_row <- next_row + 1L
  rows[[next_row]] <- run_betadisper_summary(dist_all, meta_df$Group, "all_groups")
  next_row <- next_row + 1L

  paired_ids <- meta_df %>%
    filter(.data$Group %in% c("untreated", "Treated"), !is.na(.data$Subject)) %>%
    group_by(.data$Subject) %>%
    filter(n_distinct(.data$Group) == 2) %>%
    ungroup() %>%
    pull(.data$SampleID)

  if (length(paired_ids) > 1) {
    meta_paired <- meta_df[match(paired_ids, meta_df$SampleID), , drop = FALSE]
    dist_paired <- as.dist(as.matrix(dist_all)[paired_ids, paired_ids])
    perms <- build_paired_subject_permutations(meta_paired$Subject, n_perm = n_perm)

    adonis_paired <- vegan::adonis2(
      dist_paired ~ Group + Age,
      data = meta_paired,
      permutations = perms
    )
    rows[[next_row]] <- summarize_adonis(adonis_paired, "paired_treated_vs_untreated", "adonis2")
    next_row <- next_row + 1L
    rows[[next_row]] <- run_betadisper_summary(dist_paired, meta_paired$Group, "paired_treated_vs_untreated")
    next_row <- next_row + 1L
  }

  untreated_control_ids <- meta_df$SampleID[meta_df$Group %in% c("untreated", "Control")]
  if (length(untreated_control_ids) > 1) {
    meta_uc <- meta_df[match(untreated_control_ids, meta_df$SampleID), , drop = FALSE]
    dist_uc <- as.dist(as.matrix(dist_all)[untreated_control_ids, untreated_control_ids])
    adonis_uc <- vegan::adonis2(
      dist_uc ~ Group + Sex + Age,
      data = meta_uc,
      permutations = n_perm
    )
    rows[[next_row]] <- summarize_adonis(adonis_uc, "untreated_vs_control", "adonis2")
    next_row <- next_row + 1L
    rows[[next_row]] <- run_betadisper_summary(dist_uc, meta_uc$Group, "untreated_vs_control")
    next_row <- next_row + 1L
  }

  treated_control_ids <- meta_df$SampleID[meta_df$Group %in% c("Treated", "Control")]
  if (length(treated_control_ids) > 1) {
    meta_tc <- meta_df[match(treated_control_ids, meta_df$SampleID), , drop = FALSE]
    dist_tc <- as.dist(as.matrix(dist_all)[treated_control_ids, treated_control_ids])
    adonis_tc <- vegan::adonis2(
      dist_tc ~ Group + Sex + Age,
      data = meta_tc,
      permutations = n_perm
    )
    rows[[next_row]] <- summarize_adonis(adonis_tc, "treated_vs_control", "adonis2")
    next_row <- next_row + 1L
    rows[[next_row]] <- run_betadisper_summary(dist_tc, meta_tc$Group, "treated_vs_control")
  }

  bind_rows(rows)
}


# MAIN

set.seed(seed)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

alpha_results <- bind_rows(
  run_alpha_robustness(readRDS(ps1_untreated_file), "untreated_vs_control", paired = FALSE),
  run_alpha_robustness(readRDS(ps1_treated_file), "treated_vs_control", paired = FALSE),
  run_alpha_robustness(readRDS(ps1_paired_file), "paired_treated_vs_untreated", paired = TRUE)
)

write.csv(
  alpha_results,
  file = file.path(out_dir, "alpha_diversity_robustness.csv"),
  row.names = FALSE
)

if (file.exists(faith_tree_file)) {
  faith_results <- run_faith_robustness(readRDS(faith_tree_file))
  write.csv(
    faith_results,
    file = file.path(out_dir, "faith_pd_robustness.csv"),
    row.names = FALSE
  )
} else {
  warning("Faith's PD tree file not found; skipping faith_pd_robustness.csv")
}

beta_results <- run_beta_robustness(readRDS(ps2_all_file))

write.csv(
  beta_results,
  file = file.path(out_dir, "beta_diversity_robustness_summary.csv"),
  row.names = FALSE
)

phylum_results <- run_phylum_robustness(
  ps_untreated = readRDS(ps2_untreated_file),
  ps_treated   = readRDS(ps2_treated_file),
  ps_paired    = readRDS(ps2_paired_file)
)

write.csv(
  phylum_results,
  file = file.path(out_dir, "phylum_abundance_robustness.csv"),
  row.names = FALSE
)

cat("16S robustness outputs written to:", out_dir, "\n")
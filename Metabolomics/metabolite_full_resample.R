suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(lmerTest)
})

options(stringsAsFactors = FALSE)


# USER CONFIG

nmr_data_file <- "C:/PhD/R code/nmr/data_CoeliacStoolCHCH_pqn_selected_bins.tsv"
nmr_meta_file <- "C:/PhD/NMR/metablo_meta.txt"
gcms_data_file <- "C:/PhD/Canterbury Dataset/GC-MS/GC_MS_metab.txt"
gcms_meta_file <- "C:/PhD/NMR/metablo_meta.txt"

out_dir <- "C:/PhD/Metabolomics/robustness"

n_boot <- 500
seed <- 27
q_cutoff <- 0.10

nmr_transform_mode <- "log2_shift"
gcms_transform_mode <- "log2_shift"

nmr_bins_keep <- c(
  "B1_0398", "B1_5414", "B2_6350", "B2_3582", "B2_7966",
  "B1_7458", "B2_2856", "B3_1092", "B1_4742", "B0_9539", "B1_0209",
  "B0_9982", "B1_9199", "B2_0511", "B3_9598"
)

gcms_features_keep <- NULL


# HELPERS

# GC-MS transform: matches log2_shift() in 'GCMS lm.R'
log2_shift <- function(x) {
  x <- as.numeric(x)
  xmin <- min(x, na.rm = TRUE)
  if (!is.finite(xmin)) return(x)

  x_shift <- x - xmin
  pos <- x_shift[is.finite(x_shift) & x_shift > 0]
  pseudocount <- if (length(pos)) min(pos, na.rm = TRUE) / 2 else 1e-6

  log2(x_shift + pseudocount)
}

# NMR transform: matches log2_shift() in 'NMR lm.R'
nmr_log2_shift <- function(x) {
  x <- as.numeric(x)
  xmin <- min(x, na.rm = TRUE)
  shift <- ifelse(xmin <= 0, -xmin + 1, 0)
  log2(x + shift)
}

safe_lmer <- function(formula, data) {
  tryCatch(lmerTest::lmer(formula, data = data, REML = FALSE), error = function(e) NULL)
}

safe_lm <- function(formula, data) {
  tryCatch(stats::lm(formula, data = data), error = function(e) NULL)
}

extract_terms <- function(fit, terms) {
  out <- data.frame(
    term = terms,
    estimate = NA_real_,
    se = NA_real_,
    statistic = NA_real_,
    p_value = NA_real_,
    stringsAsFactors = FALSE
  )

  if (is.null(fit)) return(out)

  coefs <- summary(fit)$coefficients

  # Match single-metabolite scripts: normal-approximation z-test on
  # estimate / std.error (NOT the Satterthwaite Pr(>|t|) from lmerTest).
  for (i in seq_along(terms)) {
    term <- terms[i]
    if (!(term %in% rownames(coefs))) next

    est <- as.numeric(coefs[term, "Estimate"])
    se <- as.numeric(coefs[term, "Std. Error"])
    stat <- est / se

    out$estimate[i] <- est
    out$se[i] <- se
    out$statistic[i] <- stat
    out$p_value[i] <- 2 * pnorm(-abs(stat))
  }

  out
}

rename_with_suffix <- function(df, vars, suffix) {
  for (var in vars) {
    if (var %in% colnames(df)) {
      df[[var]] <- paste0(as.character(df[[var]]), suffix)
    }
  }

  df
}

bootstrap_cluster_df <- function(df, cluster_var, strata_var = NULL, rename_vars = cluster_var) {
  out <- vector("list", 0)
  next_id <- 1L

  if (is.null(strata_var)) {
    cluster_ids <- unique(as.character(df[[cluster_var]]))
    draw <- sample(cluster_ids, length(cluster_ids), replace = TRUE)

    for (cid in draw) {
      part <- df[as.character(df[[cluster_var]]) == cid, , drop = FALSE]
      part <- rename_with_suffix(part, rename_vars, paste0("__boot", next_id))
      out[[length(out) + 1L]] <- part
      next_id <- next_id + 1L
    }
  } else {
    strata_levels <- unique(as.character(df[[strata_var]]))

    for (stratum in strata_levels) {
      sub_df <- df[as.character(df[[strata_var]]) == stratum, , drop = FALSE]
      cluster_ids <- unique(as.character(sub_df[[cluster_var]]))
      draw <- sample(cluster_ids, length(cluster_ids), replace = TRUE)

      for (cid in draw) {
        part <- sub_df[as.character(sub_df[[cluster_var]]) == cid, , drop = FALSE]
        part <- rename_with_suffix(part, rename_vars, paste0("__boot", next_id))
        out[[length(out) + 1L]] <- part
        next_id <- next_id + 1L
      }
    }
  }

  out_df <- bind_rows(out)
  rownames(out_df) <- NULL
  out_df
}

summarize_resampling <- function(observed_df, boot_df) {
  boot_summary <- if (nrow(boot_df) == 0) {
    observed_df %>%
      transmute(
        term = .data$term,
        bootstrap_n = 0L,
        bootstrap_mean = NA_real_,
        bootstrap_sd = NA_real_,
        bootstrap_ci_low = NA_real_,
        bootstrap_ci_high = NA_real_,
        bootstrap_sign_consistency = NA_real_
      )
  } else {
    left_join(
      boot_df,
      observed_df %>% transmute(term = .data$term, observed_estimate = .data$estimate),
      by = "term"
    ) %>%
      group_by(.data$term) %>%
      summarise(
        bootstrap_n = sum(is.finite(.data$estimate)),
        bootstrap_mean = mean(.data$estimate, na.rm = TRUE),
        bootstrap_sd = sd(.data$estimate, na.rm = TRUE),
        bootstrap_ci_low = as.numeric(quantile(.data$estimate, 0.025, na.rm = TRUE)),
        bootstrap_ci_high = as.numeric(quantile(.data$estimate, 0.975, na.rm = TRUE)),
        bootstrap_sign_consistency = mean(sign(.data$estimate) == sign(.data$observed_estimate), na.rm = TRUE),
        .groups = "drop"
      )
  }

  observed_df %>%
    transmute(
      term = .data$term,
      observed_estimate = .data$estimate,
      observed_se = .data$se,
      observed_statistic = .data$statistic,
      observed_p = .data$p_value
    ) %>%
    left_join(boot_summary, by = "term")
}

run_feature_block <- function(df, platform, comparison, feature_col, feature_id, formula_str, terms, paired, cluster_var, strata_var = NULL, rename_vars = cluster_var, use_lm = FALSE) {
  feature_df <- df[df[[feature_col]] == feature_id, , drop = FALSE]
  if (nrow(feature_df) == 0) return(data.frame())

  fit_fun <- function(dat) {
    fit <- if (use_lm) {
      safe_lm(as.formula(formula_str), dat)
    } else {
      safe_lmer(as.formula(formula_str), dat)
    }
    extract_terms(fit, terms)
  }

  observed_df <- fit_fun(feature_df)
  boot_rows <- vector("list", n_boot)

  for (i in seq_len(n_boot)) {
    boot_df <- if (paired) {
      bootstrap_cluster_df(feature_df, cluster_var = cluster_var, rename_vars = rename_vars)
    } else {
      bootstrap_cluster_df(feature_df, cluster_var = cluster_var, strata_var = strata_var, rename_vars = rename_vars)
    }

    boot_fit <- fit_fun(boot_df)
    boot_fit$bootstrap_id <- i
    boot_rows[[i]] <- boot_fit
  }

  out <- summarize_resampling(observed_df, bind_rows(boot_rows))
  out$platform <- platform
  out$comparison <- comparison
  out$feature <- feature_id
  out
}

standardize_samplecode <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace(x, "-\\d+$", "")
  ifelse(
    stringr::str_detect(x, "_"),
    x,
    stringr::str_replace(x, "^(\\d+)(W\\d+)$", "\\1_\\2")
  )
}

prepare_nmr_long <- function() {
  dat <- readr::read_tsv(nmr_data_file, show_col_types = FALSE)
  meta <- readr::read_tsv(nmr_meta_file, show_col_types = FALSE)

  bins <- intersect(nmr_bins_keep, colnames(dat))
  if (length(bins) == 0) stop("None of the requested NMR bins were found.")

  meta2 <- meta %>%
    mutate(
      Sample = as.character(.data$Sample),
      Samplecode = stringr::str_replace(.data$Sample, "-\\d+$", "")
    ) %>%
    distinct(.data$Samplecode, .keep_all = TRUE) %>%
    select(all_of(c("Samplecode", "Subject")))

  dat <- dat %>%
    mutate(
      Samplecode = stringr::str_replace(as.character(.data$Samplecode), "-\\d+$", "")
    ) %>%
    left_join(meta2, by = "Samplecode") %>%
    mutate(
      Treatment = factor(.data$Treatment, levels = c("Control", "ACeD", "TCeD")),
      Gender = factor(.data$Gender),
      Subject = factor(.data$Subject),
      Samplecode = factor(.data$Samplecode),
      Age = suppressWarnings(as.numeric(.data$Age))
    )

  if (identical(nmr_transform_mode, "log2_shift")) {
    dat <- dat %>% mutate(across(all_of(bins), ~ nmr_log2_shift(.x)))
  }

  dat %>%
    select(all_of(c("Samplecode", "Subject", "Treatment", "Age", "Gender")), all_of(bins)) %>%
    pivot_longer(cols = all_of(bins), names_to = "feature", values_to = "y") %>%
    mutate(
      Treatment = droplevels(.data$Treatment),
      Gender = droplevels(.data$Gender),
      Subject = droplevels(.data$Subject),
      Samplecode = droplevels(.data$Samplecode)
    )
}

prepare_gcms_long <- function() {
  gc <- readr::read_tsv(gcms_data_file, show_col_types = FALSE)
  names(gc)[1] <- "Samplecode"
  gc$Samplecode <- standardize_samplecode(gc$Samplecode)

  meta <- readr::read_tsv(gcms_meta_file, show_col_types = FALSE) %>%
    mutate(
      Sample = as.character(.data$Sample),
      Samplecode = standardize_samplecode(.data$Sample),
      Subject = as.character(.data$Subject),
      Gender = as.character(.data$Gender),
      Age = suppressWarnings(as.numeric(.data$Age)),
      Treatment = as.character(.data$Treatment)
    )

  sub_info <- meta %>%
    group_by(.data$Subject) %>%
    summarise(
      Age = dplyr::first(.data$Age),
      Gender = dplyr::first(.data$Gender),
      .groups = "drop"
    )

  control_subjects <- meta %>%
    filter(.data$Treatment == "Control") %>%
    pull(.data$Subject) %>%
    unique()

  gc_meta <- gc %>%
    transmute(
      Samplecode = .data$Samplecode,
      Subject = stringr::str_extract(.data$Samplecode, "^\\d+"),
      Week = stringr::str_extract(.data$Samplecode, "W\\d+$")
    ) %>%
    mutate(
      Treatment = dplyr::case_when(
        .data$Subject %in% control_subjects ~ "Control",
        .data$Week == "W26" ~ "TCeD",
        .data$Week == "W0" ~ "ACeD",
        TRUE ~ NA_character_
      )
    ) %>%
    left_join(sub_info, by = "Subject") %>%
    mutate(
      Treatment = factor(.data$Treatment, levels = c("Control", "ACeD", "TCeD")),
      Gender = factor(.data$Gender),
      Subject = factor(.data$Subject),
      Age = suppressWarnings(as.numeric(.data$Age))
    )

  dat <- gc %>% left_join(gc_meta, by = "Samplecode")

  meta_cols <- c("Samplecode", "Treatment", "Age", "Gender", "Subject", "Week")
  features <- setdiff(names(dat), meta_cols)
  if (!is.null(gcms_features_keep)) {
    features <- intersect(gcms_features_keep, features)
  }
  if (length(features) == 0) stop("No GC-MS features were selected.")

  if (identical(gcms_transform_mode, "log2_shift")) {
    dat <- dat %>% mutate(across(all_of(features), ~ log2_shift(as.numeric(.x))))
  }

  dat %>%
    select(all_of(c("Samplecode", "Subject", "Treatment", "Age", "Gender")), all_of(features)) %>%
    pivot_longer(cols = all_of(features), names_to = "feature", values_to = "y") %>%
    mutate(
      Treatment = droplevels(.data$Treatment),
      Gender = droplevels(.data$Gender),
      Subject = droplevels(.data$Subject)
    )
}

run_nmr_robustness <- function(long_df) {
  feature_ids <- unique(long_df$feature)
  out <- list()
  kk <- 1L

  unpaired_df <- long_df %>%
    filter(.data$Treatment %in% c("Control", "ACeD", "TCeD"), is.finite(.data$y))

  paired_df <- long_df %>%
    filter(.data$Treatment %in% c("ACeD", "TCeD"), is.finite(.data$y)) %>%
    group_by(.data$Subject) %>%
    filter(n_distinct(.data$Treatment) == 2) %>%
    ungroup() %>%
    mutate(Treatment = factor(.data$Treatment, levels = c("ACeD", "TCeD")))

  for (feature_id in feature_ids) {
    out[[kk]] <- run_feature_block(
      df = unpaired_df,
      platform = "NMR",
      comparison = "vs_control",
      feature_col = "feature",
      feature_id = feature_id,
      formula_str = "y ~ Treatment + Age + Gender + (1 | Samplecode)",
      terms = c("TreatmentACeD", "TreatmentTCeD"),
      paired = FALSE,
      cluster_var = "Samplecode",
      strata_var = "Treatment",
      rename_vars = c("Samplecode")
    )
    kk <- kk + 1L

    out[[kk]] <- run_feature_block(
      df = paired_df,
      platform = "NMR",
      comparison = "paired_treated_vs_untreated",
      feature_col = "feature",
      feature_id = feature_id,
      formula_str = "y ~ Treatment + Age + Gender + (1 | Subject / Samplecode)",
      terms = c("TreatmentTCeD"),
      paired = TRUE,
      cluster_var = "Subject",
      rename_vars = c("Subject", "Samplecode")
    )
    kk <- kk + 1L
  }

  bind_rows(out)
}

run_gcms_robustness <- function(long_df) {
  feature_ids <- unique(long_df$feature)
  out <- list()
  kk <- 1L

  unpaired_df <- long_df %>%
    filter(.data$Treatment %in% c("Control", "ACeD", "TCeD"), is.finite(.data$y)) %>%
    # Subject-level case/control label: constant within Subject so the cluster
    # (Subject) is nested within the stratum. Stratifying by Treatment instead
    # splits a CeD subject's ACeD (W0) and TCeD (W26) rows across two strata,
    # collapsing every resampled subject to a single observation and breaking
    # the (1 | Subject) random effect.
    mutate(CaseControl = ifelse(.data$Treatment == "Control", "Control", "CeD"))

  paired_df <- long_df %>%
    filter(.data$Treatment %in% c("ACeD", "TCeD"), is.finite(.data$y)) %>%
    group_by(.data$Subject) %>%
    filter(n_distinct(.data$Treatment) == 2) %>%
    ungroup() %>%
    mutate(Treatment = factor(.data$Treatment, levels = c("ACeD", "TCeD")))

  for (feature_id in feature_ids) {
    out[[kk]] <- run_feature_block(
      df = unpaired_df,
      platform = "GCMS",
      comparison = "vs_control",
      feature_col = "feature",
      feature_id = feature_id,
      formula_str = "y ~ Treatment + Age + Gender + (1 | Subject)",
      terms = c("TreatmentACeD", "TreatmentTCeD"),
      paired = FALSE,
      cluster_var = "Subject",
      strata_var = "CaseControl",
      rename_vars = c("Subject")
    )
    kk <- kk + 1L

    out[[kk]] <- run_feature_block(
      df = paired_df,
      platform = "GCMS",
      comparison = "paired_treated_vs_untreated",
      feature_col = "feature",
      feature_id = feature_id,
      formula_str = "y ~ Treatment + Age + Gender + (1 | Subject)",
      terms = c("TreatmentTCeD"),
      paired = TRUE,
      cluster_var = "Subject",
      rename_vars = c("Subject")
    )
    kk <- kk + 1L
  }

  bind_rows(out)
}


# MAIN

set.seed(seed)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat(
  "Running metabolomics bootstrap robustness with",
  "n_boot =", n_boot,
  "out_dir =", out_dir,
  "\n"
)

nmr_results <- run_nmr_robustness(prepare_nmr_long())
gcms_results <- run_gcms_robustness(prepare_gcms_long())

all_results <- bind_rows(nmr_results, gcms_results) %>%
  group_by(.data$platform, .data$comparison, .data$term) %>%
  mutate(observed_q = p.adjust(.data$observed_p, method = "BH")) %>%
  ungroup()

all_results$robust_selected <- !is.na(all_results$observed_q) &
  all_results$observed_q <= q_cutoff &
  !is.na(all_results$bootstrap_ci_low) &
  !is.na(all_results$bootstrap_ci_high) &
  ((all_results$bootstrap_ci_low > 0 & all_results$bootstrap_ci_high > 0) |
     (all_results$bootstrap_ci_low < 0 & all_results$bootstrap_ci_high < 0))

write.csv(
  all_results,
  file = file.path(out_dir, "metabolomics_lmm_bootstrap.csv"),
  row.names = FALSE
)

cat("Metabolomics robustness output written to:", file.path(out_dir, "metabolomics_lmm_bootstrap.csv"), "\n")
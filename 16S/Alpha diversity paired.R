## Alpha diversity for paired untreated/treated CeD analysis

#Uses LMM to determine seven alpha diversity metrics while adjusting for covariants 

## Uses ps1_untreated_treated

library(phyloseq)
library(microbiome)
library(lme4)

setwd("C:/PhD/R code/cant_study")
## ── 0) Use ps1_untreated_treated and sanity-check it ────────────────────
ps1 <- readRDS("C:/PhD/Sequencing Results/Files/Cant_uf.rds")

# subjects who have BOTH an untreated and treated sample (Diagnosed.CD == TRUE)
meta <- as.data.frame(sample_data(ps1))

is_untreated_cd  <- with(meta, Diagnosed.CD == TRUE & GFD == FALSE & !is.na(Subject))
is_treated_cd <- with(meta, Diagnosed.CD == TRUE & GFD == TRUE  & !is.na(Subject))

untreated_subjects  <- unique(meta$Subject[is_untreated_cd])
treated_subjects <- unique(meta$Subject[is_treated_cd])
paired_subjects  <- intersect(untreated_subjects, treated_subjects)

paired_samples <- rownames(meta)[
  meta$Subject %in% paired_subjects &
    meta$Diagnosed.CD == TRUE &
    meta$GFD %in% c(FALSE, TRUE)
]
if (length(paired_samples) == 0) {
  stop("No paired samples found: no subjects have both untreated and treated Diagnosed.CD == TRUE samples. Check your metadata or source `Canterbury_object_prep.r`.")
}

ps1_untreated_treated <- prune_samples(paired_samples, ps1)
ps1_untreated_treated <- prune_taxa(taxa_sums(ps1_untreated_treated) > 0, ps1_untreated_treated)


sd_try <- try(sample_data(ps1_untreated_treated), silent = TRUE)
if (inherits(sd_try, "try-error")) {
  stop("sample_data(ps1_untreated_treated) is invalid. ",
       "Rebuild ps1_untreated_treated from ps1 with prune_samples(..., ps1) so metadata is preserved.")
}

sd <- as(sd_try, "data.frame")
if (nrow(sd) == 0 || ncol(sd) == 0) {
  stop("sample_data(ps1_untreated_treated) has zero rows or zero columns. ",
       "That will cause the 'Sample Data must have non-zero dimensions' error.")
}

# Required columns
needed_cols <- c("GFD", "Subject")
missing_cols <- setdiff(needed_cols, colnames(sd))
if (length(missing_cols) > 0) {
  stop("Missing required columns in sample_data(ps1_untreated_treated): ",
       paste(missing_cols, collapse = ", "))
}

# Define State (untreated vs treated) and ensure Subject is a factor
sd$State   <- factor(ifelse(sd$GFD, "Treated", "untreated"),
                     levels = c("untreated", "Treated"))  # untreated = reference
sd$Subject <- factor(sd$Subject)

GROUPING_VAR  <- "State"
ID_VAR        <- "Subject"
out_file_name <- "ps_alpha_diversity_lmm_results_paired.csv"


## ── 1) Alpha diversity + metadata merge ──────────────────────────────

# Alpha diversity metrics
alpha_div <- microbiome::alpha(
  ps1_untreated_treated,
  index = c("observed", "shannon", "fisher", "bulla",
            "simpson", "dominance", "low_abundance", "rare_abundance")
)
alpha_div$Sample.ID <- rownames(alpha_div)

# Metadata frame
sd$Sample.ID  <- rownames(sd)
sd$Read_depth <- as.numeric(sample_sums(ps1_untreated_treated))

# Merge alpha + meta
merged_df <- merge(alpha_div, sd, by = "Sample.ID")


## ── 2) Metric mapping ────────────────────────────────────────────────

metric_mapping <- list(
  observed       = "observed",
  shannon        = "diversity_shannon",
  fisher         = "diversity_fisher",
  bulla          = "evenness_bulla",
  simpson        = "evenness_simpson",
  dominance      = "dominance_simpson",
  rare_abundance = "rarity_low_abundance"
)


## ── 3) Per-metric stats (paired LMM, robust to dropped coefficients) ─

calculate_stats_for_metric <- function(metric, merged_df, mapping, pseudocount = 1) {
  actual_metric <- mapping[[metric]]
  if (is.null(actual_metric) || !(actual_metric %in% colnames(merged_df))) {
    stop("Metric ", metric, " not found in merged data.")
  }
  
  df <- merged_df
  df$log2_metric <- log2(df[[actual_metric]] + pseudocount)
  
  # Check required columns
  stopifnot(GROUPING_VAR %in% names(df),
            ID_VAR       %in% names(df),
            "Read_depth" %in% names(df))
  
  df[[GROUPING_VAR]] <- factor(df[[GROUPING_VAR]])
  df[[ID_VAR]]       <- factor(df[[ID_VAR]])
  # Rescale Read_depth to avoid warnings about very different predictor scales
  # (centered to mean 0, scaled to sd = 1). This stabilizes model fitting.
  df$Read_depth      <- as.numeric(df$Read_depth)
  if (all(is.finite(df$Read_depth))) {
    df$Read_depth <- as.numeric(scale(df$Read_depth))
  }
  
  # Need exactly two levels: untreated, Treated
  levs <- levels(df[[GROUPING_VAR]])
  if (length(levs) != 2) {
    warning("Metric ", metric, ": grouping var has levels ",
            paste(levs, collapse = ", "),
            " — skipping this metric.")
    return(NULL)
  }
  ref_lvl <- levs[1]
  alt_lvl <- levs[2]
  
  # LMM: log2(alpha) ~ State + Read_depth + (1|Subject)
  fm <- as.formula(
    paste0("log2_metric ~ ", GROUPING_VAR, " + Read_depth + (1|", ID_VAR, ")")
  )
  
  fit <- try(lmer(fm, data = df), silent = TRUE)
  if (inherits(fit, "try-error")) {
    warning("Metric ", metric, ": lmer failed, skipping.\n", fit)
    return(NULL)
  }
  
  coefs <- summary(fit)$coefficients
  coef_name <- paste0(GROUPING_VAR, alt_lvl)  # e.g. "StateTreated"
  
  if (coef_name %in% rownames(coefs)) {
    est  <- coefs[coef_name, "Estimate"]
    se   <- coefs[coef_name, "Std. Error"]

    if ("Pr(>|t|)" %in% colnames(coefs)) {
      pval <- coefs[coef_name, "Pr(>|t|)"]
    } else if ("t value" %in% colnames(coefs)) {
      tval <- coefs[coef_name, "t value"]
      pval <- 2 * pnorm(-abs(tval))
    } else {
      pval <- NA_real_
    }
  } else {
    # model ran but dropped the State effect (e.g. rank deficiency)
    warning("Metric ", metric, ": no coefficient row for ", coef_name,
            " (possibly rank-deficient); setting effect to NA.")
    est  <- NA_real_
    se   <- NA_real_
    pval <- NA_real_
  }
  sig <- ifelse(!is.na(pval) && pval < 0.05, "Significant", "ns")
  
  # overall stats
  mean_all <- mean(df$log2_metric, na.rm = TRUE)
  se_all   <- sd(df$log2_metric,   na.rm = TRUE) / sqrt(nrow(df))
  
  # group-wise stats
  get_stats <- function(level) {
    v <- df$log2_metric[df[[GROUPING_VAR]] == level]
    c(N  = length(v),
      M  = if (length(v)) mean(v, na.rm = TRUE) else NA_real_,
      SD = if (length(v)) sd(v,   na.rm = TRUE) else NA_real_)
  }
  s_ref <- get_stats(ref_lvl)
  s_alt <- get_stats(alt_lvl)
  
  data.frame(
    Metric           = metric,
    Ref_Level        = ref_lvl,
    Alt_Level        = alt_lvl,
    Mean_Log2_All    = mean_all,
    SE_Log2_All      = se_all,
    N_Ref            = s_ref["N"],
    Mean_Log2_Ref    = s_ref["M"],
    SD_Log2_Ref      = s_ref["SD"],
    N_Alt            = s_alt["N"],
    Mean_Log2_Alt    = s_alt["M"],
    SD_Log2_Alt      = s_alt["SD"],
    Coef_Alt_vs_Ref  = est,      # Treated vs untreated
    SE_Alt_vs_Ref    = se,
    P_Alt_vs_Ref     = pval,
    Significance     = sig,
    stringsAsFactors = FALSE
  )
}


## ── 4) Run across metrics and save ──────────────────────────────────

# Keep only metrics that actually exist in merged_df
valid_map <- metric_mapping[
  vapply(metric_mapping, `%in%`, logical(1), colnames(merged_df))
]
valid_names <- names(valid_map)

if (length(valid_names) == 0) {
  stop("No alpha-diversity metrics found in merged data!")
}

alpha_list <- lapply(valid_names, function(m) {
  calculate_stats_for_metric(m, merged_df, valid_map)
})

# Drop NULLs (metrics that failed / were skipped)
alpha_list <- alpha_list[!vapply(alpha_list, is.null, logical(1))]

alpha_results <- do.call(rbind, alpha_list)

write.csv(
  alpha_results,
  file      = out_file_name,  # change path if you want it somewhere else
  row.names = FALSE
)

alpha_results

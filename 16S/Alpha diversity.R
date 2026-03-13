#Alpha Diversity

#Non-paired

#Uses LMM to determine seven alpha diversity metrics while adjusting for covariants 

#Uses ps1_untreated/treated, which is unfiltered by prevalence

ps <- ps1_treated   
ps@sam_data$Diagnosed_Celiac <- factor(ps@sam_data$Diagnosed.CD)

## ── 0) Set Up ────────────────────────────

library(phyloseq)
library(microbiome)
library(meta)   # not strictly needed now, but kept in case you add meta-analysis later

# Primary biological factor in your sample_data(ps)
GROUPING_VAR <- "Diagnosed_Celiac"

# Output
out_file_name <- "ps_alpha_diversity_lmm_results_treated.csv"


## ── 1) Define Metric Mapping ────────────────────────────

metric_mapping <- list(
  observed       = "observed",
  shannon        = "diversity_shannon",
  fisher         = "diversity_fisher",
  bulla          = "evenness_bulla",
  simpson        = "evenness_simpson",
  dominance      = "dominance_simpson",
  rare_abundance = "rarity_low_abundance"
)


## ── 2) Per-metric stats function (with Age + Sex in LMM) ─────────────────

calculate_stats_for_metric <- function(metric, merged_df, mapping, pseudocount = 1) {
  actual_metric <- mapping[[metric]]
  if (is.null(actual_metric) || !(actual_metric %in% colnames(merged_df))) {
    stop("Metric ", metric, " not found in merged data.")
  }
  
  df <- merged_df
  df$log2_metric <- log2(df[[actual_metric]] + pseudocount)
  
  ## Make sure covariates are the right type
  df[[GROUPING_VAR]] <- factor(df[[GROUPING_VAR]])
  df$Sex             <- factor(df$Sex)
  df$Age             <- as.numeric(df$Age)
  
  # LMM / linear model:
  #   log2(metric) ~ Diagnosed_Celiac + Read_depth + Age + Sex
  lm_formula <- as.formula(
    paste("log2_metric ~", GROUPING_VAR, "+ Read_depth + Age + Sex")
  )
  lmcoefs <- summary(lm(lm_formula, data = df))$coefficients
  
  # overall mean / SE
  mean_all <- mean(df$log2_metric, na.rm = TRUE)
  se_all   <- sd(df$log2_metric,   na.rm = TRUE) / sqrt(nrow(df))
  
  # group‐wise mean, sd & N (still stratified by Diagnosed_Celiac)
  get_stats <- function(flag) {
    v <- df$log2_metric[df[[GROUPING_VAR]] == flag]
    c(N  = length(v),
      M  = if (length(v)) mean(v, na.rm = TRUE) else NA,
      SD = if (length(v)) sd(v,   na.rm = TRUE) else NA)
  }
  s_false <- get_stats("FALSE")
  s_true  <- get_stats("TRUE")
  
  ## Extract effect for Diagnosed_Celiac (adjusted for depth, Age, Sex)
  coef_name <- paste0(GROUPING_VAR, "TRUE")
  if (coef_name %in% rownames(lmcoefs)) {
    est  <- lmcoefs[coef_name, "Estimate"]
    se   <- lmcoefs[coef_name, "Std. Error"]
    pval <- lmcoefs[coef_name, "Pr(>|t|)"]
  } else {
    est <- se <- pval <- NA
  }
  sig <- ifelse(!is.na(pval) && pval < .05, "Significant", "ns")
  
  res <- data.frame(
    Metric             = metric,
    Mean_Log2_All      = mean_all,
    SE_Log2_All        = se_all,
    N_FALSE            = s_false["N"],
    Mean_Log2_FALSE    = s_false["M"],
    SD_Log2_FALSE      = s_false["SD"],
    N_TRUE             = s_true["N"],
    Mean_Log2_TRUE     = s_true["M"],
    SD_Log2_TRUE       = s_true["SD"],
    Coef_Celiac        = est,
    SE_Celiac          = se,
    P_Celiac           = pval,
    Significance       = sig,
    stringsAsFactors   = FALSE
  )
  names(res)[names(res) == "Coef_Celiac"] <- paste0("Coef_", GROUPING_VAR)
  names(res)[names(res) == "SE_Celiac"]   <- paste0("SE_", GROUPING_VAR)
  names(res)[names(res) == "P_Celiac"]    <- paste0("P_", GROUPING_VAR)
  
  res
}


## ── 3) Wrapper to run on your single phyloseq object ps ──────────────

run_alpha_analysis <- function(physeq_obj, pseudocount = 1) {
  # Alpha diversity
  alpha_div <- microbiome::alpha(
    physeq_obj,
    index = c("observed", "shannon", "fisher", "bulla",
              "simpson", "dominance", "low_abundance", "rare_abundance")
  )
  alpha_div$Sample.ID <- row.names(alpha_div)
  
  # Metadata from phyloseq
  meta_df <- data.frame(sample_data(physeq_obj))
  meta_df$Sample.ID  <- row.names(meta_df)
  meta_df$Read_depth <- sample_sums(physeq_obj)
  
  # Merge alpha + meta (Age, Sex, Diagnosed_Celiac should be in meta_df)
  merged_df <- merge(alpha_div, meta_df, by = "Sample.ID")
  
  # Keep only metrics that exist in merged_df
  valid_map <- metric_mapping[
    vapply(metric_mapping, `%in%`, logical(1), colnames(merged_df))
  ]
  valid_names <- names(valid_map)
  if (length(valid_names) == 0) stop("No metrics found in merged data!")
  
  do.call(
    rbind,
    lapply(valid_names, function(m) {
      calculate_stats_for_metric(m, merged_df, valid_map, pseudocount)
    })
  )
}


## ── 4) Run on ps and save ───────────────────────────────────

alpha_results <- run_alpha_analysis(ps)

write.csv(alpha_results,
          file = out_file_name,
          row.names = FALSE)

alpha_results

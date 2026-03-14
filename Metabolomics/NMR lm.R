#NMR Linear model analysis

#Runs on replicate level PQN normalised NMR bins that represent non-overalpping bins used to idenitfy metabolites see Supplementary Figure 2
#ACeD = untreated samples; TCeD = treated samples

## 0) Packages
pkgs <- c("tidyverse", "lmerTest", "emmeans", "broom.mixed", "readr")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
lapply(pkgs, library, character.only = TRUE)


## 1) Inputs

data_file <- "C:\\PhD\\R code\\nmr\\data_CoeliacStoolCHCH_pqn_selected_bins.tsv"  # replicate-level PQN data
meta_file <- "C:\\PhD\\NMR\\metablo_meta.txt"                          # replicate-level metadata with Subject

# Your final QC'd, single-metabolite bins (edit if needed)
bins_keep <- c(
  "B1_0398","B1_5414","B2_6350","B2_3582","B2_7966",
  "B1_7458","B2_2856","B3_1092","B1_4742","B0_9539","B1_0209",
  "B0_9982","B1_9199","B2_0511","B3_9598"
)


# Transformation options:
transform_mode <- "log2_shift"  # options: "none" or "log2_shift"
if (!transform_mode %in% c("none", "log2_shift")) {
  stop("transform_mode must be one of: 'none', 'log2_shift'")
}


## 2) Helpers

log2_shift <- function(x) {
  xmin <- min(x, na.rm = TRUE)
  shift <- ifelse(xmin <= 0, -xmin + 1, 0)
  log2(x + shift)
}

safe_lmer <- function(formula, data) {
  tryCatch(
    lmerTest::lmer(formula, data = data, REML = FALSE),
    error = function(e) {
      cat("Warning: lmer fit failed (", conditionMessage(e), "); returning NULL\n", sep = "")
      NULL
    }
  )
}

empty_res <- function() {
  tibble::tibble(
    Bin = NA_character_, contrast = NA_character_, estimate = NA_real_,
    std.error = NA_real_, df = NA_real_, statistic = NA_real_, p.value = NA_real_
  )[0, ]
}


## 3) Load data

dat <- readr::read_tsv(data_file, show_col_types = FALSE)
# Keep only numeric/logical bin columns
bin_cols_all <- intersect(bins_keep, colnames(dat))
if (length(bin_cols_all) == 0) stop("No bins_keep found in data file.")

bins_keep <- intersect(bins_keep, bin_cols_all)
if (length(bins_keep) == 0) stop("None of your bins_keep are present in the data file.")

meta <- readr::read_tsv(meta_file, show_col_types = FALSE)

# meta has replicate-level Sample like "1002W26-1"; extract base Samplecode and Subject
meta2 <- meta %>%
  mutate(
    Sample = as.character(Sample),
    Samplecode = stringr::str_replace(Sample, "-\\d+$", "")
  ) %>%
  distinct(Samplecode, .keep_all = TRUE) %>%
  dplyr::select(Samplecode, Subject)

# Merge Subject into replicate data
dat <- dat %>%
  mutate(
    Samplecode = as.character(Samplecode),
    # Strip replicate suffix (e.g., "1002W0-1" -> "1002W0")
    Samplecode = stringr::str_replace(Samplecode, "-\\d+$", "")
  ) %>%
  left_join(meta2, by = "Samplecode") %>%
  mutate(
    Treatment = factor(Treatment),
    Gender    = factor(Gender),
    Subject   = factor(Subject),
    Age       = suppressWarnings(as.numeric(Age))
  )

if (any(is.na(dat$Subject))) {
  warning("Some rows have missing Subject after merge. Check Samplecode naming consistency.")
}

# Optional transform
if (transform_mode == "log2_shift") {
  dat <- dat %>%
    mutate(across(all_of(bins_keep), ~ log2_shift(.x)))
}

## 4) Long format (one row per sample-bin)
long <- dat %>%
  select(Samplecode, Subject, Treatment, Age, Gender, all_of(bins_keep)) %>%
  pivot_longer(cols = all_of(bins_keep), names_to = "Bin", values_to = "y") %>%
  mutate(
    Bin = as.character(Bin),
    Treatment = droplevels(Treatment),
    Gender = droplevels(Gender),
    Subject = droplevels(Subject)
  )

## 5) Per-bin models with explicit (named) emmeans contrasts

# --- Active/Treated vs Control (per bin)
fit_vs_control <- function(df_bin, bin_name) {

  df_bin <- df_bin %>% filter(is.finite(y))

  if (n_distinct(df_bin$Treatment) < 2) return(empty_res())
  if (!("Control" %in% unique(df_bin$Treatment))) return(empty_res())

  df_bin <- df_bin %>%
    mutate(
      Treatment = droplevels(factor(Treatment)),
      Gender = droplevels(factor(Gender)),
      Subject = droplevels(factor(Subject))
    )

  df_bin$Treatment <- relevel(df_bin$Treatment, ref = "Control")

  # Fit mixed model with random intercept per Samplecode (accounts for replicate-level variation)
  # Note: No Subject random effect since each subject appears in only one treatment group
  m <- tryCatch(
    lmerTest::lmer(y ~ Treatment + Age + Gender + (1 | Samplecode), data = df_bin, REML = FALSE),
    error = function(e) {
      cat("Error in lmer for bin", bin_name, ":", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(m)) return(empty_res())

  # Extract fixed effects and vcov
  fe <- lme4::fixef(m)
  vc <- stats::vcov(m)
  
  # Find Treatment coefficients (non-Control)
  treat_coefs <- grep("^Treatment", names(fe), value = TRUE)
  
  if (length(treat_coefs) == 0) return(empty_res())
  
  # Build results (normal reference; df not used)
  t_stats <- as.numeric(fe[treat_coefs]) / sqrt(diag(vc)[treat_coefs])
  tibble::tibble(
    Bin = bin_name,
    contrast = stringr::str_remove(treat_coefs, "^Treatment"),
    estimate = as.numeric(fe[treat_coefs]),
    std.error = sqrt(diag(vc)[treat_coefs]),
    df = NA_real_,
    statistic = t_stats,
    p.value = 2 * stats::pnorm(-abs(t_stats))
  )
}

# --- Paired TCeD vs ACeD (per bin)
fit_paired <- function(df_bin, bin_name) {

  df_bin <- df_bin %>%
    filter(Treatment %in% c("ACeD", "TCeD")) %>%
    filter(is.finite(y)) %>%
    mutate(
      Treatment = factor(Treatment, levels = c("ACeD", "TCeD")),
      Gender = droplevels(factor(Gender)),
      Subject = droplevels(factor(Subject))
    )

  # Keep only subjects with both treatments
  df_bin <- df_bin %>%
    group_by(Subject) %>%
    filter(n_distinct(Treatment) == 2) %>%
    ungroup()

  if (nrow(df_bin) < 4) return(empty_res())

  # Use nested random effects: Subject + Samplecode within Subject (for replicate-level data)
  m <- safe_lmer(y ~ Treatment + Age + Gender + (1 | Subject / Samplecode), df_bin)
  if (is.null(m)) return(empty_res())

  # Extract fixed effects and vcov
  fe <- lme4::fixef(m)
  vc <- stats::vcov(m)
  
  # For paired design, TCeD coefficient vs ACeD (reference)
  treat_coef <- grep("^TreatmentTCeD", names(fe), value = TRUE)
  
  if (length(treat_coef) == 0) return(empty_res())
  
  t_stat <- as.numeric(fe[treat_coef]) / sqrt(diag(vc)[treat_coef])
  
  tibble::tibble(
    Bin = bin_name,
    contrast = "TCeD - ACeD",
    estimate = as.numeric(fe[treat_coef]),
    std.error = sqrt(diag(vc)[treat_coef]),
    df = NA_real_,
    statistic = t_stat,
    p.value = 2 * stats::pnorm(-abs(t_stat))
  )
}


## 6) Run models across bins + FDR (combined BH across all contrasts)


# Run once for all bins (returns 2 contrasts per bin when present: ACeD and TCeD vs Control)
res_vs_ctrl_all <- split(long, long$Bin) %>%
  purrr::imap_dfr(~ fit_vs_control(.x, .y))

# Paired TCeD - ACeD
res_paired <- split(long, long$Bin) %>%
  purrr::imap_dfr(~ fit_paired(.x, .y))

# Single BH across ALL contrasts/bins (less conservative)
all_results <- dplyr::bind_rows(
  res_vs_ctrl_all %>% dplyr::mutate(contrast = dplyr::case_when(
    contrast == "ACeD" ~ "ACeD - Control",
    contrast == "TCeD" ~ "TCeD - Control",
    TRUE ~ contrast
  ), family = "vs_ctrl"),
  res_paired %>% dplyr::mutate(contrast = "TCeD - ACeD", family = "paired")
)

all_results <- all_results %>%
  dplyr::mutate(
    q.value = p.adjust(p.value, method = "BH"),
    fold_change = 2^estimate
  )

# Split back
res_vs_ctrl_all <- all_results %>% dplyr::filter(family == "vs_ctrl") %>% dplyr::select(-family)
res_paired <- all_results %>% dplyr::filter(family == "paired") %>% dplyr::select(-family)

# Split tables
res_ACeD_vs_ctrl <- res_vs_ctrl_all %>% dplyr::filter(contrast == "ACeD - Control")
res_TCeD_vs_ctrl <- res_vs_ctrl_all %>% dplyr::filter(contrast == "TCeD - Control")

# Combined file for reference
res_vs_ctrl_combined <- dplyr::bind_rows(
  res_ACeD_vs_ctrl,
  res_TCeD_vs_ctrl
)

# Write files
readr::write_csv(res_ACeD_vs_ctrl,   "C:\\PhD\\NMR\\LMM_bins_ACeD_vs_Control.csv")
readr::write_csv(res_TCeD_vs_ctrl,   "C:\\PhD\\NMR\\LMM_bins_TCeD_vs_Control.csv")
readr::write_csv(res_vs_ctrl_combined,"C:\\PhD\\NMR\\LMM_bins_contrasts_vs_control.csv")
readr::write_csv(res_paired, "C:\\PhD\\NMR\\LMM_bins_paired_TCeD_vs_ACeD.csv")


## 7) Console summaries

cat("\nSaved:\n  - LMM_bins_contrasts_vs_control.csv\n  - LMM_bins_paired_TCeD_vs_ACeD.csv\n\n")

cat("Top (vs Control):\n")
print(res_vs_ctrl_combined %>% arrange(q.value) %>% head(10))

cat("\nTop (paired TCeD - ACeD):\n")
print(res_paired %>% arrange(q.value) %>% head(10))


## 8) SVG boxplots per bin (Treatment + sample points)


# Output folder for plots
plot_dir <- "C:\\PhD\\NMR\\bin_boxplots_svg"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Ensure Treatment order (edit if you want a different order)
long <- long %>%
  mutate(Treatment = factor(Treatment, levels = c("Control", "ACeD", "TCeD")))

# Function to save one plot per bin
save_bin_boxplot_svg <- function(bin_name) {

  dfp <- long %>%
    dplyr::filter(Bin == bin_name, is.finite(y)) %>%
    dplyr::filter(!is.na(Treatment))

  # Skip if no data
  if (nrow(dfp) == 0) return(invisible(NULL))

  # Optional: show n per group in x-axis labels
  n_tbl <- dfp %>%
    dplyr::count(Treatment, name = "n") %>%
    dplyr::mutate(lbl = paste0(as.character(Treatment), "\n(n=", n, ")"))

  dfp <- dfp %>%
    dplyr::left_join(n_tbl %>% dplyr::select(Treatment, lbl), by = "Treatment")

  p <- ggplot(dfp, aes(x = lbl, y = y)) +
    geom_boxplot(width = 0.55, outlier.shape = NA) +
    geom_jitter(width = 0.12, height = 0, alpha = 0.85, size = 2) +
    labs(
      title = paste0(bin_name, " (PQN bin intensity)"),
      x = "Treatment",
      y = "Intensity"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(size = 10)
    )

  out_file <- file.path(plot_dir, paste0(bin_name, "_boxplot.svg"))

  # Use base R SVG device (no extra packages required)
  ggsave(filename = out_file, plot = p,
         device = grDevices::svg, width = 5.5, height = 4.2, units = "in")
}

# Generate plots for all bins you’re modelling
purrr::walk(unique(long$Bin), save_bin_boxplot_svg)

cat("\nSVG boxplots written to: ", plot_dir, "\n", sep = "")


## 9) PCA plot (SVG) + 95% ellipses + PERMANOVA


# Ensure vegan is available for PERMANOVA
if (!requireNamespace("vegan", quietly = TRUE)) {
  install.packages("vegan", repos = "https://cloud.r-project.org")
}
library(vegan)

# Output paths
pca_svg   <- "C:\\PhD\\NMR\\PCA_by_Treatment_95ellipse.svg"
perm_csv  <- "C:\\PhD\\NMR\\PERMANOVA_Treatment_Age_Gender.csv"

# Choose scaling for PCA: "none", "pareto", or "auto"
pca_scale <- "pareto"

# Optional override if you want PCA log2 even when transform_mode == "none"
pca_force_log2 <- FALSE

# Placeholder label (updated after PERMANOVA below)
perm_label <- "PERMANOVA pending"

# Build PCA matrix (rows = samples, cols = bins)
# Average replicates to sample-level first to match LMM analysis
pca_df <- dat %>%
  group_by(Samplecode, Treatment, Age, Gender, Subject) %>%
  summarise(across(all_of(bins_keep), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  dplyr::filter(!is.na(Treatment))

# Keep a metadata table in the same row order
perm_meta <- pca_df %>%
  dplyr::select(Samplecode, Treatment, Age, Gender, Subject) %>%
  dplyr::mutate(
    Treatment = factor(Treatment, levels = c("Control", "ACeD", "TCeD")),
    Gender = factor(Gender),
    Subject = factor(Subject),
    Age = suppressWarnings(as.numeric(Age))
  )

X <- as.matrix(pca_df[, bins_keep, drop = FALSE])
X <- apply(X, 2, function(z) as.numeric(z))

# Handle NA by feature median
for (j in seq_len(ncol(X))) {
  med <- median(X[, j], na.rm = TRUE)
  X[is.na(X[, j]), j] <- med
}

# Optional log2 transform (if desired)
if (pca_force_log2 && transform_mode == "none") {
  for (j in seq_len(ncol(X))) {
    xmin <- min(X[, j], na.rm = TRUE)
    X[, j] <- X[, j] - xmin
    pos <- X[X[, j] > 0, j]
    pc  <- if (length(pos)) min(pos, na.rm = TRUE) / 2 else 1e-6
    X[, j] <- log2(X[, j] + pc)
  }
}

# Center
Xc <- scale(X, center = TRUE, scale = FALSE)

# Scale
if (pca_scale == "auto") {
  Xc <- scale(Xc, center = FALSE, scale = apply(Xc, 2, sd))
} else if (pca_scale == "pareto") {
  Xc <- scale(Xc, center = FALSE, scale = sqrt(apply(Xc, 2, sd)))
}

# PCA
p <- prcomp(Xc, center = FALSE, scale. = FALSE)

# Variance explained
ve <- (p$sdev^2) / sum(p$sdev^2)
pc1_lab <- paste0("PC1 (", round(ve[1] * 100, 1), "%)")
pc2_lab <- paste0("PC2 (", round(ve[2] * 100, 1), "%)")

# --- scores table (add Subject and a time variable for ordering)
scores <- as.data.frame(p$x[, 1:2, drop = FALSE])
scores$Samplecode <- perm_meta$Samplecode
scores$Treatment  <- perm_meta$Treatment
scores$Subject    <- perm_meta$Subject

scores <- scores %>%
  mutate(
    # robust: works for "1002W0", "1002W26", "1002_W0", "1002_W26", etc.
    Time = case_when(
      stringr::str_detect(Samplecode, "W0")  ~ 0,
      stringr::str_detect(Samplecode, "W26") ~ 26,
      TRUE ~ NA_real_
    )
  )

# 95% ellipse data: only groups with >=3 points
ell_df <- scores %>%
  dplyr::group_by(Treatment) %>%
  dplyr::filter(dplyr::n() >= 3) %>%
  dplyr::ungroup()

# --- subject linking lines: only subjects with >1 point
link_df <- scores %>%
  dplyr::filter(!is.na(Subject)) %>%
  dplyr::group_by(Subject) %>%
  dplyr::filter(dplyr::n() >= 2) %>%
  dplyr::arrange(Time, .by_group = TRUE) %>%
  dplyr::ungroup()

## -------------------------
## Plot with 95% ellipses + annotation
## -------------------------
p_pca <- ggplot(scores, aes(x = PC1, y = PC2)) +
  geom_path(
    data = link_df,
    aes(group = Subject),
    color = "grey50",
    alpha = 0.6,
    linewidth = 0.6
  ) +
  geom_point(aes(color = Treatment), size = 3, alpha = 0.9) +
  stat_ellipse(
    data = ell_df,
    aes(color = Treatment),
    level = 0.95,
    type = "t",
    linewidth = 0.8
  ) +
  scale_color_manual(
    values = c(ACeD = "red", TCeD = "blue", Control = "green"),
    breaks = c("Control", "ACeD", "TCeD"),
    drop = FALSE
  ) +
  labs(
    title = "PCA of NMR bins",
    x = pc1_lab,
    y = pc2_lab
  ) +
  annotate(
    "text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.2,
    label = perm_label, size = 3.3
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(filename = pca_svg, plot = p_pca, device = grDevices::svg, width = 7, height = 5.5, units = "in")
cat("\nPCA SVG saved to: ", pca_svg, "\n", sep = "")
## If subjects repeat, restrict permutations within Subject (strata)
## -------------------------
dist_mat <- dist(Xc, method = "euclidean")

set.seed(123)

use_strata <- !all(is.na(perm_meta$Subject)) && any(duplicated(perm_meta$Subject))
if (use_strata) {
  ad <- vegan::adonis2(
    dist_mat ~ Treatment + Age + Gender,
    data = perm_meta,
    permutations = 9999,
    by = "margin",
    strata = perm_meta$Subject
  )
} else {
  ad <- vegan::adonis2(
    dist_mat ~ Treatment + Age + Gender,
    data = perm_meta,
    permutations = 9999,
    by = "margin"
  )
}

ad_df <- as.data.frame(ad)
ad_df$term <- rownames(ad_df)
readr::write_csv(ad_df, perm_csv)

# Extract Treatment R2 and p-value for annotation
tr_row <- ad_df %>% dplyr::filter(term == "Treatment")
t_r2 <- if (nrow(tr_row) == 1) tr_row$R2 else NA_real_
t_p  <- if (nrow(tr_row) == 1) tr_row$`Pr(>F)` else NA_real_

perm_label <- if (is.finite(t_r2) && is.finite(t_p)) {
  paste0("PERMANOVA (Treatment | Age, Gender): R²=", sprintf("%.3f", t_r2),
         ", p=", format(t_p, digits = 3))
} else {
  "PERMANOVA: Treatment row not found"
}

cat("\nPERMANOVA table saved to: ", perm_csv, "\n", sep = "")
print(ad)

## -------------------------
## Plot with 95% ellipses + annotation
## -------------------------
p_pca <- ggplot(scores, aes(x = PC1, y = PC2)) +
  geom_path(
    data = link_df,
    aes(group = Subject),
    color = "grey50",
    alpha = 0.6,
    linewidth = 0.6
  ) +
  geom_point(aes(color = Treatment), size = 3, alpha = 0.9) +
  stat_ellipse(
    data = ell_df,
    aes(color = Treatment),
    level = 0.95,
    type = "t",
    linewidth = 0.8
  ) +
  scale_color_manual(
    values = c(ACeD = "red", TCeD = "blue", Control = "green"),
    breaks = c("Control", "ACeD", "TCeD"),
    drop = FALSE
  ) +
  labs(
    title = "PCA of NMR bins",
    x = pc1_lab,
    y = pc2_lab
  ) +
  annotate(
    "text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.2,
    label = perm_label, size = 3.3
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(filename = pca_svg, plot = p_pca, device = grDevices::svg, width = 7, height = 5.5, units = "in")
cat("\nPCA SVG saved to: ", pca_svg, "\n", sep = "")

## GC-MS LMM + Boxplots + PCA + PERMANOVA (adapted from NMR script)

## ACeD = untreated samples; TCeD = treated

## 0) Packages

pkgs <- c("tidyverse", "lmerTest", "readr", "vegan")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
lapply(pkgs, library, character.only = TRUE)


## 1) Inputs

data_file <- "C:\\PhD\\Canterbury Dataset\\GC-MS\\GC_MS_metab.txt"     # GC-MS matrix (samples in rows, metabolites in columns)
meta_file <- "C:\\PhD\\NMR\\metablo_meta.txt"     # replicate-level metadata (used to derive Subject/Age/Gender/controls)

out_dir  <- "C:\\PhD\\GCMS\\outputs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Feature selection:

features_keep <- NULL

# Transformation options:
# - "none": use values as-is
# - "log2_shift": log2 after shifting per feature
transform_mode <- "log2_shift"

# PCA scaling: "none", "pareto", or "auto"
pca_scale <- "pareto"


## 2) Helpers

log2_shift <- function(x) {
  xmin <- min(x, na.rm = TRUE)
  x_shift <- x - xmin
  pos <- x_shift[is.finite(x_shift) & (x_shift > 0)]
  pc <- if (length(pos)) min(pos, na.rm = TRUE) / 2 else 1e-6
  log2(x_shift + pc)
}

safe_lmer <- function(formula, data) {
  tryCatch(
    lmerTest::lmer(formula, data = data, REML = FALSE),
    error = function(e) NULL
  )
}

empty_res <- function() {
  tibble(
    Metabolite = character(),
    contrast   = character(),
    estimate   = double(),
    std.error  = double(),
    df         = double(),
    statistic  = double(),
    p.value    = double()
  )
}

# Standardize sample codes:
#  - remove "-rep" suffix like "-1"
#  - convert "1002W26" -> "1002_W26"
#  - keep "1002_W26" as-is
standardize_samplecode <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace(x, "-\\d+$", "")
  x <- ifelse(stringr::str_detect(x, "_"),
              x,
              stringr::str_replace(x, "^(\\d+)(W\\d+)$", "\\1_\\2"))
  x
}

# Safe file names for metabolites (spaces, punctuation)
safe_filename <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "[^A-Za-z0-9]+", "_")
  x <- stringr::str_replace_all(x, "_+", "_")
  x <- stringr::str_replace(x, "^_|_$", "")
  substr(x, 1, 120)
}


## 3) Read GC-MS data + metadata

gc <- readr::read_tsv(data_file, show_col_types = FALSE)

# Rename first column to Samplecode
names(gc)[1] <- "Samplecode"
gc <- gc %>%
  mutate(Samplecode = standardize_samplecode(Samplecode))

# Read metadata and standardize meta samplecodes
meta <- readr::read_tsv(meta_file, show_col_types = FALSE) %>%
  mutate(
    Sample     = as.character(Sample),
    Samplecode = standardize_samplecode(Sample),
    Subject    = as.character(Subject),
    Gender     = as.character(Gender),
    Age        = suppressWarnings(as.numeric(Age)),
    Treatment  = as.character(Treatment)
  )

# Subject-level Age/Gender 
sub_info <- meta %>%
  group_by(Subject) %>%
  summarise(
    Age    = dplyr::first(Age),
    Gender = dplyr::first(Gender),
    .groups = "drop"
  )

control_subjects <- meta %>%
  filter(Treatment == "Control") %>%
  pull(Subject) %>%
  unique()

# Build GC-MS sample-level metadata 
gc_meta <- gc %>%
  transmute(
    Samplecode = Samplecode,
    Subject = stringr::str_extract(Samplecode, "^\\d+"),
    Week    = stringr::str_extract(Samplecode, "W\\d+$")
  ) %>%
  mutate(
    Treatment = dplyr::case_when(
      Subject %in% control_subjects ~ "Control",
      Week == "W26" ~ "TCeD",
      Week == "W0"  ~ "ACeD",
      TRUE ~ NA_character_
    )
  ) %>%
  left_join(sub_info, by = "Subject") %>%
  mutate(
    Treatment = factor(Treatment, levels = c("Control","ACeD","TCeD")),
    Gender    = factor(Gender),
    Subject   = factor(Subject),
    Age       = suppressWarnings(as.numeric(Age))
  )

# Merge metadata into GC-MS data
dat <- gc %>%
  left_join(gc_meta, by = "Samplecode")

if (any(is.na(dat$Treatment))) {
  warning("Some rows have missing Treatment after assignment. Check sample naming or control_subjects definition.")
}

# Identify feature columns (all except metadata columns)
meta_cols <- c("Samplecode","Treatment","Age","Gender","Subject","Week")
feat_cols_all <- setdiff(names(dat), meta_cols)

if (!is.null(features_keep)) {
  features_keep <- intersect(features_keep, feat_cols_all)
  if (length(features_keep) == 0) stop("None of features_keep are present in the GC-MS data file.")
} else {
  features_keep <- feat_cols_all
}

# Optional transform
if (transform_mode == "log2_shift") {
  dat <- dat %>%
    mutate(across(all_of(features_keep), ~ log2_shift(as.numeric(.x))))
}


## 4) Long format (one row per sample-metabolite)

long <- dat %>%
  select(Samplecode, Subject, Treatment, Age, Gender, all_of(features_keep)) %>%
  pivot_longer(cols = all_of(features_keep), names_to = "Metabolite", values_to = "y") %>%
  mutate(
    Metabolite = as.character(Metabolite),
    Treatment  = droplevels(Treatment),
    Gender     = droplevels(Gender),
    Subject    = droplevels(Subject)
  )


## 5) Per-metabolite models

# --- ACeD/TCeD vs Control
fit_vs_control <- function(df_feat, feat_name) {

  df_feat <- df_feat %>% filter(is.finite(y))

  if (n_distinct(df_feat$Treatment) < 2) return(empty_res())
  if (!("Control" %in% unique(df_feat$Treatment))) return(empty_res())

  df_feat <- df_feat %>%
    mutate(
      Treatment = droplevels(factor(Treatment)),
      Gender    = droplevels(factor(Gender)),
      Subject   = droplevels(factor(Subject))
    )

  df_feat$Treatment <- relevel(df_feat$Treatment, ref = "Control")

  m <- safe_lmer(y ~ Treatment + Age + Gender + (1 | Subject), df_feat)
  if (is.null(m)) return(empty_res())

  fe <- lme4::fixef(m)
  vc <- stats::vcov(m)

  treat_coefs <- grep("^Treatment", names(fe), value = TRUE)
  if (length(treat_coefs) == 0) return(empty_res())

  tibble::tibble(
    Metabolite = feat_name,
    contrast   = stringr::str_remove(treat_coefs, "^Treatment"),  # "ACeD" or "TCeD"
    estimate   = as.numeric(fe[treat_coefs]),
    std.error  = sqrt(diag(vc)[treat_coefs]),
    df         = NA_real_,
    statistic  = as.numeric(fe[treat_coefs]) / sqrt(diag(vc)[treat_coefs]),
    p.value    = 2 * stats::pnorm(-abs(statistic))
  )
}

# --- Paired TCeD vs ACeD within subjects (requires both timepoints)
fit_paired <- function(df_feat, feat_name) {

  df_feat <- df_feat %>%
    filter(Treatment %in% c("ACeD", "TCeD")) %>%
    filter(is.finite(y)) %>%
    mutate(
      Treatment = factor(Treatment, levels = c("ACeD", "TCeD")),
      Gender    = droplevels(factor(Gender)),
      Subject   = droplevels(factor(Subject))
    )

  df_feat <- df_feat %>%
    group_by(Subject) %>%
    filter(n_distinct(Treatment) == 2) %>%
    ungroup()

  if (nrow(df_feat) < 4) return(empty_res())

  m <- safe_lmer(y ~ Treatment + Age + Gender + (1 | Subject), df_feat)
  if (is.null(m)) return(empty_res())

  fe <- lme4::fixef(m)
  vc <- stats::vcov(m)

  treat_coef <- grep("^TreatmentTCeD", names(fe), value = TRUE)
  if (length(treat_coef) == 0) return(empty_res())

  tibble::tibble(
    Metabolite = feat_name,
    contrast   = "TCeD - ACeD",
    estimate   = as.numeric(fe[treat_coef]),
    std.error  = sqrt(diag(vc)[treat_coef]),
    df         = NA_real_,
    statistic  = as.numeric(fe[treat_coef]) / sqrt(diag(vc)[treat_coef]),
    p.value    = 2 * stats::pnorm(-abs(statistic))
  )
}

## 6) Run models across metabolites + FDR (split by comparison)


res_vs_ctrl_all <- split(long, long$Metabolite) %>%
  purrr::imap_dfr(~ fit_vs_control(.x, .y))

res_ACeD_vs_ctrl <- res_vs_ctrl_all %>% filter(contrast == "ACeD") %>%
  mutate(q.value = p.adjust(p.value, method = "BH"))

res_TCeD_vs_ctrl <- res_vs_ctrl_all %>% filter(contrast == "TCeD") %>%
  mutate(q.value = p.adjust(p.value, method = "BH"))

res_vs_ctrl_combined <- bind_rows(
  res_ACeD_vs_ctrl %>% mutate(contrast = "ACeD - Control"),
  res_TCeD_vs_ctrl %>% mutate(contrast = "TCeD - Control")
)

write_csv(res_ACeD_vs_ctrl,    file.path(out_dir, "GCMS_LMM_ACeD_vs_Control.csv"))
write_csv(res_TCeD_vs_ctrl,    file.path(out_dir, "GCMS_LMM_TCeD_vs_Control.csv"))
write_csv(res_vs_ctrl_combined,file.path(out_dir, "GCMS_LMM_contrasts_vs_control.csv"))

res_paired <- split(long, long$Metabolite) %>%
  purrr::imap_dfr(~ fit_paired(.x, .y)) %>%
  mutate(q.value = p.adjust(p.value, method = "BH"))

write_csv(res_paired, file.path(out_dir, "GCMS_LMM_paired_TCeD_vs_ACeD.csv"))

cat("\nSaved LMM results to: ", out_dir, "\n", sep = "")
cat("Top (ACeD vs Control):\n"); print(res_ACeD_vs_ctrl %>% arrange(q.value) %>% head(10))
cat("\nTop (TCeD vs Control):\n"); print(res_TCeD_vs_ctrl %>% arrange(q.value) %>% head(10))
cat("\nTop (paired TCeD - ACeD):\n"); print(res_paired %>% arrange(q.value) %>% head(10))


## 7) SVG boxplots per metabolite (Treatment + sample points)

plot_dir <- file.path(out_dir, "metabolite_boxplots_svg")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

long <- long %>%
  mutate(Treatment = factor(Treatment, levels = c("Control", "ACeD", "TCeD")))

save_met_boxplot_svg <- function(met_name) {

  dfp <- long %>%
    filter(Metabolite == met_name, is.finite(y)) %>%
    filter(!is.na(Treatment))

  if (nrow(dfp) == 0) return(invisible(NULL))

  n_tbl <- dfp %>%
    count(Treatment, name = "n") %>%
    mutate(lbl = paste0(as.character(Treatment), "\n(n=", n, ")"))

  dfp <- dfp %>% left_join(n_tbl %>% select(Treatment, lbl), by = "Treatment")

  p <- ggplot(dfp, aes(x = lbl, y = y)) +
    geom_boxplot(width = 0.55, outlier.shape = NA) +
    geom_jitter(width = 0.12, height = 0, alpha = 0.85, size = 2) +
    labs(
      title = paste0(met_name, " (GC-MS intensity)"),
      x = "Treatment",
      y = "Intensity"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(size = 10)
    )

  out_file <- file.path(plot_dir, paste0(safe_filename(met_name), "_boxplot.svg"))
  ggsave(out_file, plot = p, device = grDevices::svg, width = 6.2, height = 4.4, units = "in")
}

purrr::walk(unique(long$Metabolite), save_met_boxplot_svg)
cat("\nSVG boxplots written to: ", plot_dir, "\n", sep = "")


## 8) PCA plot (SVG) + 95% ellipses + PERMANOVA

pca_svg  <- file.path(out_dir, "GCMS_PCA_by_Treatment_95ellipse.svg")
perm_csv <- file.path(out_dir, "GCMS_PERMANOVA_Treatment_Age_Gender.csv")

pca_df <- dat %>%
  select(Samplecode, Treatment, Age, Gender, Subject, all_of(features_keep)) %>%
  filter(!is.na(Treatment))

perm_meta <- pca_df %>%
  select(Samplecode, Treatment, Age, Gender, Subject) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("Control","ACeD","TCeD")),
    Gender    = factor(Gender),
    Subject   = factor(Subject),
    Age       = suppressWarnings(as.numeric(Age))
  )

X <- as.matrix(pca_df[, features_keep, drop = FALSE])
X <- apply(X, 2, function(z) as.numeric(z))

# Fill NA by feature median
for (j in seq_len(ncol(X))) {
  med <- median(X[, j], na.rm = TRUE)
  X[is.na(X[, j]), j] <- med
}

# Center
Xc <- scale(X, center = TRUE, scale = FALSE)

# Scale
if (pca_scale == "auto") {
  Xc <- scale(Xc, center = FALSE, scale = apply(Xc, 2, sd))
} else if (pca_scale == "pareto") {
  Xc <- scale(Xc, center = FALSE, scale = sqrt(apply(Xc, 2, sd)))
}

p <- prcomp(Xc, center = FALSE, scale. = FALSE)

ve <- (p$sdev^2) / sum(p$sdev^2)
pc1_lab <- paste0("PC1 (", round(ve[1] * 100, 1), "%)")
pc2_lab <- paste0("PC2 (", round(ve[2] * 100, 1), "%)")

scores <- as.data.frame(p$x[, 1:2, drop = FALSE])
scores$Samplecode <- perm_meta$Samplecode
scores$Treatment  <- perm_meta$Treatment
scores$Subject    <- perm_meta$Subject

# Extract time variable for ordering paired samples
scores <- scores %>%
  mutate(
    Time = case_when(
      stringr::str_detect(Samplecode, "W0")  ~ 0,
      stringr::str_detect(Samplecode, "W26") ~ 26,
      TRUE ~ NA_real_
    )
  )

# Ellipse data: only groups with >=3 points
ell_df <- scores %>%
  group_by(Treatment) %>%
  filter(n() >= 3) %>%
  ungroup()

# Subject linking lines: only subjects with >1 point
link_df <- scores %>%
  filter(!is.na(Subject)) %>%
  group_by(Subject) %>%
  filter(n() >= 2) %>%
  arrange(Time, .by_group = TRUE) %>%
  ungroup()

# PERMANOVA:
# Prefer Treatment + Age + Gender if enough complete cases; otherwise Treatment only.
complete_idx <- complete.cases(perm_meta[, c("Treatment","Age","Gender","Subject")])
Xc_perm <- Xc[complete_idx, , drop = FALSE]
perm_meta2 <- perm_meta[complete_idx, , drop = FALSE]

dist_mat <- dist(Xc_perm, method = "euclidean")

set.seed(123)
use_strata <- any(duplicated(perm_meta2$Subject))

if (nrow(perm_meta2) >= 6) {
  if (use_strata) {
    ad <- adonis2(dist_mat ~ Treatment + Age + Gender,
                  data = perm_meta2, permutations = 9999, by = "margin",
                  strata = perm_meta2$Subject)
  } else {
    ad <- adonis2(dist_mat ~ Treatment + Age + Gender,
                  data = perm_meta2, permutations = 9999, by = "margin")
  }
} else {
  if (use_strata) {
    ad <- adonis2(dist_mat ~ Treatment,
                  data = perm_meta2, permutations = 9999,
                  strata = perm_meta2$Subject)
  } else {
    ad <- adonis2(dist_mat ~ Treatment,
                  data = perm_meta2, permutations = 9999)
  }
}

ad_df <- as.data.frame(ad)
ad_df$term <- rownames(ad_df)
write_csv(ad_df, perm_csv)

tr_row <- ad_df %>% filter(term == "Treatment")
t_r2 <- if (nrow(tr_row) == 1) tr_row$R2 else NA_real_
t_p  <- if (nrow(tr_row) == 1) tr_row$`Pr(>F)` else NA_real_

perm_label <- if (is.finite(t_r2) && is.finite(t_p)) {
  paste0("PERMANOVA (Treatment): R²=", sprintf("%.3f", t_r2), ", p=", format(t_p, digits = 3))
} else {
  "PERMANOVA: Treatment row not found"
}

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
    title = "PCA of GC-MS metabolites",
    x = pc1_lab,
    y = pc2_lab
  ) +
  annotate(
    "text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.2,
    label = perm_label, size = 3.3
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(pca_svg, plot = p_pca, device = grDevices::svg, width = 6.4, height = 5.0, units = "in")

cat("\nPERMANOVA saved to: ", perm_csv, "\n", sep = "")
cat("PCA SVG saved to: ", pca_svg, "\n", sep = "")

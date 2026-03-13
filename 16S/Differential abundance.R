#Differential abundance 

#Non-paired

#Uses ANCOM-BC2; this example runs treated vs controls comparison

library(phyloseq)
library(ANCOMBC)
library(dplyr)
library(tidyr)

treated_ps1 <- readRDS("C:/PhD/Sequencing Results/Files/Cant_5%_treated.rds")

outdir <- "C:/PhD/Canterbury Dataset/16S/ancom/treated"
dir.create(outdir, showWarnings = FALSE)

# ----------------------------
set.seed(7)

ps <- treated_ps1

# Number of taxa
n_taxa <- ntaxa(ps)

# Create ASV IDs
new_ids <- paste0("ASV", seq_len(n_taxa))

# Assign as taxa names
taxa_names(ps) <- new_ids

#Set factor level
sd <- sample_data(ps)
sd$CD <- factor(sd$Diagnosed.CD, levels = c(FALSE, TRUE))
sample_data(ps) <- sd

comb_res <- ancombc2(
  data         = ps,
  fix_formula  = "CD + Age + Sex",
  p_adj_method = "fdr",
  group = "CD",
  struc_zero   = TRUE,
  pseudo_sens  = TRUE,
  prv_cut      = 0.15,  
  lib_cut      = 1000,
  alpha        = 0.05,
  n_cl         = 4,
  verbose      = TRUE
)$res

lfc_col <- grep("^lfc_.*CD", names(comb_res), value = TRUE)[1]
se_col  <- grep("^se_.*CD",  names(comb_res), value = TRUE)[1]
p_col   <- sub("^lfc_", "p_",  lfc_col)
q_col   <- sub("^lfc_", "q_",  lfc_col)

if (any(is.na(c(lfc_col, se_col, p_col, q_col))))
  stop("Could not locate the expected lfc / se / p / q columns in pooled ANCOM-BC² output.")

lfc_thresh <- 0.50   # ≈ 1.65-fold change (natural-log scale)

comb_tidy <- comb_res %>%
  transmute(
    taxon        = taxon,
    beta_ancom   = .data[[lfc_col]],
    se_ancom     = .data[[se_col]],
    p_ancom      = .data[[p_col]],
    q_ancom      = .data[[q_col]],
    sig_ancom    = p_ancom < 0.05 & abs(beta_ancom) >= lfc_thresh,
    sign_ancom   = sign(beta_ancom)
  )

comb_sig <- filter(comb_tidy, sig_ancom)


## ── 15) Write output ───────────────────────────────────────────────────────
write.csv(comb_tidy,
          file = file.path(outdir, "combined_ancombc2_results.csv"),
          row.names = FALSE)

write.csv(comb_sig,
          file = file.path(outdir, "combined_ancombc2_sig_only.csv"),
          row.names = FALSE)

message("✓ Pooled ANCOM-BC² complete:",
        "\n  • Full table    → ", file.path(outdir, "combined_ancombc2_results.csv"),
        "\n  • Sig-only taxa → ", file.path(outdir, "combined_ancombc2_sig_only.csv"))

pooled_keep <- comb_sig %>%                              # from pooled ANCOM-BC²
  transmute(
    taxon,
    lfc_pooled = beta_ancom,
    se_pooled  = se_ancom,
    p_pooled   = p_ancom,
    q_pooled   = q_ancom
  )

# Union of significant taxa, keeping all columns from both methods
sig_tbl <- pooled_keep

# Add taxonomy from the original phyloseq object
tax_tbl <- as.data.frame(tax_table(ps))          # matrix → data.frame
tax_tbl$taxon <- rownames(tax_tbl)

sig_tbl <- left_join(sig_tbl, tax_tbl, by = "taxon")

##  ── 17) Add mean relative-abundance columns (overall & by Coeliac status) ─────────────────────

## -- Relative-abundance matrix ------------------------------------------
rel_abund <- transform_sample_counts(ps, function(x) x / sum(x))

abund_mat <- as(otu_table(rel_abund), "matrix")
if (!taxa_are_rows(rel_abund))
  abund_mat <- t(abund_mat)

## Grouping
diag_vec     <- sample_data(rel_abund)$CD
diag_logical <- as.logical(diag_vec)

is_true  <- diag_logical == TRUE
is_false <- diag_logical == FALSE

## Means (with na.rm = TRUE just to be safe)
overall_mean <- rowMeans(abund_mat, na.rm = TRUE)

mean_true <- if (any(is_true)) {
  rowMeans(abund_mat[, is_true, drop = FALSE], na.rm = TRUE)
} else {
  rep(NA_real_, nrow(abund_mat))
}

mean_false <- if (any(is_false)) {
  rowMeans(abund_mat[, is_false, drop = FALSE], na.rm = TRUE)
} else {
  rep(NA_real_, nrow(abund_mat))
}

## Prevalence
overall_prev <- rowMeans(abund_mat > 0, na.rm = TRUE)

prev_true <- if (any(is_true)) {
  rowMeans(abund_mat[, is_true, drop = FALSE] > 0, na.rm = TRUE)
} else {
  rep(NA_real_, nrow(abund_mat))
}

prev_false <- if (any(is_false)) {
  rowMeans(abund_mat[, is_false, drop = FALSE] > 0, na.rm = TRUE)
} else {
  rep(NA_real_, nrow(abund_mat))
}

abund_tbl <- tibble(
  taxon                           = rownames(abund_mat),
  mean_rel_abund_overall          = overall_mean,
  mean_rel_abund_celiac_TRUE      = mean_true,
  mean_rel_abund_celiac_FALSE     = mean_false,
  prev_overall                    = overall_prev,
  prev_celiac_TRUE                = prev_true,
  prev_celiac_FALSE               = prev_false
)

##  -- Merge abundances into sig_tbl ----------------------------
sig_tbl <- left_join(
  sig_tbl,
  abund_tbl,
  by     = "taxon",
  suffix = c("", ".abund")   # ← everything from abund_tbl gets ".abund"
)

##  -- (e) Coalesce & remove duplicates ------------------------------------
for (base in c("mean_rel_abund_overall",
               "mean_rel_abund_celiac_TRUE",
               "mean_rel_abund_celiac_FALSE")) {

  dup <- paste0(base, ".abund")

  if (dup %in% names(sig_tbl)) {              # only if the duplicate exists
    sig_tbl[[base]] <- dplyr::coalesce(sig_tbl[[base]], sig_tbl[[dup]])
    sig_tbl[[dup]]  <- NULL                   # drop the suffixed column
  }
}


## -- 18) Write the final CSV -------------------------------------------------
out_csv <- file.path(outdir, "sig_taxa_combined.csv")
write.csv(sig_tbl, file = out_csv, row.names = FALSE)

message("✓ Combined significant-taxa table (with relative abundances) written to ",
        out_csv)

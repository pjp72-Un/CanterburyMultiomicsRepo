# Aitchison Distance Beta-Diversity Analysis with PERMANOVA & Post-Hoc Tests
# Computes Aitchison distance (CLR-transformed) between samples, runs PERMANOVA group

library(phyloseq)
library(vegan)
library(dplyr)
library(compositions)  # for CLR transformation
library(pairwiseAdonis) # for pairwise post-hoc tests

set.seed(7)

# ── 0) Configuration ────────────────────────────────────────────────

# Output directory
out_dir <- "C:/PhD/Canterbury Dataset/16S/beta_diversity_aitchison"
dir.create(out_dir, recursive = TRUE)

# Load phyloseq object (ps2_untreated_treated or your full object)
ps2_file <- "C:/PhD/Sequencing Results/Files/Cant_5%.rds"
stopifnot(file.exists(ps2_file))
ps <- readRDS(ps2_file)

message("Loaded phyloseq object with ", nsamples(ps), " samples and ", ntaxa(ps), " taxa")

# ── 1) Extract and Prepare Data ────────────────────────────────────

# Extract OTU table and metadata
otu_mat <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) otu_mat <- t(otu_mat)

meta_df <- as(sample_data(ps), "data.frame")
meta_df$SampleID <- rownames(meta_df)

message("OTU matrix: ", nrow(otu_mat), " samples × ", ncol(otu_mat), " taxa")
message("Metadata: ", nrow(meta_df), " rows × ", ncol(meta_df), " columns")

# Ensure required columns exist
required_cols <- c("Subject", "GFD", "Diagnosed.CD", "Age", "Sex")
missing <- setdiff(required_cols, colnames(meta_df))
if (length(missing) > 0) {
  stop("Missing required columns in sample_data: ", paste(missing, collapse = ", "))
}

# ── 2) CLR Transformation and Aitchison Distance ─────────────────────

message("Applying CLR transformation...")

# Add pseudocount to handle zeros (standard approach)
pseudocount <- 1e-6
otu_pseudo <- otu_mat + pseudocount

# CLR transformation: log-ratio of each component to geometric mean
clr_mat <- t(apply(otu_pseudo, 1, function(x) {
  log(x) - mean(log(x))
}))

message("CLR-transformed matrix: ", nrow(clr_mat), " samples × ", ncol(clr_mat), " features")

# Compute Aitchison distance (Euclidean distance on CLR-transformed data)
aitchison_dist <- dist(clr_mat, method = "euclidean")

message("Computed Aitchison distance (", length(aitchison_dist), " pairwise distances)")

# ── 3) Define Groups for PERMANOVA ───────────────────────────────────

# Create a 'Group' variable for PERMANOVA:
#  - "Control" (Diagnosed.CD == FALSE)
#  - "untreated" (Diagnosed.CD == TRUE & GFD == FALSE)
#  - "Treated" (Diagnosed.CD == TRUE & GFD == TRUE)
meta_df <- meta_df %>%
  mutate(Group = case_when(
    Diagnosed.CD == FALSE ~ "Control",
    Diagnosed.CD == TRUE & GFD == FALSE ~ "untreated",
    Diagnosed.CD == TRUE & GFD == TRUE ~ "Treated",
    TRUE ~ "Unknown"
  ))

# Ensure Subject, Group, and Sex are factors
meta_df$Subject <- factor(meta_df$Subject)
meta_df$Group <- factor(meta_df$Group)
meta_df$Sex <- factor(meta_df$Sex)

# Ensure Age is numeric
meta_df$Age <- as.numeric(meta_df$Age)

message("Group breakdown:\n", paste(table(meta_df$Group), collapse = ", "))

# ── 4) PERMANOVA: Global Test by Group ────────────────────────────────

message("\n=== PERMANOVA: Global test by Group ===")

permanova_global <- adonis2(
  aitchison_dist ~ Group + Sex + Age,
  data = meta_df,
  permutations = 999,
  method = "euclidean"  # distances already computed, but specified for clarity
)

print(permanova_global)
permanova_global_df <- as.data.frame(permanova_global)
permanova_global_df$term <- rownames(permanova_global_df)
write.csv(permanova_global_df, file.path(out_dir, "permanova_global_by_group.csv"), row.names = FALSE)
message("Saved PERMANOVA global results")

# ── 5) Pairwise Post-Hoc: Paired untreated vs. Treated ────────────────────

message("\n=== Post-hoc: Paired untreated vs. Treated ===")

# Subset to paired samples (subjects with both untreated and Treated)
paired_samples <- meta_df %>%
  filter(Group %in% c("untreated", "Treated"), !is.na(Subject)) %>%
  group_by(Subject) %>%
  filter(n_distinct(Group) == 2) %>%  # must have both states
  pull(SampleID)

if (length(paired_samples) > 0) {
  meta_paired <- meta_df %>% filter(SampleID %in% paired_samples)
  # Re-factor to remove unused levels
  meta_paired$Subject <- factor(meta_paired$Subject)
  meta_paired$Group <- factor(meta_paired$Group)
  dist_paired <- as.dist(as.matrix(aitchison_dist)[paired_samples, paired_samples])

  # Use blocks parameter to condition on Subject (stratified permutations within subjects)
  # This is equivalent to Condition() but uses the blocks argument
  perm_strata <- how(blocks = meta_paired$Subject, nperm = 999)
  permanova_paired <- adonis2(
    dist_paired ~ Group + Age,
    data = meta_paired,
    permutations = perm_strata
  )
  print(permanova_paired)
  permanova_paired_df <- as.data.frame(permanova_paired)
  permanova_paired_df$term <- rownames(permanova_paired_df)
  write.csv(permanova_paired_df, file.path(out_dir, "permanova_paired_untreated_vs_treated.csv"), row.names = FALSE)
  message("Saved paired PERMANOVA results")

  # Pairwise post-hoc for paired comparison
  pairwise_paired <- pairwise.adonis(
    dist_paired,
    meta_paired$Group,
    p.adjust.m = "BH",
    perm = 999
  )
  write.csv(pairwise_paired, file.path(out_dir, "posthoc_pairwise_paired_untreated_vs_treated.csv"), row.names = FALSE)
  message("Saved paired pairwise post-hoc results")
} else {
  message("Warning: No paired samples found (subjects with both untreated and Treated)")
}

# ── 6) Pairwise Post-Hoc: Unpaired untreated vs. Controls ──────────────────

message("\n=== Post-hoc: Unpaired untreated vs. Controls ===")

untreated_control_samples <- meta_df %>%
  filter(Group %in% c("untreated", "Control")) %>%
  pull(SampleID)

if (length(untreated_control_samples) > 1) {
  meta_ac <- meta_df %>% filter(SampleID %in% untreated_control_samples)
  dist_ac <- as.dist(as.matrix(aitchison_dist)[untreated_control_samples, untreated_control_samples])

  permanova_ac <- adonis2(
    dist_ac ~ Group + Sex + Age,
    data = meta_ac,
    permutations = 999
  )
  print(permanova_ac)
  permanova_ac_df <- as.data.frame(permanova_ac)
  permanova_ac_df$term <- rownames(permanova_ac_df)
  write.csv(permanova_ac_df, file.path(out_dir, "permanova_unpaired_untreated_vs_control.csv"), row.names = FALSE)
  message("Saved untreated vs. Control PERMANOVA results")

  # Pairwise post-hoc
  pairwise_ac <- pairwise.adonis(
    dist_ac,
    meta_ac$Group,
    p.adjust.m = "BH",
    perm = 999
  )
  write.csv(pairwise_ac, file.path(out_dir, "posthoc_pairwise_unpaired_untreated_vs_control.csv"), row.names = FALSE)
  message("Saved untreated vs. Control pairwise post-hoc results")
} else {
  message("Warning: Insufficient samples for untreated vs. Control comparison")
}

# ── 7) Pairwise Post-Hoc: Unpaired Treated vs. Controls ─────────────────

message("\n=== Post-hoc: Unpaired Treated vs. Controls ===")

treated_control_samples <- meta_df %>%
  filter(Group %in% c("Treated", "Control")) %>%
  pull(SampleID)

if (length(treated_control_samples) > 1) {
  meta_tc <- meta_df %>% filter(SampleID %in% treated_control_samples)
  dist_tc <- as.dist(as.matrix(aitchison_dist)[treated_control_samples, treated_control_samples])

  permanova_tc <- adonis2(
    dist_tc ~ Group + Sex+ Age,
    data = meta_tc,
    permutations = 999
  )
  print(permanova_tc)
  permanova_tc_df <- as.data.frame(permanova_tc)
  permanova_tc_df$term <- rownames(permanova_tc_df)
  write.csv(permanova_tc_df, file.path(out_dir, "permanova_unpaired_treated_vs_control.csv"), row.names = FALSE)
  message("Saved Treated vs. Control PERMANOVA results")

  # Pairwise post-hoc
  pairwise_tc <- pairwise.adonis(
    dist_tc,
    meta_tc$Group,
    p.adjust.m = "BH",
    perm = 999
  )
  write.csv(pairwise_tc, file.path(out_dir, "posthoc_pairwise_unpaired_treated_vs_control.csv"), row.names = FALSE)
  message("Saved Treated vs. Control pairwise post-hoc results")
} else {
  message("Warning: Insufficient samples for Treated vs. Control comparison")
}

# ── 8) PCoA Plots (First Two Dimensions) ────────────────────────────────

message("\n=== Generating PCoA Plots ===")

# Create plot directory
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Register Arial font for use in ggplot2
windowsFonts(Arial = windowsFont("Arial"))

library(ggplot2)

# Perform PCoA (Classical Multidimensional Scaling)
pcoa <- cmdscale(aitchison_dist, k = 2, eig = TRUE)
pcoa_coords <- as.data.frame(pcoa$points)
colnames(pcoa_coords) <- c("PC1", "PC2")
pcoa_coords$SampleID <- rownames(pcoa_coords)

# Merge with metadata
pcoa_df <- merge(pcoa_coords, meta_df, by = "SampleID", all.x = TRUE)

# Calculate variance explained
eigenvalues <- pcoa$eig[pcoa$eig > 0]
var_exp <- eigenvalues / sum(eigenvalues) * 100
pc1_var <- var_exp[1]
pc2_var <- var_exp[2]

# Plot 1: All samples colored by Group (with paired lines)
paired_lines <- pcoa_df %>%
  filter(!is.na(Subject)) %>%
  group_by(Subject) %>%
  filter(n_distinct(Group) > 1) %>%  # only subjects with multiple groups
  arrange(Subject, PC1) %>%
  ungroup()

p_all <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
  geom_line(data = paired_lines, aes(x = PC1, y = PC2, group = Subject), color = "gray60", alpha = 0.4, linewidth = 0.5, inherit.aes = FALSE) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal(base_family = "Arial", base_size = 12) +
  labs(
    title = "PCoA of Aitchison Distance - All Samples (Lines Connect Paired Samples)",
    x = paste0("PC1 (", round(pc1_var, 1), "%)"),
    y = paste0("PC2 (", round(pc2_var, 1), "%)")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(plot_dir, "PCoA_all_samples_by_group.svg"), p_all, width = 8, height = 6, device = "svg")
message("Saved: PCoA_all_samples_by_group.svg")

# Plot 2: All samples colored by GFD status (with paired lines)
p_gfd <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = factor(GFD), fill = factor(GFD))) +
  geom_line(data = paired_lines, aes(x = PC1, y = PC2, group = Subject), color = "gray60", alpha = 0.4, linewidth = 0.5, inherit.aes = FALSE) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "#1f77b4", "TRUE" = "#ff7f0e"), 
                      labels = c("FALSE" = "GFD No", "TRUE" = "GFD Yes")) +
  scale_fill_manual(values = c("FALSE" = "#1f77b4", "TRUE" = "#ff7f0e"),
                     labels = c("FALSE" = "GFD No", "TRUE" = "GFD Yes")) +
  theme_minimal(base_family = "Arial", base_size = 12) +
  labs(
    title = "PCoA of Aitchison Distance - By GFD Status (Lines Connect Paired Samples)",
    x = paste0("PC1 (", round(pc1_var, 1), "%)"),
    y = paste0("PC2 (", round(pc2_var, 1), "%)"),
    color = "GFD Status",
    fill = "GFD Status"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(plot_dir, "PCoA_all_samples_by_gfd.svg"), p_gfd, width = 8, height = 6, device = "svg")
message("Saved: PCoA_all_samples_by_gfd.svg")

# Plot 3: All samples colored by CD diagnosis status (with paired lines)
p_cd <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = factor(Diagnosed.CD), fill = factor(Diagnosed.CD))) +
  geom_line(data = paired_lines, aes(x = PC1, y = PC2, group = Subject), color = "gray60", alpha = 0.4, linewidth = 0.5, inherit.aes = FALSE) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "#2ca02c", "TRUE" = "#d62728"),
                      labels = c("FALSE" = "No CD", "TRUE" = "Diagnosed CD")) +
  scale_fill_manual(values = c("FALSE" = "#2ca02c", "TRUE" = "#d62728"),
                     labels = c("FALSE" = "No CD", "TRUE" = "Diagnosed CD")) +
  theme_minimal(base_family = "Arial", base_size = 12) +
  labs(
    title = "PCoA of Aitchison Distance - By CD Diagnosis (Lines Connect Paired Samples)",
    x = paste0("PC1 (", round(pc1_var, 1), "%)"),
    y = paste0("PC2 (", round(pc2_var, 1), "%)"),
    color = "CD Status",
    fill = "CD Status"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(plot_dir, "PCoA_all_samples_by_cd_status.svg"), p_cd, width = 8, height = 6, device = "svg")
message("Saved: PCoA_all_samples_by_cd_status.svg")

message("All PCoA plots saved to: ", plot_dir)

# ── 9) Summary ──────────────────────────────────────────────────────────

message("\n=== Analysis Complete ===")
message("Results saved to: ", out_dir)
message("Files produced:")
message("  - permanova_global_by_group.csv")
message("  - permanova_paired_untreated_vs_treated.csv (if paired samples available)")
message("  - posthoc_pairwise_paired_untreated_vs_treated.csv (if paired samples available)")
message("  - permanova_unpaired_untreated_vs_control.csv (if samples available)")
message("  - posthoc_pairwise_unpaired_untreated_vs_control.csv (if samples available)")
message("  - permanova_unpaired_treated_vs_control.csv (if samples available)")
message("  - posthoc_pairwise_unpaired_treated_vs_control.csv (if samples available)")

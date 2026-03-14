#Prevalence filtering of phyloseq objects for use in downstreat machine learing

library(phyloseq)

# Minimal: make 15% and 30% prevalence-filtered species CLR+scaled tables

in_rds  <- "C://PhD//Sequencing Results//Files//Cant_5%.rds"
out_dir <- "C://PhD//Sequencing Results//Files"

Cant_16S <- readRDS(in_rds)

make_species_table <- function(ps, prev_prop = 0.15, pseudocount = 1e-6) {

  # 1) Agglomerate to Species
  ps_sp <- tax_glom(ps, taxrank = "Species")

  # 2) Create Genus_Species name (only when Species exists)
  tax <- as.data.frame(tax_table(ps_sp))
  tax$Genus_Species <- ifelse(
    is.na(tax$Species) | tax$Species == "",
    NA_character_,
    paste(tax$Genus, tax$Species, sep = "_")
  )
  tax_table(ps_sp) <- as.matrix(tax)

  # 3) Keep only taxa with a valid Species name
  keep_named <- !is.na(tax_table(ps_sp)[, "Genus_Species"]) &
                tax_table(ps_sp)[, "Genus_Species"] != ""
  taxa_to_keep <- taxa_names(ps_sp)[keep_named]
  ps_sp <- prune_taxa(taxa_to_keep, ps_sp)

  # 4) Prevalence filter (fraction of samples with abundance > 0)
  otu <- as(otu_table(ps_sp), "matrix")
  if (!taxa_are_rows(ps_sp)) otu <- t(otu)

  prev <- colMeans(otu > 0)
  keep_taxa <- prev >= prev_prop
  taxa_to_keep <- colnames(otu)[keep_taxa]
  ps_sp <- prune_taxa(taxa_to_keep, ps_sp)

  # 5) OTU matrix with samples as rows
  otu <- as(otu_table(ps_sp), "matrix")
  if (taxa_are_rows(ps_sp)) otu <- t(otu)

  # 6) Set taxa names = Genus_Species (ensure uniqueness)
  gs <- as.character(tax_table(ps_sp)[, "Genus_Species"])
  colnames(otu) <- make.unique(gs)

  # 7) Relative abundance -> CLR -> scale
  rel <- sweep(otu, 1, rowSums(otu), "/")
  clr <- log(rel + pseudocount) - rowMeans(log(rel + pseudocount))
  scaled <- scale(clr, center = TRUE, scale = TRUE)

  as.matrix(scaled)
}

# Run for 15% and 30%
for (p in c(0.15, 0.30)) {
  mat <- make_species_table(Cant_16S, prev_prop = p)

  tag <- paste0("prev", as.integer(p * 100))
  saveRDS(mat, file.path(out_dir, paste0("Cant_DIABLO_species_", tag, ".rds")))
  write.csv(mat, file.path(out_dir, paste0("Cant_DIABLO_species_", tag, ".csv")), row.names = TRUE)

  cat("Wrote:", tag, "| samples =", nrow(mat), "| species =", ncol(mat), "\n")
}

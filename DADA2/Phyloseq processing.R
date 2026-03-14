
#Preperation of phyloseq objects

#CantUF = unfiltered dada2 phyloseq (the exported phyloseq object produced from Dada2 pipeline.R)

library(phyloseq)

Cant <- readRDS("C:/PhD/Sequencing Results/Files/CantUF.rds")

## 1. Read metadata from CSV
meta <- read.csv(
  "C:/PhD/Sequencing Results/Files/meta.csv",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# use sample_id as rownames (must be a column in meta.csv)
rownames(meta) <- meta$sample_id

## 2. Pull current sample_data from Cant
sd <- as.data.frame(sample_data(Cant))

## 3. Add Age, Sex, Subject from meta, matched by sample (rownames)
sd[c("Age", "Sex", "Subject")] <- meta[rownames(sd), c("Age", "Sex", "Subject")]

## 4. Set types
sd$Age     <- as.numeric(sd$Age)
sd$Sex     <- factor(sd$Sex)
sd$Subject <- factor(sd$Subject)

## 5. Put back into Cant
sample_data(Cant) <- sample_data(sd)


table(tax_table(Cant)[, "Phylum"], exclude = NULL)
ps0 <- subset_taxa(Cant, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
prevdf = apply(X = otu_table(ps0),
                 MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(ps0),
                      tax_table(ps0))
        
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# Define phyla to filter
filterPhyla = c("Chlorobi", "Chloroflexi","CPR2","Cyanobacteria","FBP","Miscellaneous_Euryarchaeotic_Group(MEG)",
                "GN02","OD1","RBG-1_(Zixibacteria)","SR1_(Absconditabacteria)","TM6_(Dependentiae)"
                ,"WS2","WS6")

# Filter entries with unidentified Phylum.
ps1 = subset_taxa(ps0, !Phylum %in% filterPhyla)
ps1

prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps0),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

#  Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps0)
prevalenceThreshold

## [1] 18

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps0)

# untreated diet status (not on GFD), both controls and CeD cases
ps1_untreated <- subset_samples(ps1,
                             GFD == FALSE & !is.na(Diagnosed.CD))

ps1_untreated <- prune_taxa(taxa_sums(ps1_untreated) > 0, ps1_untreated)

ps1_treated <- subset_samples(
  ps1,
  (GFD == TRUE  & Diagnosed.CD == TRUE)  |   # treated CeD
  (GFD == FALSE & Diagnosed.CD == FALSE)     # controls not on GFD
)
ps1_treated <- prune_taxa(taxa_sums(ps1_treated) > 0, ps1_treated)


# pull metadata
meta <- as.data.frame(sample_data(ps1))

# define untreated and treated CD samples
is_untreated_cd  <- with(meta, Diagnosed.CD == TRUE & GFD == FALSE & !is.na(Subject))
is_treated_cd <- with(meta, Diagnosed.CD == TRUE & GFD == TRUE  & !is.na(Subject))

# subjects that have BOTH an untreated and a treated sample
untreated_subjects  <- unique(meta$Subject[is_untreated_cd])
treated_subjects <- unique(meta$Subject[is_treated_cd])
paired_subjects  <- intersect(untreated_subjects, treated_subjects)

# keep only CD samples from those subjects, both states (untreated + treated)
paired_samples <- rownames(meta)[
  meta$Subject %in% paired_subjects &
    meta$Diagnosed.CD == TRUE &
    meta$GFD %in% c(FALSE, TRUE)
]

ps1_untreated_treated <- prune_samples(paired_samples, ps1)
ps1_untreated_treated <- prune_taxa(taxa_sums(ps1_untreated_treated) > 0, ps1_untreated_treated)


#5% objects

# untreated diet status (not on GFD), both controls and CD cases
ps2_untreated <- subset_samples(ps2,
                             GFD == FALSE & !is.na(Diagnosed.CD))

ps2_untreated <- prune_taxa(taxa_sums(ps2_untreated) > 0, ps2_untreated)

ps2_treated <- subset_samples(
  ps2,
  (GFD == TRUE  & Diagnosed.CD == TRUE)  |   # treated CD
  (GFD == FALSE & Diagnosed.CD == FALSE)     # controls not on GFD
)
ps2_treated <- prune_taxa(taxa_sums(ps2_treated) > 0, ps2_treated)


library(phyloseq)

# pull metadata
meta <- as.data.frame(sample_data(ps2))

# define untreated and treated CD samples
is_untreated_cd  <- with(meta, Diagnosed.CD == TRUE & GFD == FALSE & !is.na(Subject))
is_treated_cd <- with(meta, Diagnosed.CD == TRUE & GFD == TRUE  & !is.na(Subject))

# subjects that have BOTH an untreated and a treated sample
untreated_subjects  <- unique(meta$Subject[is_untreated_cd])
treated_subjects <- unique(meta$Subject[is_treated_cd])
paired_subjects  <- intersect(untreated_subjects, treated_subjects)

# keep only CD samples from those subjects, both states (untreated + treated)
paired_samples <- rownames(meta)[
  meta$Subject %in% paired_subjects &
    meta$Diagnosed.CD == TRUE &
    meta$GFD %in% c(FALSE, TRUE)
]

ps2_untreated_treated <- prune_samples(paired_samples, ps2)
ps2_untreated_treated <- prune_taxa(taxa_sums(ps2_untreated_treated) > 0, ps2_untreated_treated)



saveRDS(ps2, "C:/PhD/Sequencing Results/Files/Cant_5%.rds")
saveRDS(ps1, "C:/PhD/Sequencing Results/Files/Cant_uf.rds")
saveRDS(ps1_untreated, "C:/PhD/Sequencing Results/Files/Cant_untreated.rds")
saveRDS(ps1_treated, "C:/PhD/Sequencing Results/Files/Cant_treated.rds")
saveRDS(ps1_untreated_treated, "C:/PhD/Sequencing Results/Files/Cant_untreated_treated.rds")

saveRDS(ps2_untreated, "C:/PhD/Sequencing Results/Files/Cant_5%_untreated.rds")
saveRDS(ps2_treated, "C:/PhD/Sequencing Results/Files/Cant_5%_treated.rds")
saveRDS(ps2_untreated_treated, "C:/PhD/Sequencing Results/Files/Cant_5%_untreated_treated.rds")

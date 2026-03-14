#Phylogentic Diversity

#Faiths PD

library(phyloseq)
library(picante)
library(lme4)
library(lmerTest)
library(dplyr)

ps <- readRDS("C:/PhD/Canterbury Dataset/objects/Cant_UF_tree.RDS")

# untreated diet status (not on GFD), both controls and CD cases
ps0_untreated <- subset_samples(ps,
                             GFD == FALSE & !is.na(Diagnosed.CD))

ps0_treated <- subset_samples(
  ps,
  (GFD == TRUE  & Diagnosed.CD == TRUE)  |   # treated CD
  (GFD == FALSE & Diagnosed.CD == FALSE)     # controls not on GFD
)
ps0_treated <- prune_taxa(taxa_sums(ps0_treated) > 0, ps0_treated)

# pull metadata
meta <- as.data.frame(sample_data(ps))

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

ps0_untreated_treated <- prune_samples(paired_samples, ps)
ps1_untreated_treated <- prune_taxa(taxa_sums(ps1_untreated_treated) > 0, ps1_untreated_treated)

## 1. Compute Faith's PD + metadata for a phyloseq object ----
compute_faith_pd_df <- function(ps) {
    comm <- as(otu_table(ps), "matrix")
    if (taxa_are_rows(ps)) comm <- t(comm)
    
    tree <- phy_tree(ps)
    
    pd_res <- picante::pd(comm, tree, include.root = FALSE)
    
    meta <- as(sample_data(ps), "data.frame")
    meta$SampleID <- rownames(meta)
    meta$FaithPD  <- pd_res$PD[match(rownames(meta), rownames(pd_res))]
    
    meta
}

## 2. Build data frames for each object & write CSVs ----
df_untreated  <- compute_faith_pd_df(ps0_untreated)
df_treated <- compute_faith_pd_df(ps0_treated)
df_pair    <- compute_faith_pd_df(ps0_untreated_treated)

write.csv(df_untreated,  "ps0_untreated_FaithPD_metadata.csv",         row.names = FALSE)
write.csv(df_treated, "ps0_treated_FaithPD_metadata.csv",        row.names = FALSE)
write.csv(df_pair,    "ps0_untreated_treated_FaithPD_metadata.csv", row.names = FALSE)

## 3. Prep for modelling (set types) ----
prep_for_model <- function(df) {
    df %>%
        mutate(
            Diagnosed.CD = factor(Diagnosed.CD),GFD = factor(GFD),
            Sex          = factor(Sex),
            Subject      = factor(Subject),
            Age          = as.numeric(Age)
        )
}

df_untreated_m  <- prep_for_model(df_untreated)
df_treated_m <- prep_for_model(df_treated)
df_pair_m    <- prep_for_model(df_pair)

## 4. Models ----
# Unpaired: simple linear models
m_untreated  <- lm(FaithPD ~ Diagnosed.CD + Age + Sex, data = df_untreated_m)
m_treated <- lm(FaithPD ~ Diagnosed.CD + Age + Sex, data = df_treated_m)

# Paired: mixed model with random intercept for Subject
m_pair    <- lmer(FaithPD ~ GFD + Age + (1 | Subject),
                  data = df_pair_m)

## 5. Extract estimates, SE, p-values ----
tidy_model <- function(mod, model_name) {
    s <- summary(mod)
    coefs <- as.data.frame(s$coefficients)
    coefs$term  <- rownames(coefs)
    coefs$model <- model_name
    
    if (!"df" %in% colnames(coefs)) {
        df_resid <- s$df[2]
        coefs$df <- df_resid
    }
    
    coefs <- coefs[, c("model", "term", "Estimate", "Std. Error", "df",
                       "t value", "Pr(>|t|)")]
    rownames(coefs) <- NULL
    coefs
}

res_untreated  <- tidy_model(m_untreated,  "ps0_untreated")
res_treated <- tidy_model(m_treated, "ps0_treated")
res_pair    <- tidy_model(m_pair,    "ps0_untreated_treated")

model_results <- rbind(res_untreated, res_treated, res_pair)

write.csv(model_results,
          "FaithPD_models_all_objects.csv",
          row.names = FALSE)

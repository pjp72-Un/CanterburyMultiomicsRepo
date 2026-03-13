# Phylum top-10 models for ps2_untreated / ps2_treated / ps2_untreated_treated
# Computes relative abundance at the Phylum level, finds the top-10 phyla by mean relative abundance
# Runs the following models for each phylum:
#  - unpaired objects (ps2_untreated, ps2_treated): lm(Abund ~ Diagnosed.CD + Age + Sex)
#  - paired object (ps2_untreated_treated): lmer(Abund ~ GFD + Age + (1 | Subject))
#
# Outputs: "Phylum_top10_models_all_objects.csv"

library(phyloseq)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(broom)

# Set working directory 
setwd("C:/PhD/R code/cant_study")

# --- 0) Load objects ---
ps2_untreated_file        <- "C:/PhD/Sequencing Results/Files/Cant_5%_untreated.rds"
ps2_treated_file       <- "C:/PhD/Sequencing Results/Files/Cant_5%_treated.rds"
ps2_untreated_treated_file <- "C:/PhD/Sequencing Results/Files/Cant_5%_untreated_treated.rds"

stopifnot(file.exists(ps2_untreated_file), file.exists(ps2_treated_file), file.exists(ps2_untreated_treated_file))

ps2_untreated         <- readRDS(ps2_untreated_file)
ps2_treated        <- readRDS(ps2_treated_file)
ps2_untreated_treated <- readRDS(ps2_untreated_treated_file)

# --- helper: Phylum-level relative abundance (phyloseq -> data.frame with samples x phylum) ---
get_phylum_rel_abund_df <- function(ps) {
  # Agglomerate to Phylum
  ps_ph <- tax_glom(ps, taxrank = "Phylum")

  # drop taxa with zero counts
  ps_ph <- prune_taxa(taxa_sums(ps_ph) > 0, ps_ph)

  # relative abundance per sample
  ps_ph_rel <- transform_sample_counts(ps_ph, function(x) x / sum(x))

  # turn into long data frame: SampleID, Phylum, Abundance
  comm <- as.data.frame(otu_table(ps_ph_rel))
  if (taxa_are_rows(ps_ph_rel)) {
    comm <- t(comm)
  }
  comm <- as.data.frame(comm)

  # tax_glom already aggregates to Phylum level; ensure taxa names are readable
  # (Phylum names) so downstream pivot_longer produces meaningful labels.
  tt <- as.data.frame(tax_table(ps_ph_rel))
  if ("Phylum" %in% colnames(tt)) {
    phylum_vec <- as.character(tt[, "Phylum"])
    phylum_vec[is.na(phylum_vec) | phylum_vec == ""] <- "Unidentified"
    # make unique names when duplicates exist (tax_glom should already combine by Phylum,
    # but this is a safe guard)
    taxa_names(ps_ph_rel) <- make.unique(phylum_vec, sep = "_")
    # update comm to use the newly-assigned taxa names
    comm <- as.data.frame(otu_table(ps_ph_rel))
    if (taxa_are_rows(ps_ph_rel)) comm <- t(comm)
    comm <- as.data.frame(comm)
  } else {
    warning("tax_table does not contain a 'Phylum' column; using OTU/taxa IDs as column names")
  }

  # add SampleID column
  comm$SampleID <- rownames(comm)
  long <- tidyr::pivot_longer(comm, cols = -all_of("SampleID"), names_to = "Phylum", values_to = "Abundance")

  # Merge metadata
  md <- as.data.frame(sample_data(ps_ph_rel))
  md$SampleID <- rownames(md)

  df <- dplyr::left_join(long, md, by = "SampleID")
  df
}

# --- choose top 10 phyla by mean abundance in the paired object (ps2_untreated_treated) ---
phylum_df_pair <- get_phylum_rel_abund_df(ps2_untreated_treated)

## If you want a custom set of phyla to test, list them here. If left NULL,
## the script will fall back to selecting the top-10 by mean rel-abundance.
selected_phyla <- c(
  "Verrucomicrobia",
  "Actinobacteria",
  "Proteobacteria",
  "Firmicutes",
  "Bacteroidetes"
)

# compute top-10 as a fallback (kept for compatibility)
top10_phyla <- phylum_df_pair %>%
  group_by(.data$Phylum) %>%
  summarise(mean_abund = mean(.data$Abundance, na.rm = TRUE)) %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = 10) %>%
  pull(.data$Phylum)

# Decide which phyla to test: prefer explicitly selected_phyla when provided
if (!is.null(selected_phyla) && length(selected_phyla) > 0) {
  # only keep those that exist in the paired object
  phyla_to_test <- intersect(selected_phyla, unique(phylum_df_pair$Phylum))
  missing_phyla <- setdiff(selected_phyla, phyla_to_test)
  if (length(missing_phyla) > 0) {
    message("Warning: the following requested phyla were not found in ps2_untreated_treated and will be skipped: ", paste(missing_phyla, collapse = ", "))
  }
} else {
  phyla_to_test <- top10_phyla
}

if (length(phyla_to_test) == 0) stop("No phyla available to test — check your selected_phyla and taxonomy.")

message("Phyla selected for testing: ", paste(phyla_to_test, collapse = ", "))

## Subset the ps objects to only keep taxa in phyla_to_test, so downstream
if (!is.null(phyla_to_test) && length(phyla_to_test) > 0) {
  for (obj_name in c("ps2_untreated", "ps2_treated", "ps2_untreated_treated")) {
    ps_obj <- get(obj_name)
    tt <- as.data.frame(tax_table(ps_obj))
    if (!"Phylum" %in% colnames(tt)) {
      warning("Object ", obj_name, " does not have Phylum in tax_table; skipping taxonomic subset for this object")
      next
    }
    # Subset taxa to only those phyla
    ps_sub <- subset_taxa(ps_obj, Phylum %in% phyla_to_test)
    ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
    assign(obj_name, ps_sub)
    message("Subset ", obj_name, " -> kept taxa in phyla: ", paste(unique(as.data.frame(tax_table(ps_sub))$Phylum), collapse = ", "))
  }
  # Rebuild phylum_df_pair to reflect the subset
  phylum_df_pair <- get_phylum_rel_abund_df(ps2_untreated_treated)
}

# Diagnostics collector: record why models fail / skip
diagnostics <- list()

# --- modeling helper: fit model for one phylum for a given data frame and object type ---
fit_phylum_model <- function(df, object_type = c("untreated","treated","paired"), phylum) {
  object_type <- match.arg(object_type)

  subdf <- df %>% filter(.data$Phylum == phylum)

  # guard against zero-variance / empty
  if (nrow(subdf) == 0) {
    diagnostics <<- append(diagnostics, list(list(Phylum = phylum, Object = object_type,
                                                  reason = "no rows for this phylum in object")))
    return(NULL)
  }

  # Transform abundance for modelling. Proportions (0..1): use arcsin-sqrt transform
  subdf <- subdf %>% mutate(Abund_tr = asin(sqrt(.data$Abundance)))

  # Prepare variables
  subdf$Age <- as.numeric(subdf$Age)
  if ("Diagnosed.CD" %in% names(subdf)) subdf$Diagnosed.CD <- factor(subdf$Diagnosed.CD)
  if ("Sex" %in% names(subdf)) subdf$Sex <- factor(subdf$Sex)
  if ("GFD" %in% names(subdf)) subdf$GFD <- factor(subdf$GFD)
  if ("Subject" %in% names(subdf)) subdf$Subject <- factor(subdf$Subject)

  if (object_type %in% c("untreated", "treated")) {
    # Unpaired linear model: Abund_tr ~ Diagnosed.CD + Age + Sex
    if (!("Diagnosed.CD" %in% names(subdf))) {
      diagnostics <<- append(diagnostics, list(list(Phylum = phylum, Object = object_type,
                                                    reason = "Diagnosed.CD not present in sample_data for this object")))
      return(NULL)
    }

    fm <- try(lm(Abund_tr ~ Diagnosed.CD + Age + Sex, data = subdf), silent = TRUE)
    if (inherits(fm, "try-error")) {
      diagnostics <<- append(diagnostics, list(list(Phylum = phylum, Object = object_type,
                                                    reason = paste0("lm failed: ", conditionMessage(attr(fm, "condition"))))))
      return(NULL)
    }

    # Use broom::tidy() for a consistent, safe tidy output
    td <- broom::tidy(fm)
    df_resid <- try(df.residual(fm), silent = TRUE)
    if (inherits(df_resid, "try-error")) df_resid <- NA
    td$df <- df_resid

    coefs <- td %>%
      transmute(
        model = object_type,
        term = .data$term,
        Estimate = .data$estimate,
        Std.Error = .data$std.error,
        df = .data$df,
        t.value = .data$statistic,
        P_value = .data$p.value
      ) %>%
      mutate(Phylum = phylum, Object = object_type)
    rownames(coefs) <- NULL
    return(coefs)

  } else if (object_type == "paired") {
    # Paired: lmer(Abund_tr ~ GFD + Age + (1 | Subject))
    if (!("GFD" %in% names(subdf)) || !("Subject" %in% names(subdf))) {
      diagnostics <<- append(diagnostics, list(list(Phylum = phylum, Object = object_type,
                                                    reason = "GFD or Subject missing in sample_data for this object")))
      return(NULL)
    }

    if (length(unique(subdf$Subject)) < 2L) {
      diagnostics <<- append(diagnostics, list(list(Phylum = phylum, Object = object_type,
                                                    reason = "Subject has fewer than 2 unique levels; cannot fit random intercept")))
      return(NULL)
    }
    if (length(unique(subdf$GFD)) < 2L) {
      diagnostics <<- append(diagnostics, list(list(Phylum = phylum, Object = object_type,
                                                    reason = "GFD has fewer than 2 levels; no contrast to test")))
      return(NULL)
    }

    fm <- try(lmer(Abund_tr ~ GFD + Age + (1 | Subject), data = subdf), silent = TRUE)
    if (inherits(fm, "try-error")) {
      diagnostics <<- append(diagnostics, list(list(Phylum = phylum, Object = object_type,
                                                    reason = paste0("lmer failed: ", conditionMessage(attr(fm, "condition"))))))
      return(NULL)
    }

    td <- tryCatch(
      broom::tidy(fm, effects = "fixed"),
      error = function(e) {
        s <- summary(fm)
        coefs <- as.data.frame(s$coefficients)
        # coefs often has rownames = terms and columns such as
        # 'Estimate', 'Std. Error', 't value', 'Pr(>|t|)'
        coefs$term <- rownames(coefs)
        stat <- if ("t value" %in% colnames(coefs)) coefs[["t value"]] else NA_real_
        pval <- if ("Pr(>|t|)" %in% colnames(coefs)) coefs[["Pr(>|t|)"]] else NA_real_
        data.frame(term = coefs$term,
                   estimate = coefs[["Estimate"]],
                   std.error = coefs[["Std. Error"]],
                   statistic = stat,
                   p.value = pval,
                   stringsAsFactors = FALSE)
      }
    )

    # If tidy output lacks p-values but has a test statistic, approximate via normal
    if (!"p.value" %in% names(td) && "statistic" %in% names(td)) {
      td$p.value <- 2 * pnorm(-abs(td$statistic))
    }
    if (!"p.value" %in% names(td)) td$p.value <- NA_real_
    if (!"statistic" %in% names(td)) td$statistic <- NA_real_
    td$df <- NA

    coefs <- td %>%
      transmute(
        model = object_type,
        term = .data$term,
        Estimate = .data$estimate,
        Std.Error = .data$std.error,
        df = .data$df,
        t.value = .data$statistic,
        P_value = .data$p.value
      ) %>%
      mutate(Phylum = phylum, Object = object_type)

    rownames(coefs) <- NULL
    return(coefs)
  }

  NULL
}

# ---  Build per-object dataframes ---
df_untreated  <- get_phylum_rel_abund_df(ps2_untreated)
df_treated <- get_phylum_rel_abund_df(ps2_treated)
df_pair    <- get_phylum_rel_abund_df(ps2_untreated_treated)

# --- Run models across top10 and return tidy results ---
results_list <- list()
for (ph in top10_phyla) {
  # Fit on untreated object
  r_untreated <- fit_phylum_model(df_untreated,  object_type = "untreated",  phylum = ph)
  r_treated <- fit_phylum_model(df_treated, object_type = "treated", phylum = ph)
  r_pair <- fit_phylum_model(df_pair, object_type = "paired", phylum = ph)

  results_list <- c(results_list, list(r_untreated, r_treated, r_pair))
}

# combine and clean
results_df <- do.call(rbind, lapply(results_list, function(x) { if (is.null(x)) NULL else x }))

if (is.null(results_df) || nrow(results_df) == 0) {
  # If all models failed, write a diagnostics table so the user can see why
  diag_df <- do.call(rbind, lapply(diagnostics, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
  if (!is.null(diag_df) && nrow(diag_df) > 0) {
    diag_file <- "C:/PhD/Canterbury Dataset/16S/phylum/Phylum_top10_models_diagnostics.csv"
    write.csv(diag_df, diag_file, row.names = FALSE)
    stop(sprintf("No model results produced — models may have failed for all phyla. Diagnostics written to: %s", diag_file))
  } else {
    stop("No model results produced — models may have failed for all phyla. No diagnostics available.")
  }
}

# reorder columns
results_df <- results_df %>%
  relocate(Phylum, Object, model, term, Estimate, Std.Error, df, t.value, P_value)

# Add multiple-testing corrections
# - P_adj_BH_global: BH correction across all model-terms
# - P_adj_BH_by_term: BH correction applied within each model term (e.g. all 'GFD' tests across phyla)
# - P_adj_bonf_global: global Bonferroni correction across all tests
results_df <- results_df %>%
  mutate(P_value = as.numeric(P_value)) %>%
  # global corrections
  mutate(
    P_adj_BH_global = p.adjust(P_value, method = "BH"),
    P_adj_bonf_global = p.adjust(P_value, method = "bonferroni")
  ) %>%
  group_by(term) %>%
  mutate(P_adj_BH_by_term = p.adjust(P_value, method = "BH")) %>%
  ungroup()

## Add a user-friendly significance column based on the term-wise BH adjustment
# Mark 'sig' when P_adj_BH_by_term < 0.05, otherwise 'ns'
results_df <- results_df %>%
  mutate(sig_by_term = ifelse(!is.na(P_adj_BH_by_term) & P_adj_BH_by_term < 0.05, "sig", "ns"))

# write outputs
write.csv(results_df, "C:/PhD/Canterbury Dataset/16S/phylum/Phylum_top10_models_all_objects.csv", row.names = FALSE)

message("Done. Results saved to Phylum_top10_models_all_objects.csv")

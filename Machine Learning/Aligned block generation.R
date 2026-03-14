# Prepare aligned blocks for downstream glmnet pipeline

# This script keeps only data preparation steps:
# - Load microbiome, metabolome, metadata
# - Canonicalize IDs and intersect samples
# - Optional transform
# - Impute + near-zero variance filter
# - Save aligned blocks as RDS for glmnet workflow

options(stringsAsFactors = FALSE)


# USER CONFIG

#Example run at species level

project_dir <- "C:/PhD/Canterbury Dataset/DIABLO"
setwd(project_dir)

microbiome_file <- "Cant_DIABLO_species_Transposed_30.csv"
metabolome_file <- "metabolome_merged_gcms_primary.csv"
metadata_file   <- "Metadata_good.csv"

metadata_sample_id_col <- NULL
metadata_group_col     <- "Status"

microbiome_transform <- "none" # options: none, clr, auto
metabolome_transform <- "none" # options: none, log1p
clr_pseudocount      <- 1e-6

max_na_prop            <- 0.20
min_feature_nonmissing <- 5
near_zero_var_tol      <- 1e-8

out_dir <- "glmnet_inputs"
out_rds <- "aligned_blocks_for_glmnet.rds"


# HELPERS

canon_id <- function(x) {
  toupper(gsub("[^A-Za-z0-9]+", "", trimws(as.character(x))))
}

safe_as_numeric <- function(v) {
  raw <- v
  v <- trimws(as.character(v))
  v_lower <- tolower(v)

  miss <- is.na(v_lower) | v_lower %in% c(
    "", "na", "n/a", "nan", "null",
    "nd", "n.d.", "not detected", "below detection", "belowdetection",
    "bdl", "trace",
    "#num!", "#div/0!", "#value!", "#ref!"
  )
  v[miss] <- NA
  v <- gsub("^[<>]=?\\s*", "", v)
  v <- gsub(",", "", v)
  v <- gsub("[^0-9eE+\\-\\.]", "", v)

  num <- suppressWarnings(as.numeric(v))
  num[!is.finite(num)] <- NA

  nonmissing_raw <- !is.na(raw) & trimws(as.character(raw)) != ""
  coerced <- nonmissing_raw & is.na(num)

  attr(num, "n_coercion") <- sum(coerced, na.rm = TRUE)
  attr(num, "n_nonmissing_raw") <- sum(nonmissing_raw, na.rm = TRUE)
  num
}

read_omics_csv <- function(path,
                           block_name = NULL,
                           max_na_prop = 0.20,
                           min_feature_nonmissing = 5) {
  if (!file.exists(path)) stop("File not found: ", path)
  if (is.null(block_name)) block_name <- basename(path)

  df_raw <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(df_raw) < 2) stop("Expected at least 2 columns in: ", path)

  ids_can <- canon_id(df_raw[[1]])
  X_raw <- df_raw[-1]
  rownames(X_raw) <- ids_can

  conv <- lapply(X_raw, safe_as_numeric)
  coercion_counts <- sapply(conv, function(x) attr(x, "n_coercion"))
  nonmissing_counts <- sapply(conv, function(x) attr(x, "n_nonmissing_raw"))

  total_coercion <- sum(coercion_counts, na.rm = TRUE)
  total_nonmissing <- sum(nonmissing_counts, na.rm = TRUE)
  total_coercion_prop <- if (total_nonmissing > 0) total_coercion / total_nonmissing else 0

  if (total_coercion > 0) {
    cat("\n--- PARSE DIAGNOSTICS:", block_name, "---\n")
    cat("Total coerced-to-NA entries:", total_coercion,
        sprintf("(%.2f%% of non-missing raw entries)\n", 100 * total_coercion_prop))
  }

  X <- as.data.frame(conv, check.names = FALSE)
  rownames(X) <- ids_can

  na_prop <- colMeans(is.na(X))
  nonmiss <- colSums(!is.na(X))
  drop <- (na_prop > max_na_prop) | (nonmiss < min_feature_nonmissing)

  cat("\n--- NA PROFILE AFTER PARSING:", block_name, "---\n")
  cat("Samples:", nrow(X), "Features:", ncol(X),
      "| NA% median:", round(100 * median(na_prop), 2),
      "| max:", round(100 * max(na_prop), 2),
      "| drop:", sum(drop), "\n")

  if (sum(drop) > 0) X <- X[, !drop, drop = FALSE]
  if (anyDuplicated(rownames(X))) stop("Non-unique canonical sample IDs after canon_id() in ", block_name)

  as.matrix(X)
}

median_impute <- function(X) {
  X2 <- X
  for (j in seq_len(ncol(X2))) {
    v <- X2[, j]
    if (anyNA(v)) {
      med <- median(v, na.rm = TRUE)
      if (!is.finite(med)) med <- 0
      v[is.na(v)] <- med
      X2[, j] <- v
    }
  }
  X2
}

drop_near_zero_var <- function(X, tol = 1e-8) {
  sds <- apply(X, 2, sd)
  keep <- is.finite(sds) & (sds > tol)
  list(X = X[, keep, drop = FALSE], dropped = colnames(X)[!keep])
}

clr_transform <- function(X, pseudocount = 1e-6) {
  X2 <- X
  X2[!is.finite(X2)] <- NA
  X2[X2 < 0] <- NA
  X2 <- median_impute(X2)
  X2 <- X2 + pseudocount
  logX <- log(X2)
  logX - rowMeans(logX)
}

choose_micro_transform <- function(X, mode = "auto") {
  if (mode == "none") return("none")
  if (mode == "clr") return("clr")
  if (any(X < 0, na.rm = TRUE)) return("none")
  "clr"
}

auto_pick_col <- function(df, candidates) {
  nms <- colnames(df)
  low <- tolower(nms)
  for (c in candidates) {
    hit <- which(low == tolower(c))
    if (length(hit) == 1) return(nms[hit])
  }
  NULL
}

read_metadata <- function(path, sample_col = NULL, group_col = "Group") {
  if (!file.exists(path)) stop("Metadata file not found: ", path)
  md <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)

  if (is.null(sample_col)) {
    sample_col <- auto_pick_col(md, c("sampleid", "sample_id", "sample", "id", "subject", "participant"))
    if (is.null(sample_col)) sample_col <- colnames(md)[1]
  }
  if (!(sample_col %in% colnames(md))) stop("Sample ID column not found in metadata: ", sample_col)

  if (!(group_col %in% colnames(md))) {
    g2 <- auto_pick_col(md, c("group", "condition", "status", "class"))
    if (!is.null(g2)) group_col <- g2
  }
  if (!(group_col %in% colnames(md))) stop("Group column not found in metadata: ", group_col)

  ids_can <- canon_id(md[[sample_col]])
  grp <- droplevels(as.factor(md[[group_col]]))

  cat("\n--- METADATA COLUMN SELECTION ---\n")
  cat("Sample ID column:", sample_col, "\n")
  cat("Group column    :", group_col, "\n")
  cat("Group levels    :", paste(levels(grp), collapse = ", "), "\n")

  data.frame(sample_can = ids_can, Group = grp, stringsAsFactors = FALSE)
}

                              
# MAIN

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Preparing aligned blocks for glmnet ===\n")
cat("microbiome_file:", microbiome_file, "\n")
cat("metabolome_file:", metabolome_file, "\n")
cat("metadata_file  :", metadata_file, "\n")
cat("out_dir        :", out_dir, "\n")

micro_X <- read_omics_csv(
  microbiome_file,
  block_name = "Microbiome",
  max_na_prop = max_na_prop,
  min_feature_nonmissing = min_feature_nonmissing
)

metab_X <- read_omics_csv(
  metabolome_file,
  block_name = "Metabolome(merged)",
  max_na_prop = max_na_prop,
  min_feature_nonmissing = min_feature_nonmissing
)

md <- read_metadata(
  metadata_file,
  sample_col = metadata_sample_id_col,
  group_col = metadata_group_col
)

ids <- Reduce(intersect, list(rownames(micro_X), rownames(metab_X), md$sample_can))
cat("\n--- INTERSECTION (Micro+Metab+Meta) ---\n")
cat("Count:", length(ids), "\n")
if (length(ids) < 6) stop("Too few matched samples after ID alignment.")

micro_X <- micro_X[ids, , drop = FALSE]
metab_X <- metab_X[ids, , drop = FALSE]
group <- droplevels(md$Group[match(ids, md$sample_can)])

cat("\n--- GROUP COUNTS ---\n")
print(table(group))

micro_mode <- choose_micro_transform(micro_X, microbiome_transform)
cat("\n--- PREPROCESSING ---\n")
cat("Microbiome transform:", micro_mode, "\n")
cat("Metabolome transform:", metabolome_transform, "\n")

if (micro_mode == "clr") {
  micro_X <- clr_transform(micro_X, pseudocount = clr_pseudocount)
}
if (metabolome_transform == "log1p") {
  metab_X <- log1p(pmax(metab_X, 0))
}

micro_X <- median_impute(micro_X)
metab_X <- median_impute(metab_X)

mzv <- drop_near_zero_var(micro_X, tol = near_zero_var_tol)
micro_X <- mzv$X

nzv <- drop_near_zero_var(metab_X, tol = near_zero_var_tol)
metab_X <- nzv$X

out_path <- file.path(out_dir, out_rds)
saveRDS(
  list(Microbiome = micro_X, Metabolome = metab_X, Group = group),
  out_path
)

cat("\nSaved aligned blocks to:", out_path, "\n")
cat("Final dimensions:\n")
cat("  Microbiome:", nrow(micro_X), "x", ncol(micro_X), "\n")
cat("  Metabolome:", nrow(metab_X), "x", ncol(metab_X), "\n")
cat("  Group levels:", paste(levels(group), collapse = ", "), "\n")

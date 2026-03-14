# Merge NMR + GCMS into a single "Metabolome" block (small-n safe)

# - Drops high-missing features, median-imputes the rest


# OUTPUTS (in out_dir):
#   - metabolome_merged_gcms_primary.csv (samples x features; GCMS preferred on overlaps)
#   - metabolome_qc_summary.txt     (counts, NA stats, dropped features)

options(stringsAsFactors = FALSE)


# USER CONFIG

nmr_file  <- "nmr_metabolites.csv"
gcms_file <- "GC_MS_metabolites.csv"

# How to merge samples across platforms:
#   "intersection": keep only samples present in BOTH NMR and GCMS
#   "union": keep all samples, leaving NAs for missing platform measurements
merge_mode <- "intersection"

# Feature missingness filtering 
max_na_prop          <- 0.20  # drop features with >20% missing
min_feature_nonmissing <- 5   # drop features with fewer than 5 non-missing samples

# Imputation
impute_method <- "median"     

# Optional transforms (use if your inputs are NOT already log2)
assume_log2_inputs <- TRUE    # set FALSE to apply log2(x + pseudo)

# Output folder
out_dir <- "metabolome_merge_outputs"


# HELPERS

canon_id <- function(x) {
  # Canonicalise sample IDs so e.g. "WK0-1002" and "wk0_1002" match
  toupper(gsub("[^A-Za-z0-9]+", "", as.character(x)))
}

is_missing_token <- function(x_lower) {
  x_lower %in% c(
    "", "na", "n/a", "nan", "null",
    "nd", "n.d.", "notdetected", "not detected", "belowdetection", "below detection",
    "bdl", "trace",
    "#num!", "#div/0!", "#value!", "#ref!"
  )
}

safe_numeric <- function(x) {

  x0 <- as.character(x)
  x1 <- trimws(x0)

  # track which entries were genuinely present pre-clean
  present_raw <- !(is.na(x1) | x1 == "")

  x_lower <- tolower(x1)
  miss <- is.na(x_lower) | is_missing_token(x_lower)
  x1[miss] <- NA

  # comparator numeric strings: "<=0.1", "> 5" -> "0.1", "5"
  x1 <- gsub("^[<>]=?\\s*", "", x1)

  # remove thousands separators
  x1 <- gsub(",", "", x1)

  # keep only numeric characters (incl exponent)
  x_clean <- gsub("[^0-9eE+.-]", "", x1)

  num <- suppressWarnings(as.numeric(x_clean))

  # coercion = raw was present, but we ended with NA after parsing
  coerced <- present_raw & !is.na(x0) & !(tolower(trimws(x0)) %in% tolower(c(
    "", "na", "n/a", "nan", "null",
    "nd", "n.d.", "not detected", "below detection", "bdl", "trace",
    "#num!", "#div/0!", "#value!", "#ref!"
  ))) & is.na(num)

  list(num = num, coerced = coerced, raw = x0)
}

read_metabolite_csv <- function(path,
                               max_na_prop = 0.20,
                               min_feature_nonmissing = 5) {

  if (!file.exists(path)) stop("File not found: ", path)

  df <- read.csv(path, check.names = FALSE)
  if (ncol(df) < 2) stop("Expected at least 2 columns in: ", path)

  # First column is assumed sample IDs (common after write.csv with row.names=TRUE)
  sid_col <- 1
  sample_ids <- as.character(df[[sid_col]])
  sample_ids <- trimws(sample_ids)

  # Handle accidental duplicates early
  if (anyDuplicated(sample_ids)) {
    dups <- unique(sample_ids[duplicated(sample_ids)])
    stop("Duplicate sample IDs in ", path, ": ", paste(head(dups, 10), collapse = ", "))
  }

  df_feat <- df[, -sid_col, drop = FALSE]

  # Robust numeric conversion + diagnostics
  coerced_total <- 0L
  coerced_by_col <- integer(ncol(df_feat))
  names(coerced_by_col) <- colnames(df_feat)

  out_cols <- vector("list", ncol(df_feat))
  names(out_cols) <- colnames(df_feat)

  for (j in seq_len(ncol(df_feat))) {
    res <- safe_numeric(df_feat[[j]])
    out_cols[[j]] <- res$num
    coerced_by_col[[j]] <- sum(res$coerced, na.rm = TRUE)
    coerced_total <- coerced_total + coerced_by_col[[j]]
  }

  mat <- as.data.frame(out_cols, check.names = FALSE)
  rownames(mat) <- sample_ids

  # Use raw feature names so cross-platform overlaps can be resolved explicitly
  colnames(mat) <- trimws(colnames(mat))

  if (anyDuplicated(colnames(mat))) {
    dups <- unique(colnames(mat)[duplicated(colnames(mat))])
    stop("Duplicate metabolite names in ", basename(path), ": ", paste(head(dups, 10), collapse = ", "))
  }

  # Parse diagnostics
  nonmissing_raw <- 0L
  for (j in seq_len(ncol(df_feat))) {
    x <- as.character(df_feat[[j]])
    nonmissing_raw <- nonmissing_raw + sum(!(is.na(x) | trimws(x) == ""), na.rm = TRUE)
  }
  cat("\n--- PARSE DIAGNOSTICS:", basename(path), "---\n")
  cat("Total coerced-to-NA entries:", coerced_total,
      sprintf("(%.2f%% of non-missing raw entries)\n", 100 * coerced_total / max(1, nonmissing_raw)))

  if (coerced_total > 0) {
    top <- sort(coerced_by_col, decreasing = TRUE)
    top <- top[top > 0]
    show_n <- min(10, length(top))
    if (show_n > 0) {
      cat("Top columns by coercion count (up to 10):\n")
      for (k in seq_len(show_n)) {
        cat("  -", names(top)[k], ":", top[[k]], "\n")
      }
      warning("Minor non-numeric values were coerced to NA in ", basename(path),
              " (top column: ", names(top)[1], ", n=", top[[1]], "). Proceeding with NA + imputation.")
    }
  }

  # NA profile + feature filtering
  na_prop <- colMeans(is.na(mat))
  nonmiss <- colSums(!is.na(mat))

  drop <- (na_prop > max_na_prop) | (nonmiss < min_feature_nonmissing)
  cat("\n--- NA PROFILE AFTER PARSING:", basename(path), "---\n")
  cat("Samples:", nrow(mat), "Features:", ncol(mat),
      "| NA% median:", round(100 * median(na_prop), 2),
      "| max:", round(100 * max(na_prop), 2),
      "| drop:", sum(drop), "\n")

  dropped_features <- names(drop)[drop]
  if (sum(drop) > 0) {
    mat <- mat[, !drop, drop = FALSE]
  }

  list(
    mat = as.matrix(mat),
    dropped_features = dropped_features,
    sample_ids = sample_ids,
    sample_ids_canon = canon_id(sample_ids)
  )
}

align_to_canon <- function(mat, canon_ids, target_canon) {
  # Align matrix rows to target canonical IDs (fills missing rows with NA)
  stopifnot(nrow(mat) == length(canon_ids))
  idx <- match(target_canon, canon_ids)
  out <- matrix(NA_real_, nrow = length(target_canon), ncol = ncol(mat),
                dimnames = list(target_canon, colnames(mat)))
  hit <- which(!is.na(idx))
  out[hit, ] <- mat[idx[hit], , drop = FALSE]
  out
}

impute_median <- function(mat) {
  # Column-wise median imputation
  out <- mat
  for (j in seq_len(ncol(out))) {
    v <- out[, j]
    if (anyNA(v)) {
      med <- median(v, na.rm = TRUE)
      if (!is.finite(med)) med <- 0
      v[is.na(v)] <- med
      out[, j] <- v
    }
  }
  out
}

drop_zero_variance <- function(mat, tol = 1e-8) {
  sds <- apply(mat, 2, sd)
  keep <- is.finite(sds) & (sds > tol)
  list(mat = mat[, keep, drop = FALSE], dropped = colnames(mat)[!keep])
}

maybe_log2 <- function(mat) {
  # Apply log2 transform if not already log2.
  # Uses a small pseudo-count based on the smallest positive value.
  if (assume_log2_inputs) return(mat)

  vmin <- suppressWarnings(min(mat, na.rm = TRUE))
  if (!is.finite(vmin)) stop("Cannot log-transform: all values are NA.")
  if (vmin <= 0) {
    # Shift up so minimum is slightly above 0
    shift <- abs(vmin) + 1e-6
    mat2 <- mat + shift
  } else {
    mat2 <- mat
  }

  pos <- mat2[is.finite(mat2) & mat2 > 0]
  pseudo <- if (length(pos) > 0) min(pos) / 2 else 1e-6
  log2(mat2 + pseudo)
}


# MAIN

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
qc_path <- file.path(out_dir, "metabolome_qc_summary.txt")
qc_con <- file(qc_path, open = "wt")
sink(qc_con, split = TRUE)

cat("=== Metabolome merge QC summary ===\n")
cat("NMR file :", nmr_file, "\n")
cat("GCMS file:", gcms_file, "\n")
cat("merge_mode:", merge_mode, "\n\n")

nmr  <- read_metabolite_csv(nmr_file,
                            max_na_prop = max_na_prop,
                            min_feature_nonmissing = min_feature_nonmissing)

gcms <- read_metabolite_csv(gcms_file,
                            max_na_prop = max_na_prop,
                            min_feature_nonmissing = min_feature_nonmissing)

# Check canonical ID uniqueness (critical for safe alignment)
if (anyDuplicated(nmr$sample_ids_canon)) stop("NMR has non-unique canonical sample IDs after canon_id().")
if (anyDuplicated(gcms$sample_ids_canon)) stop("GCMS has non-unique canonical sample IDs after canon_id().")

if (merge_mode == "intersection") {
  target_canon <- intersect(nmr$sample_ids_canon, gcms$sample_ids_canon)
} else if (merge_mode == "union") {
  target_canon <- union(nmr$sample_ids_canon, gcms$sample_ids_canon)
} else {
  stop("merge_mode must be 'intersection' or 'union'")
}

cat("\n--- SAMPLE ALIGNMENT ---\n")
cat("NMR samples :", length(nmr$sample_ids_canon), "\n")
cat("GCMS samples:", length(gcms$sample_ids_canon), "\n")
cat("Target samples (", merge_mode, "):", length(target_canon), "\n", sep = "")

# Align each platform to the target set
nmr_aligned  <- align_to_canon(nmr$mat,  nmr$sample_ids_canon,  target_canon)
gcms_aligned <- align_to_canon(gcms$mat, gcms$sample_ids_canon, target_canon)

# Merge features with GCMS priority on overlaps
nmr_feats <- colnames(nmr_aligned)
gcms_feats <- colnames(gcms_aligned)

overlap_feats <- intersect(nmr_feats, gcms_feats)
nmr_only_feats <- setdiff(nmr_feats, gcms_feats)
gcms_only_feats <- setdiff(gcms_feats, nmr_feats)

cat("\n--- FEATURE OVERLAP (GCMS PRIORITY) ---\n")
cat("NMR features :", length(nmr_feats), "\n")
cat("GCMS features:", length(gcms_feats), "\n")
cat("Overlapping features (using GCMS values):", length(overlap_feats), "\n")
cat("NMR-only features :", length(nmr_only_feats), "\n")
cat("GCMS-only features:", length(gcms_only_feats), "\n")

metab_raw <- cbind(
  gcms_aligned[, c(gcms_only_feats, overlap_feats), drop = FALSE],
  nmr_aligned[, nmr_only_feats, drop = FALSE]
)

# Optional transform (if inputs not already log2)
metab_raw <- maybe_log2(metab_raw)

# Feature NA filter after merge (optional but useful for union mode)
na_prop_merged <- colMeans(is.na(metab_raw))
nonmiss_merged <- colSums(!is.na(metab_raw))
drop_m <- (na_prop_merged > max_na_prop) | (nonmiss_merged < min_feature_nonmissing)
cat("\n--- MERGED METABOLOME FEATURE FILTER ---\n")
cat("Merged features:", ncol(metab_raw), "\n")
cat("Dropping merged features:", sum(drop_m), "\n")
metab_filt <- metab_raw[, !drop_m, drop = FALSE]

# Impute
if (impute_method == "median") {
  metab_imp <- impute_median(metab_filt)
} else {
  stop("Only median imputation is implemented in this basic script.")
}

# Drop zero-variance
zv <- drop_zero_variance(metab_imp)
metab_final <- zv$mat

cat("\n--- FINAL METABOLOME BLOCK ---\n")
cat("Samples :", nrow(metab_final), "\n")
cat("Features:", ncol(metab_final), "\n")
cat("Dropped zero-variance features:", length(zv$dropped), "\n")

# Restore readable rownames:

rownames(metab_final) <- target_canon

out_csv <- file.path(out_dir, "metabolome_merged_gcms_primary.csv")
write.csv(metab_final, out_csv, quote = TRUE)
cat("\nWrote:", out_csv, "\n")
cat("Wrote:", qc_path, "\n")

sink()
close(qc_con)

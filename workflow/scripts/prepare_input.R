# Snakemake shim: prepare PAM50 + intrinsic inputs from the expression matrix and
# sample sheet. Applies the (generalised) sample filter and optional replicate
# deduplication, maps to the 50 PAM50 genes, and subsets + median-centres the
# intrinsic gene matrix. Scientific logic is in pam50perou.

suppressPackageStartupMessages(library(pam50perou))

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message"); sink(log, type = "output")

cfg <- snakemake@config
id_col <- cfg$inputs$sample_id_column

# --- read inputs ---
expr <- as.matrix(read.table(snakemake@input[["expression"]],
                             sep = "\t", header = TRUE, row.names = 1, check.names = FALSE))
samples <- read.table(snakemake@input[["samples"]],
                      sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
if (!id_col %in% colnames(samples)) {
  stop(sprintf("sample_id_column '%s' not found in sample sheet", id_col))
}
rownames(samples) <- samples[[id_col]]
message(sprintf("inputs: %d genes x %d samples; %d sample-sheet rows",
                nrow(expr), ncol(expr), nrow(samples)))

# --- generalised sample filter (default: all samples) ---
filt <- cfg$filter
if (!is.null(filt$column) && !is.na(filt$column) && nzchar(filt$column)) {
  if (!filt$column %in% colnames(samples)) {
    stop(sprintf("filter.column '%s' not in sample sheet", filt$column))
  }
  include <- unlist(filt$include)
  keep <- as.character(samples[[filt$column]]) %in% as.character(include)
  samples <- samples[keep, , drop = FALSE]
  message(sprintf("filter %s in {%s}: %d samples retained",
                  filt$column, paste(include, collapse = ","), nrow(samples)))
} else {
  message("no filter: using all samples")
}

# --- optional replicate deduplication ---
dd <- cfg$deduplicate
if (!is.null(dd) && !is.null(dd$by) && !is.na(dd$by) && nzchar(dd$by)) {
  if (!dd$by %in% colnames(samples)) stop(sprintf("deduplicate.by '%s' not in sample sheet", dd$by))
  sel <- dd$selections
  grp <- as.character(samples[[dd$by]])
  keep_ids <- character(0)
  for (g in unique(grp)) {
    members <- rownames(samples)[grp == g]
    pref <- if (!is.null(sel)) sel[[g]] else NULL
    if (!is.null(pref) && pref %in% members) {
      keep_ids <- c(keep_ids, pref)
    } else {
      if (!is.null(pref)) message(sprintf("  dedup: selection %s=%s not eligible; using default", g, pref))
      keep_ids <- c(keep_ids, sort(members)[1])
    }
  }
  samples <- samples[keep_ids, , drop = FALSE]
  message(sprintf("deduplicated by %s: %d samples retained", dd$by, nrow(samples)))
}

# --- align expression to the retained samples ---
keep_samples <- intersect(rownames(samples), colnames(expr))
if (length(keep_samples) < nrow(samples)) {
  message(sprintf("WARNING: %d sample(s) in sheet not found in expression matrix",
                  nrow(samples) - length(keep_samples)))
}
samples <- samples[keep_samples, , drop = FALSE]
expr <- expr[, keep_samples, drop = FALSE]

# --- PAM50 gene matrix (HGNC symbol -> pcrID, subset 50) ---
pam50_expr <- map_to_pam50(expr)
message(sprintf("PAM50 genes: %d/50", nrow(pam50_expr)))

# --- intrinsic gene matrix (subset + median-centre) ---
intr <- intrinsic_genes()
int_mask <- rownames(expr) %in% intr$Name
int_expr <- expr[int_mask, , drop = FALSE]
int_expr <- int_expr[!duplicated(rownames(int_expr)), , drop = FALSE]
message(sprintf("intrinsic genes: %d", nrow(int_expr)))
int_expr_medctr <- median_center(int_expr)

# --- write outputs ---
write_tsv <- function(mat, path, key = "gene") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  df <- data.frame(key = rownames(mat), mat, check.names = FALSE)
  colnames(df)[1] <- key
  write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE)
}
write_tsv(pam50_expr, snakemake@output[["pam50_expr"]])
write_tsv(int_expr, snakemake@output[["int_expr"]])
write_tsv(int_expr_medctr, snakemake@output[["int_expr_medctr"]])

dir.create(dirname(snakemake@output[["samples_filtered"]]), recursive = TRUE, showWarnings = FALSE)
write.table(samples, snakemake@output[["samples_filtered"]],
            sep = "\t", row.names = FALSE, quote = FALSE)
message("prepare_input: done")

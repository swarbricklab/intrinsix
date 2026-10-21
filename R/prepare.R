# Input preparation helpers: map an expression matrix from HGNC symbols to the
# classifier's pcrID identifiers, and median-centre a matrix (for the intrinsic
# heatmap path). Mirrors the subset/relabel logic of prepare_pam50_input.py but
# operates on a plain matrix (the standalone workflow's TSV entry contract).

#' Map an expression matrix to the 50 PAM50 classifier genes
#'
#' Relabels rows from HGNC symbols to the classifier's `pcrID`, subsets to the 50
#' PAM50 genes, and resolves duplicate symbols by keeping the first occurrence
#' (matching the current pipeline).
#'
#' @param expr Numeric matrix/data.frame, genes x samples. Row names are HGNC
#'   gene symbols.
#' @param annotation Annotation table with `pcrID` and `GeneName` columns;
#'   defaults to [pam50_annotation()].
#' @return Numeric matrix of the PAM50 genes (pcrID row names) x samples. Warns
#'   if fewer than 50 of the PAM50 genes are present.
#' @export
map_to_pam50 <- function(expr, annotation = pam50_annotation()) {
  expr <- as.matrix(expr)
  # Map on the modern HGNC GeneName, and also accept the legacy pcrID itself as a
  # key, so both modern- and legacy-symbol inputs resolve.
  sym2pcr <- stats::setNames(annotation$pcrID, annotation$GeneName)
  sym2pcr <- c(sym2pcr, stats::setNames(annotation$pcrID, annotation$pcrID))
  sym2pcr <- sym2pcr[!duplicated(names(sym2pcr))]

  present <- rownames(expr) %in% names(sym2pcr)
  sub <- expr[present, , drop = FALSE]
  rownames(sub) <- sym2pcr[rownames(sub)]

  # keep-first on duplicate pcrID
  sub <- sub[!duplicated(rownames(sub)), , drop = FALSE]

  n <- nrow(sub)
  if (n < 50L) {
    missing <- setdiff(annotation$pcrID, rownames(sub))
    warning(sprintf("map_to_pam50: only %d/50 PAM50 genes present; missing: %s",
                    n, paste(missing, collapse = ", ")), call. = FALSE)
  }
  storage.mode(sub) <- "double"
  sub
}

#' Select the PAM50 genes using precomputed classifier IDs
#'
#' The drift-proof alternative to [map_to_pam50()]: instead of re-matching HGNC
#' symbols (which rename over time, e.g. ORC6L -> ORC6), use a per-gene classifier
#' ID assigned upstream from a stable identifier (e.g. Ensembl gene IDs). This is
#' how the source pipeline fed the classifier: an h5ad `var$pam50_symbol` column
#' holding the `pcrID` for the 50 PAM50 genes (and empty for all others).
#'
#' @param expr Numeric matrix, genes x samples (any row names).
#' @param pam50_ids Character vector, one entry per row of `expr` (same order),
#'   giving each gene's classifier `pcrID`; empty/`NA` for non-PAM50 genes.
#' @return Numeric matrix of the PAM50 genes (pcrID row names) x samples. Warns if
#'   fewer than 50 are present.
#' @export
select_pam50_by_id <- function(expr, pam50_ids) {
  expr <- as.matrix(expr)
  if (length(pam50_ids) != nrow(expr)) {
    stop("pam50_ids must have one entry per row of expr.", call. = FALSE)
  }
  ids <- as.character(pam50_ids)
  keep <- !is.na(ids) & nzchar(ids) & !(tolower(ids) %in% c("nan", "na"))
  sub <- expr[keep, , drop = FALSE]
  rownames(sub) <- ids[keep]
  sub <- sub[!duplicated(rownames(sub)), , drop = FALSE]
  if (nrow(sub) < 50L) {
    warning(sprintf("select_pam50_by_id: only %d/50 PAM50 genes present.", nrow(sub)),
            call. = FALSE)
  }
  storage.mode(sub) <- "double"
  sub
}

#' Median-centre a matrix (row-wise)
#'
#' Subtracts each row's median. Used to prepare the intrinsic expression matrix
#' for clustering.
#'
#' @param x Numeric matrix, genes x samples.
#' @return Row-median-centred matrix.
#' @export
median_center <- function(x) {
  x <- as.matrix(x)
  medians <- apply(x, 1, stats::median, na.rm = TRUE)
  sweep(x, 1, medians, "-")
}

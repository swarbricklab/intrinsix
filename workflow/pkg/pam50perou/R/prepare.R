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
  sym2pcr <- stats::setNames(annotation$pcrID, annotation$GeneName)

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

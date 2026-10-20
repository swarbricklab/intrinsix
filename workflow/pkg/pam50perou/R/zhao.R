# Zhao4 ER/HER2 subtype-specific quantile centering. Faithful port of
# vendor/PAM50_zhao_normalization.R (Alec Lobanov, Perou Lab). For each clinical
# (ER x HER2) subtype group, each PAM50 gene is centred by subtracting its
# SIGMA-defined quantile computed within that group.

#' Zhao4 subtype-specific quantile centering
#'
#' @param expr Numeric matrix, PAM50 genes x samples (pcrID row names), log2-scale.
#' @param subtypes Named character vector mapping sample ID -> clinical subtype.
#'   Names must match `colnames(expr)`; values must match columns of `sigma`
#'   (e.g. "ERpos_HER2neg", "ERneg_HER2neg", "HER2pos_ERpos", "HER2pos_ERneg").
#' @param sigma SIGMA quantile table (genes x subtypes); defaults to [sigma_table()].
#' @return Zhao-normalised expression matrix (genes x samples), columns in the
#'   same grouped order the pipeline produces.
#' @export
zhao_center <- function(expr, subtypes, sigma = sigma_table()) {
  expr <- as.matrix(expr)
  storage.mode(expr) <- "double"

  quantile_centering <- function(expr_matrix, gene_quantile) {
    res <- expr_matrix
    for (name in rownames(expr_matrix)) {
      q <- gene_quantile[name, ]
      q_value <- stats::quantile(as.numeric(expr_matrix[name, ]), probs = q, na.rm = TRUE)
      res[name, ] <- expr_matrix[name, ] - q_value
    }
    res
  }

  centered_matrix <- list()
  for (subtype in unique(subtypes)) {
    samples <- names(subtypes)[subtypes == subtype]
    gene_expr_subtype <- expr[, colnames(expr) %in% samples, drop = FALSE]
    percentile_group <- sigma[, subtype, drop = FALSE]
    centered_matrix[[subtype]] <- quantile_centering(gene_expr_subtype, percentile_group)
  }

  normalized_matrix <- do.call(cbind, centered_matrix)
  # do.call(cbind, <named list>) prefixes column names with "<subtype>."; strip it
  # back to the bare sample ID (matches the pipeline's tstrsplit(.,"\\.")[[2]]).
  colnames(normalized_matrix) <- sub("^[^.]*\\.", "", colnames(normalized_matrix))
  normalized_matrix
}

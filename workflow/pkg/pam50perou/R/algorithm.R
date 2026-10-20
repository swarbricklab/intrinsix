# Algorithm functions. These are STUBS pending the port of the vendored Perou
# scripts (see vendor/perou-lab/ and vendor/*.R) into explicit-argument functions.
# Tracked as the next milestone in DESIGN.md. Each stub documents the target
# signature and its source so the port is mechanical.

.not_ported <- function(src) {
  stop(sprintf(
    "pam50perou: not yet ported from %s. See DESIGN.md (algorithm packaging milestone).",
    src
  ), call. = FALSE)
}

#' Map an expression matrix to PAM50 classifier gene identifiers
#'
#' Relabels rows from HGNC symbols to the classifier's pcrID via
#' [pam50_annotation()], subsets to the 50 PAM50 genes, and resolves duplicate
#' symbols (keep-first, matching the current pipeline).
#'
#' @param expr Numeric matrix/data.frame, genes x samples, log2-scale.
#' @param annotation Annotation table; defaults to [pam50_annotation()].
#' @return Expression matrix restricted to the 50 PAM50 genes, pcrID-labelled.
#' @export
map_to_pam50 <- function(expr, annotation = pam50_annotation()) {
  .not_ported("prepare_pam50_input.py (subset/relabel logic)")
}

#' Zhao4 subtype-specific quantile centering
#'
#' Per clinical (ER x HER2) subtype, subtracts the SIGMA-defined quantile of each
#' gene. Source: vendor/PAM50_zhao_normalization.R.
#'
#' @param expr PAM50 expression matrix (genes x samples).
#' @param subtypes Named vector mapping sample -> clinical subtype (SIGMA columns).
#' @param sigma SIGMA quantile table; defaults to [sigma_table()].
#' @return Zhao-normalised expression matrix.
#' @export
zhao_center <- function(expr, subtypes, sigma = sigma_table()) {
  .not_ported("vendor/PAM50_zhao_normalization.R")
}

#' PAM50 subtype prediction
#'
#' Spearman rank correlation to the published centroids, with ROR-S / ROR-P risk
#' scores. Wraps the Perou `sspPredict`/`readarray`/`standardize` helpers
#' (vendor/perou-lab/subtypePrediction_{functions,distributed}.R), refactored to
#' take explicit arguments instead of reading globals.
#'
#' @param expr Normalised PAM50 expression matrix (genes x samples).
#' @param calibration Perou calibrationParameters (-1 = none).
#' @param collapse "mean" or "iqr".
#' @param has_clinical Logical; TRUE enables ROR-C (requires T-stage).
#' @param centroids Centroid table; defaults to [pam50_centroids()].
#' @return A data frame: centroid distances, Call, Confidence, ROR scores, etc.
#' @export
predict_pam50 <- function(expr, calibration = -1, collapse = "mean",
                          has_clinical = FALSE, centroids = pam50_centroids()) {
  .not_ported("vendor/perou-lab/subtypePrediction_distributed.R")
}

#' Cluster3.0 centroid clustering of the intrinsic matrix
#'
#' Source: vendor/clustering_scripts.R (centroid_clustering, xcluster2r).
#'
#' @param expr_medctr Median-centred intrinsic expression matrix.
#' @param outdir Directory for the Cluster3.0 .cdt/.atr/.gtr files.
#' @return Paths to the cluster files.
#' @export
centroid_clustering <- function(expr_medctr, outdir) {
  .not_ported("vendor/clustering_scripts.R")
}

#' Identify claudin-low samples
#'
#' Source: vendor/identify_claudin_low.R.
#'
#' @param cluster_files Output of [centroid_clustering()].
#' @param expr_medctr Median-centred intrinsic expression matrix.
#' @return Character vector of claudin-low sample IDs.
#' @export
identify_claudin_low <- function(cluster_files, expr_medctr) {
  .not_ported("vendor/identify_claudin_low.R")
}

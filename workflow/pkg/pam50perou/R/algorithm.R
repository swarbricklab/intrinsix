# Remaining un-ported algorithm functions (intrinsic clustering / claudin-low).
# These are STUBS pending the port of the vendored clustering scripts; tracked as
# the next milestone in DESIGN.md. The PAM50 classification path (map_to_pam50 ->
# zhao_center -> predict_pam50) is fully ported and verified; see R/prepare.R,
# R/zhao.R, R/predict.R.

.not_ported <- function(src) {
  stop(sprintf(
    "pam50perou: not yet ported from %s. See DESIGN.md (clustering milestone).",
    src
  ), call. = FALSE)
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

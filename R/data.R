#' Accessors for bundled PAM50 reference data
#'
#' These read the reference files shipped in `inst/extdata`, so no external
#' checkout (git submodule, DVC pull) is required at runtime.
#'
#' @name reference_data
#' @return A data frame (or named list) of the requested reference data.
NULL

#' @describeIn reference_data PAM50 gene annotation (HGNC symbol -> classifier pcrID).
#' @export
pam50_annotation <- function() {
  utils::read.delim(
    system.file("extdata", "pam50_annotation.txt", package = "intrinsix"),
    stringsAsFactors = FALSE
  )
}

#' @describeIn reference_data Published PAM50 subtype centroids (50 genes x 5 subtypes).
#' @export
pam50_centroids <- function() {
  utils::read.table(
    system.file("extdata", "pam50_centroids.txt", package = "intrinsix"),
    sep = "\t", header = TRUE, row.names = 1
  )
}

#' @describeIn reference_data SIGMA quantile table for Zhao4 normalisation.
#' @export
sigma_table <- function() {
  utils::read.table(
    system.file("extdata", "SIGMA.txt", package = "intrinsix"),
    sep = "\t", header = TRUE, row.names = 1
  )
}

#' @describeIn reference_data Intrinsic subtype gene list.
#' @export
intrinsic_genes <- function() {
  utils::read.delim(
    system.file("extdata", "intrinsic_subtype_gene_list.txt", package = "intrinsix"),
    stringsAsFactors = FALSE
  )
}

#' @describeIn reference_data Training set (visualisation only; 220 arrays + 12 normal).
#' @export
training_set <- function() {
  system.file("extdata", "trainingset_220arrays_50g.txt", package = "intrinsix")
}

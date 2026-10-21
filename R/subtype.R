# High-level interactive interface: one call from a raw expression matrix to a
# rich result object, plus S3 methods (print/summary/plot/as.data.frame).

#' PAM50 subtyping, end to end
#'
#' One-call interface: maps an expression matrix to the PAM50 genes, optionally
#' applies Zhao4 subtype-specific normalisation, classifies, and (optionally) runs
#' intrinsic clustering + claudin-low identification. Returns a [pam50_result]
#' object you can `print()`, `summary()`, `plot()`, or coerce with
#' `as.data.frame()`.
#'
#' @param object The expression input, as any of: a numeric matrix (genes x
#'   samples); a [SummarizedExperiment::SummarizedExperiment] /
#'   `SingleCellExperiment` (use `layer` to pick the assay); an `anndataR`
#'   `AnnData`; or a path to an `.h5ad` file (needs the `anndataR` package, R >=
#'   4.5). For objects carrying sample annotations (SE/AnnData/h5ad), the sample
#'   sheet is taken from them unless `samples` is given. Row names must be HGNC
#'   gene symbols (set `map_symbols = FALSE` if already pcrID). For `cluster =
#'   TRUE` the input must contain the intrinsic genes (a full matrix).
#' @param samples Optional data frame, one row per sample, carrying the clinical
#'   subtype column used for Zhao normalisation. Overrides annotations from an
#'   SE/AnnData/h5ad input.
#' @param subtype_col Name of the clinical (ER x HER2) subtype column in
#'   `samples`. Required when `zhao = TRUE`.
#' @param layer For SE/AnnData/h5ad inputs, the assay/layer to use (e.g.
#'   `"log2_uqn"`); defaults to the first assay / `X`.
#' @param zhao Logical; apply Zhao4 normalisation. Defaults to `TRUE` when
#'   `subtype_col` is supplied, else `FALSE`.
#' @param map_symbols Logical; map HGNC symbols to pcrID and subset to the 50
#'   PAM50 genes first (default `TRUE`).
#' @param calibration Passed to [predict_pam50()]. Defaults to `-1` when `zhao`
#'   is on (input already centred) and `NA` (median-centre) otherwise.
#' @param collapse Duplicate-gene collapse method, `"mean"` or `"iqr"`.
#' @param cluster Logical; also run intrinsic clustering + claudin-low (adds a
#'   `Call_CLOW` column and an `$intrinsic` component). Default `FALSE`.
#' @param id_col Sample-ID column in `samples` (default `"sample_id"`).
#' @return A [pam50_result] object: `$calls` (the scores/calls data frame),
#'   `$classify` (prediction intermediates for plotting), `$params`, and, when
#'   `cluster = TRUE`, `$intrinsic`.
#' @examples
#' \dontrun{
#' res <- pam50_subtype(expr, samples = meta, subtype_col = "clinical_subtype")
#' res
#' summary(res)
#' plot(res)
#' }
#' @export
pam50_subtype <- function(object, samples = NULL, subtype_col = NULL, layer = NULL,
                          zhao = !is.null(subtype_col), map_symbols = TRUE,
                          calibration = if (zhao) -1 else NA,
                          collapse = "mean", cluster = FALSE,
                          id_col = "sample_id") {
  resolved <- .resolve_input(object, samples, layer, id_col)
  expr <- resolved$expr
  samples <- resolved$samples
  pam50_expr <- if (map_symbols) map_to_pam50(expr) else expr

  used <- pam50_expr
  if (zhao) {
    if (is.null(samples) || is.null(subtype_col)) {
      stop("zhao = TRUE requires `samples` and `subtype_col`.", call. = FALSE)
    }
    if (!id_col %in% colnames(samples) || !subtype_col %in% colnames(samples)) {
      stop(sprintf("`samples` must contain columns '%s' and '%s'.", id_col, subtype_col),
           call. = FALSE)
    }
    subt <- stats::setNames(as.character(samples[[subtype_col]]), samples[[id_col]])
    subt <- subt[colnames(pam50_expr)]
    used <- zhao_center(pam50_expr, subt)
  }

  core <- .classify_core(used, calibration = calibration, collapse = collapse)

  res <- list(
    calls = core$scores,
    classify = core,
    params = list(zhao = zhao, calibration = calibration, collapse = collapse,
                  cluster = cluster)
  )

  if (isTRUE(cluster)) {
    intr <- intrinsic_genes()
    int_mask <- rownames(expr) %in% intr$Name
    int_expr <- expr[int_mask, , drop = FALSE]
    int_expr <- int_expr[!duplicated(rownames(int_expr)), , drop = FALSE]
    int_medctr <- median_center(int_expr)
    cl <- centroid_clustering(int_medctr, write_files = FALSE)
    clow <- identify_claudin_low(cl, int_medctr)
    res$intrinsic <- list(int_expr = int_expr, int_medctr = int_medctr,
                          clusters = cl, claudin_low = clow)
    res$calls$Call_CLOW <- ifelse(rownames(res$calls) %in% clow, "CLOW", res$calls$Call)

    # silhouette widths (clusters = CLOW-adjusted calls; 1 - Pearson distance)
    if (requireNamespace("cluster", quietly = TRUE)) {
      ids <- rownames(res$calls)
      cl_idx <- as.numeric(as.factor(res$calls$Call_CLOW))
      names(cl_idx) <- ids
      d <- stats::as.dist(1 - stats::cor(int_expr[, ids, drop = FALSE], method = "pearson"))
      sil <- cluster::silhouette(cl_idx, d)
      res$calls$silhouette_width <- sil[, 3]
      res$intrinsic$silhouette <- sil
    }
  }

  structure(res, class = "pam50_result")
}

#' The pam50_result object
#'
#' Returned by [pam50_subtype()]. A list with `$calls` (the scores/calls data
#' frame), `$classify` (prediction intermediates), `$params`, and optionally
#' `$intrinsic`. Has `print`, `summary`, `plot`, and `as.data.frame` methods.
#'
#' @param x,object A `pam50_result`.
#' @param ... Passed on (e.g. to [plot_pam50_scores()]).
#' @name pam50_result
NULL

#' @rdname pam50_result
#' @export
print.pam50_result <- function(x, ...) {
  n <- nrow(x$calls)
  cat(sprintf("<pam50_result>  %d samples\n", n))
  call_col <- if (!is.null(x$calls$Call_CLOW)) "Call_CLOW" else "Call"
  tab <- table(x$calls[[call_col]])
  cat(sprintf("  %s: %s\n", call_col,
              paste(names(tab), tab, sep = "=", collapse = "  ")))
  if (!is.null(x$intrinsic)) {
    cat(sprintf("  claudin-low: %d\n", length(x$intrinsic$claudin_low)))
  }
  cat("  $calls = table; summary() for ROR groups; plot() for score boxplots\n")
  invisible(x)
}

#' @rdname pam50_result
#' @export
summary.pam50_result <- function(object, ...) {
  calls <- object$calls
  structure(list(
    n = nrow(calls),
    subtype = table(calls$Call),
    ror_s_group = table(calls[["ROR-S Group (Subtype Only)"]]),
    ror_p_group = table(calls[["ROR-P Group (Subtype + Proliferation)"]])
  ), class = "summary.pam50_result")
}

#' @export
print.summary.pam50_result <- function(x, ...) {
  cat(sprintf("pam50_result summary (%d samples)\n\nSubtype:\n", x$n))
  print(x$subtype)
  cat("\nROR-S group:\n"); print(x$ror_s_group)
  cat("\nROR-P group:\n"); print(x$ror_p_group)
  invisible(x)
}

#' @rdname pam50_result
#' @export
as.data.frame.pam50_result <- function(x, ...) x$calls

#' @rdname pam50_result
#' @export
plot.pam50_result <- function(x, ...) plot_pam50_scores(x, ...)

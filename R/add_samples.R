# Incremental calling: add new samples to an already-called cohort without
# shifting the existing calls. PAM50 normalisation (Zhao4 / median-centring) is
# cohort-dependent, so simply re-calling the whole enlarged cohort can move
# borderline calls for prior samples. The Perou Lab recommendation is to recall
# ONE new sample at a time against the fixed reference cohort and keep only that
# sample's call. pam50_add_samples() automates that loop.

#' Add new samples to an already-called cohort (one-at-a-time recall)
#'
#' Calls each new sample against the reference cohort plus *only that sample*,
#' keeping the reference cohort's existing calls fixed. This follows the Perou Lab
#' recommendation: because PAM50 normalisation depends on the cohort "denominator",
#' adding samples in a single batch can shift prior borderline calls, so each new
#' sample is recalled individually against the stable reference.
#'
#' Scope: this returns PAM50 subtype calls + scores. Intrinsic clustering /
#' claudin-low are whole-cohort operations and are **not** done incrementally here;
#' recompute them on the final cohort with [pam50_subtype()] if needed.
#'
#' @param reference_expr Reference expression matrix (genes x samples), columns =
#'   the already-called samples.
#' @param reference_samples Sample sheet for the reference cohort (an `id_col` and,
#'   for Zhao, a `subtype_col`).
#' @param new_expr Expression matrix of the new samples (genes x samples).
#' @param new_samples Sample sheet for the new samples (same columns as
#'   `reference_samples`).
#' @param reference_calls Optional pre-computed calls for the reference cohort (a
#'   data frame as returned in `pam50_subtype()$calls`). If `NULL`, computed once
#'   with the same settings.
#' @param subtype_col,zhao,calibration,collapse,id_col,map_symbols Passed to
#'   [pam50_subtype()]; the defaults match it (Zhao on when `subtype_col` is given).
#' @return A data frame of calls for the full cohort: the reference calls
#'   (unchanged) followed by one row per new sample, with a logical `added` column
#'   (`FALSE` = reference, `TRUE` = newly added).
#' @seealso [pam50_subtype()]
#' @export
pam50_add_samples <- function(reference_expr, reference_samples,
                              new_expr, new_samples,
                              reference_calls = NULL, subtype_col = NULL,
                              zhao = !is.null(subtype_col),
                              calibration = if (zhao) -1 else NA,
                              collapse = "mean", id_col = "sample_id",
                              map_symbols = TRUE) {
  reference_expr <- as.matrix(reference_expr)
  new_expr <- as.matrix(new_expr)

  # Expression columns are the source of truth; align the sheets to them.
  ref_ids <- colnames(reference_expr)
  new_ids <- colnames(new_expr)
  if (any(new_ids %in% ref_ids)) {
    stop("new sample id(s) collide with reference sample ids: ",
         paste(intersect(new_ids, ref_ids), collapse = ", "), call. = FALSE)
  }
  rownames(reference_samples) <- reference_samples[[id_col]]
  rownames(new_samples) <- new_samples[[id_col]]
  reference_samples <- reference_samples[ref_ids, , drop = FALSE]
  new_samples <- new_samples[new_ids, , drop = FALSE]

  # Common genes + sample-sheet columns.
  genes <- intersect(rownames(reference_expr), rownames(new_expr))
  if (length(genes) == 0L) stop("reference_expr and new_expr share no genes.", call. = FALSE)
  reference_expr <- reference_expr[genes, , drop = FALSE]
  new_expr <- new_expr[genes, , drop = FALSE]
  cols <- intersect(colnames(reference_samples), colnames(new_samples))

  sub_args <- list(subtype_col = subtype_col, zhao = zhao, calibration = calibration,
                   collapse = collapse, cluster = FALSE, id_col = id_col,
                   map_symbols = map_symbols)

  # Reference calls (kept fixed). Compute once if not supplied.
  if (is.null(reference_calls)) {
    ref_res <- do.call(pam50_subtype,
                       c(list(reference_expr, samples = reference_samples), sub_args))
    reference_calls <- ref_res$calls
  }

  # One new sample at a time against the fixed reference.
  new_rows <- vector("list", length(new_ids))
  for (j in seq_along(new_ids)) {
    nid <- new_ids[j]
    aug_expr <- cbind(reference_expr, new_expr[, nid, drop = FALSE])
    aug_samples <- rbind(reference_samples[, cols, drop = FALSE],
                         new_samples[nid, cols, drop = FALSE])
    res_j <- do.call(pam50_subtype, c(list(aug_expr, samples = aug_samples), sub_args))
    new_rows[[j]] <- res_j$calls[nid, , drop = FALSE]
    message(sprintf("added %s -> %s", nid, res_j$calls[nid, "Call"]))
  }
  new_calls <- do.call(rbind, new_rows)

  # Combine on shared columns, flag provenance.
  shared <- intersect(colnames(reference_calls), colnames(new_calls))
  out <- rbind(reference_calls[, shared, drop = FALSE], new_calls[, shared, drop = FALSE])
  out$added <- c(rep(FALSE, nrow(reference_calls)), rep(TRUE, nrow(new_calls)))
  out
}

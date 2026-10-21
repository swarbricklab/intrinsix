# Input normalisation. intrinsix accepts the expression as a plain matrix, a
# SummarizedExperiment / SingleCellExperiment, an AnnData object (anndataR), or a
# path to an .h5ad file. .resolve_input() turns any of these into a plain
# (expression matrix [genes x samples], sample-sheet data.frame) pair. The
# SummarizedExperiment and anndataR packages are Suggests, loaded only on demand.

# @noRd
.resolve_input <- function(object, samples = NULL, layer = NULL, id_col = "sample_id",
                           gene_col = NULL, pam50_id_col = NULL) {
  # .h5ad path -> anndataR AnnData
  if (is.character(object) && length(object) == 1L && grepl("\\.h5ad$", object, ignore.case = TRUE)) {
    if (!requireNamespace("anndataR", quietly = TRUE)) {
      stop("Reading an .h5ad file needs the 'anndataR' package (R >= 4.5).", call. = FALSE)
    }
    object <- anndataR::read_h5ad(object)
  }

  # anndataR AnnData (layers/X are obs x var -> transpose to genes x samples)
  if (inherits(object, "AbstractAnnData")) {
    mat <- if (!is.null(layer)) object$layers[[layer]] else object$X
    if (is.null(mat)) {
      stop(sprintf("AnnData has no %s.",
                   if (is.null(layer)) "X matrix" else sprintf("layer '%s'", layer)), call. = FALSE)
    }
    expr <- t(as.matrix(mat))
    if (!is.null(gene_col)) {
      if (!gene_col %in% colnames(object$var)) {
        stop(sprintf("gene_col '%s' not in var.", gene_col), call. = FALSE)
      }
      gnames <- as.character(object$var[[gene_col]])
    } else {
      gnames <- rownames(object$var); if (is.null(gnames)) gnames <- object$var_names
    }
    snames <- rownames(object$obs); if (is.null(snames)) snames <- object$obs_names
    rownames(expr) <- gnames; colnames(expr) <- snames
    pam50_ids <- NULL
    if (!is.null(pam50_id_col)) {
      if (!pam50_id_col %in% colnames(object$var)) {
        stop(sprintf("pam50_id_col '%s' not in var.", pam50_id_col), call. = FALSE)
      }
      pam50_ids <- as.character(object$var[[pam50_id_col]])  # positional, = expr rows
    }
    if (is.null(samples)) {
      samples <- as.data.frame(object$obs)
      if (!id_col %in% colnames(samples)) samples[[id_col]] <- snames
    }
    return(list(expr = expr, samples = samples, pam50_ids = pam50_ids))
  }

  # SummarizedExperiment / SingleCellExperiment (assays are genes x samples)
  if (inherits(object, "SummarizedExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("A SummarizedExperiment input needs the 'SummarizedExperiment' package.", call. = FALSE)
    }
    expr <- as.matrix(SummarizedExperiment::assay(object, if (is.null(layer)) 1L else layer))
    cd <- as.data.frame(SummarizedExperiment::colData(object))
    rd <- SummarizedExperiment::rowData(object)
    if (!is.null(gene_col)) {
      if (!gene_col %in% colnames(rd)) stop(sprintf("gene_col '%s' not in rowData.", gene_col), call. = FALSE)
      rownames(expr) <- as.character(rd[[gene_col]])
    } else if (is.null(rownames(expr))) {
      rownames(expr) <- rownames(object)
    }
    pam50_ids <- NULL
    if (!is.null(pam50_id_col)) {
      if (!pam50_id_col %in% colnames(rd)) {
        stop(sprintf("pam50_id_col '%s' not in rowData.", pam50_id_col), call. = FALSE)
      }
      pam50_ids <- as.character(rd[[pam50_id_col]])  # positional, = expr rows
    }
    if (is.null(colnames(expr))) colnames(expr) <- rownames(cd)
    if (is.null(samples)) {
      samples <- cd
      if (!id_col %in% colnames(samples)) samples[[id_col]] <- rownames(cd)
    }
    return(list(expr = expr, samples = samples, pam50_ids = pam50_ids))
  }

  # plain matrix / data.frame (samples supplied separately)
  list(expr = as.matrix(object), samples = samples, pam50_ids = NULL)
}

# Config-driven whole-pipeline runner. Collapses what the (now retired) Snakemake
# rules did into one call: read config -> filter/dedup samples -> pam50_subtype()
# -> write the output files. Reachable as pam50_run() or `intrinsix run --config ...`,
# and intended to be the command a DVC stage invokes (in a conda env or container).

`%||%` <- function(a, b) if (is.null(a)) b else a

# @noRd
.validate_run_config <- function(cfg) {
  need <- function(path, x) if (is.null(x)) stop(sprintf("config: missing '%s'", path), call. = FALSE)
  need("inputs$expression", cfg$inputs$expression)
  need("inputs$samples", cfg$inputs$samples)
  need("inputs$sample_id_column", cfg$inputs$sample_id_column)
  need("outs$output_dir", cfg$outs$output_dir)
  invisible(TRUE)
}

# @noRd
.apply_filter <- function(samples, cfg) {
  filt <- cfg$filter
  if (is.null(filt$column) || is.na(filt$column) || !nzchar(filt$column)) return(samples)
  if (!filt$column %in% colnames(samples)) stop(sprintf("filter.column '%s' not in sample sheet", filt$column), call. = FALSE)
  include <- unlist(filt$include)
  samples[as.character(samples[[filt$column]]) %in% as.character(include), , drop = FALSE]
}

# @noRd
.apply_dedup <- function(samples, cfg, id_col) {
  dd <- cfg$deduplicate
  if (is.null(dd) || is.null(dd$by) || is.na(dd$by) || !nzchar(dd$by)) return(samples)
  if (!dd$by %in% colnames(samples)) stop(sprintf("deduplicate.by '%s' not in sample sheet", dd$by), call. = FALSE)
  sel <- dd$selections
  grp <- as.character(samples[[dd$by]])
  keep <- character(0)
  for (g in unique(grp)) {
    members <- rownames(samples)[grp == g]
    pref <- if (!is.null(sel)) sel[[g]] else NULL
    keep <- c(keep, if (!is.null(pref) && pref %in% members) pref else sort(members)[1])
  }
  samples[keep, , drop = FALSE]
}

#' Run the PAM50 pipeline from a config file
#'
#' End-to-end, config-driven run that replaces the Snakemake workflow: reads the
#' expression matrix + sample sheet, applies the (optional) sample filter and
#' replicate deduplication, runs [pam50_subtype()] (map -> optional Zhao ->
#' classify -> optional intrinsic clustering), and writes the outputs. This is the
#' command a DVC stage invokes (via `intrinsix run --config config.yaml`), in a conda
#' environment or the container.
#'
#' @param config Path to a YAML config, or an already-parsed config list. See
#'   the bundled template: `system.file("extdata", "config.example.yaml",
#'   package = "intrinsix")`.
#' @return The [pam50_result], invisibly. Side effect: writes to
#'   `config$outs$output_dir`.
#' @export
pam50_run <- function(config) {
  cfg <- if (is.character(config)) yaml::read_yaml(config) else config
  .validate_run_config(cfg)

  id_col <- cfg$inputs$sample_id_column
  expr_in <- cfg$inputs$expression
  samples_in <- if (!is.null(cfg$inputs$samples)) {
    utils::read.table(cfg$inputs$samples, sep = "\t", header = TRUE,
                      check.names = FALSE, stringsAsFactors = FALSE)
  } else NULL
  # .h5ad -> pass the path; SE/AnnData -> rare in a config; else read the TSV matrix.
  object <- if (grepl("\\.h5ad$", expr_in, ignore.case = TRUE)) {
    expr_in
  } else {
    as.matrix(utils::read.table(expr_in, sep = "\t", header = TRUE,
                                row.names = 1, check.names = FALSE))
  }
  resolved <- .resolve_input(object, samples_in, cfg$inputs$layer, id_col)
  expr <- resolved$expr
  samples <- resolved$samples
  if (is.null(samples)) {
    stop("No sample sheet: set inputs$samples, or use an annotated input (SE/h5ad).", call. = FALSE)
  }
  rownames(samples) <- samples[[id_col]]

  samples <- .apply_filter(samples, cfg)
  samples <- .apply_dedup(samples, cfg, id_col)
  keep <- intersect(rownames(samples), colnames(expr))
  samples <- samples[keep, , drop = FALSE]
  message(sprintf("pam50_run: %d samples after filter/dedup", length(keep)))

  cluster <- isTRUE(cfg$intrinsic_clustering$enabled)
  res <- pam50_subtype(
    expr[, keep, drop = FALSE],
    samples = samples, subtype_col = cfg$zhao_normalization$subtype_column,
    zhao = isTRUE(cfg$zhao_normalization$enabled),
    calibration = cfg$classifier$calibration %||% -1,
    collapse = cfg$classifier$collapse_method %||% "mean",
    cluster = cluster, id_col = id_col, map_symbols = TRUE
  )

  out <- cfg$outs$output_dir
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  prefix <- cfg$classifier$prefix %||% "PAM50scores"

  # raw scores (same layout as the reference *_pam50scores.txt)
  score_cols <- setdiff(colnames(res$calls), c("Call_CLOW", "silhouette_width"))
  utils::write.table(res$calls[, score_cols], file.path(out, paste0(prefix, "_pam50scores.txt")),
                     sep = "\t", col.names = NA, quote = TRUE)

  # enriched call table (scores + sample-sheet metadata)
  meta_cols <- setdiff(colnames(samples), id_col)
  calls <- data.frame(sample = rownames(res$calls),
                      samples[rownames(res$calls), meta_cols, drop = FALSE],
                      res$calls, check.names = FALSE, row.names = NULL)
  colnames(calls)[1] <- id_col
  utils::write.table(calls, file.path(out, "pam50_calls.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)

  if (cluster && !is.null(res$intrinsic)) {
    hdir <- file.path(out, "intrinsic_heatmap")
    dir.create(hdir, recursive = TRUE, showWarnings = FALSE)

    clow_df <- data.frame(sample = rownames(res$calls),
                          PAM50_Call = res$calls$Call,
                          PAM50_Call_CLOW = res$calls$Call_CLOW,
                          check.names = FALSE, row.names = NULL)
    if (!is.null(res$calls$silhouette_width)) {
      clow_df$silhouette_width <- res$calls$silhouette_width
      clow_df$sil_width_bin <- ifelse(res$calls$silhouette_width < 0, "Negative", "Positive")
    }
    utils::write.table(clow_df, file.path(out, "pam50_calls_clow.tsv"),
                       sep = "\t", row.names = FALSE, quote = FALSE)

    if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
      ann <- clow_df[, intersect(c("PAM50_Call_CLOW", "sil_width_bin"), colnames(clow_df)), drop = FALSE]
      rownames(ann) <- clow_df$sample
      ids <- colnames(res$intrinsic$int_medctr)
      pam50_heatmap(res$intrinsic$int_medctr, annotations = ann[ids, , drop = FALSE],
                    title = sprintf("Intrinsic subtype heatmap (n=%d)", length(ids)),
                    file = file.path(hdir, "intrinsic_medctr_heatmap.pdf"))
    }
    if (!is.null(res$intrinsic$silhouette)) {
      grDevices::pdf(file.path(hdir, "silhouette_plot.pdf"), width = 9, height = 8)
      graphics::plot(res$intrinsic$silhouette, border = NA)
      graphics::abline(v = 0)
      grDevices::dev.off()
    }
  }

  message(sprintf("pam50_run: wrote outputs to %s", out))
  invisible(res)
}

# Config-driven whole-pipeline runner. Collapses what the (now retired) Snakemake
# rules did into one call: read config -> filter/dedup samples -> pam50_subtype()
# -> write the output files. Reachable as pam50_run() or `intrinsix run --config ...`,
# and intended to be the command a DVC stage invokes (in a conda env or container).

`%||%` <- function(a, b) if (is.null(a)) b else a

# @noRd
.validate_run_config <- function(cfg) {
  need <- function(path, x) if (is.null(x)) stop(sprintf("config: missing '%s'", path), call. = FALSE)
  need("inputs$expression", cfg$inputs$expression)
  # inputs$samples is optional: an h5ad / SummarizedExperiment carries its own obs/colData.
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
#' @section Output:
#' Writes a single well-formed call table, `pam50_calls.tsv`, to
#' `outs$output_dir`: one row per sample, a named id column
#' (`inputs$sample_id_column`), then the carried-through sample metadata, then the
#' classifier columns under clean, unique, snake_case names — `call`,
#' `confidence`, `score_basal`..`score_normal`, `ror_s`, `ror_s_group`,
#' `proliferation_score`, `ror_p`, `ror_p_group`, `er_signature`,
#' `her2_signature`, plus `call_clow`, `silhouette_width`, `sil_width_bin` when
#' `intrinsic_clustering` is enabled. Values are unquoted. (The in-memory
#' [pam50_result] keeps the raw Perou-native column names; only the file output is
#' renamed.)
#'
#' Sample metadata is inherited from the input's `obs` / `colData` (or the
#' `inputs$samples` sheet). Every metadata column is carried through by default;
#' set `inputs$metadata_columns` to an allowlist to keep only those. If a carried
#' metadata column name collides with a classifier output column, the run stops
#' with an error (rename it upstream, or drop it via `inputs$metadata_columns`).
#'
#' Set `classifier$legacy_perou_scores: true` to additionally emit the legacy
#' Perou-format `"<prefix>_pam50scores.txt"` (quoted cells, unnamed first column,
#' duplicate `Her2`) for downstream Perou Lab scripts that expect it. Off by
#' default.
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
  resolved <- .resolve_input(object, samples_in, cfg$inputs$layer, id_col,
                             gene_col = cfg$inputs$gene_col,
                             pam50_id_col = cfg$inputs$pam50_id_col,
                             entrez_col = cfg$inputs$entrez_col)
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
    cluster = cluster, id_col = id_col, map_symbols = TRUE,
    pam50_ids = resolved$pam50_ids   # drift-proof IDs when inputs$pam50_id_col set
  )

  out <- cfg$outs$output_dir
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  prefix <- cfg$classifier$prefix %||% "PAM50scores"

  # Optional legacy Perou-format scores file (off by default). Reproduces the
  # historical "<prefix>_pam50scores.txt" exactly: native classifier column names
  # (which include a duplicate 'Her2' and the ER/HER2 signature columns), quoted
  # cells, unnamed first column — for Perou Lab scripts that parse that format.
  if (isTRUE(cfg$classifier$legacy_perou_scores)) {
    score_cols <- setdiff(colnames(res$calls), c("Call_CLOW", "silhouette_width"))
    utils::write.table(res$calls[, score_cols], file.path(out, paste0(prefix, "_pam50scores.txt")),
                       sep = "\t", col.names = NA, quote = TRUE)
  }

  # Clean, well-formed call table (the default, and only, tabular output). The
  # classifier core keeps the Perou-native column names; we rename them to a
  # unique, snake_case schema here (res$calls itself is left untouched). The first
  # `n_core` columns are the fixed classifier block (see .classify_core); any
  # trailing columns are the optional clustering annotations.
  expected_native <- c(
    "Basal", "Her2", "LumA", "LumB", "Normal", "Call", "Confidence",
    "ROR-S (Subtype Only)", "ROR-S Group (Subtype Only)", "Proliferation Score",
    "ROR-P (Subtype + Proliferation)", "ROR-P Group (Subtype + Proliferation)",
    "ER", "Her2"
  )
  clean_names <- c(
    "score_basal", "score_her2", "score_luma", "score_lumb", "score_normal",
    "call", "confidence", "ror_s", "ror_s_group", "proliferation_score",
    "ror_p", "ror_p_group", "er_signature", "her2_signature"
  )
  n_core <- length(expected_native)
  if (ncol(res$calls) < n_core ||
      !identical(colnames(res$calls)[seq_len(n_core)], expected_native)) {
    stop("intrinsix: unexpected classifier column layout; cannot build the clean ",
         "output. This is an internal error — please report it.", call. = FALSE)
  }
  scores_clean <- res$calls[, seq_len(n_core), drop = FALSE]
  colnames(scores_clean) <- clean_names
  if (!is.null(res$calls$Call_CLOW)) scores_clean$call_clow <- res$calls$Call_CLOW
  if (!is.null(res$calls$silhouette_width)) {
    scores_clean$silhouette_width <- res$calls$silhouette_width
    scores_clean$sil_width_bin <- ifelse(res$calls$silhouette_width < 0, "Negative", "Positive")
  }

  # Sample metadata inherited from obs/colData (or inputs$samples). Carry every
  # column by default; inputs$metadata_columns narrows it to an allowlist.
  meta_cols <- setdiff(colnames(samples), id_col)
  allow <- cfg$inputs$metadata_columns
  if (!is.null(allow)) {
    allow <- as.character(unlist(allow))
    absent <- setdiff(allow, colnames(samples))
    if (length(absent)) {
      stop(sprintf("inputs$metadata_columns: column(s) not in sample metadata: %s",
                   paste(absent, collapse = ", ")), call. = FALSE)
    }
    meta_cols <- setdiff(allow, id_col)
  }
  # Error on clash: a carried metadata column must not shadow a classifier column.
  clash <- intersect(meta_cols, colnames(scores_clean))
  if (length(clash)) {
    stop(sprintf(paste0("intrinsix: sample-metadata column(s) collide with classifier ",
                        "output columns: %s. Rename them upstream, or exclude them via ",
                        "inputs$metadata_columns."),
                 paste(clash, collapse = ", ")), call. = FALSE)
  }

  ids <- rownames(res$calls)
  calls <- data.frame(.id = ids,
                      samples[ids, meta_cols, drop = FALSE],
                      scores_clean,
                      check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE)
  colnames(calls)[1] <- id_col
  utils::write.table(calls, file.path(out, "pam50_calls.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)

  if (cluster && !is.null(res$intrinsic)) {
    hdir <- file.path(out, "intrinsic_heatmap")
    dir.create(hdir, recursive = TRUE, showWarnings = FALSE)

    if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
      ann_cols <- intersect(c("call_clow", "sil_width_bin"), colnames(scores_clean))
      ann <- scores_clean[, ann_cols, drop = FALSE]
      rownames(ann) <- ids
      hm_ids <- colnames(res$intrinsic$int_medctr)
      pam50_heatmap(res$intrinsic$int_medctr, annotations = ann[hm_ids, , drop = FALSE],
                    title = sprintf("Intrinsic subtype heatmap (n=%d)", length(hm_ids)),
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

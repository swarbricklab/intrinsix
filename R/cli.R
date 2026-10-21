# Command-line interface. Exposes the same pam50 functions that the Snakemake
# rule wrappers call, so each operation is reachable from both a rule and the shell
# (the spotty "one function, two front-ends" idea, in R form). Dependency-free:
# a minimal base-R flag parser, no optparse.
#
#   pam50 classify --expr E.tsv --out S.tsv [--calibration -1] [--collapse mean] [--map-symbols]
#   pam50 zhao     --expr E.tsv --samples M.tsv --subtype-col clinical_subtype [--id-col sample_id] --out Z.tsv
#   pam50 cluster  --expr-medctr M.tsv --outdir DIR [--out-claudin C.txt]

.pam50_parse <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (startsWith(a, "--")) {
      key <- sub("^--", "", a)
      if (i + 1L <= length(args) && !startsWith(args[[i + 1L]], "--")) {
        out[[key]] <- args[[i + 1L]]; i <- i + 2L
      } else {
        out[[key]] <- TRUE; i <- i + 1L
      }
    } else {
      i <- i + 1L
    }
  }
  out
}

.read_matrix <- function(path, sep = "\t") {
  as.matrix(utils::read.table(path, sep = sep, header = TRUE, row.names = 1, check.names = FALSE))
}

.require_opts <- function(opt, needed, cmd) {
  miss <- needed[!needed %in% names(opt)]
  if (length(miss)) {
    stop(sprintf("pam50 %s: missing required option(s): %s",
                 cmd, paste0("--", miss, collapse = ", ")), call. = FALSE)
  }
}

.cli_classify <- function(args) {
  opt <- .pam50_parse(args)
  .require_opts(opt, c("expr", "out"), "classify")
  expr <- .read_matrix(opt$expr)
  if (isTRUE(opt$`map-symbols`)) expr <- map_to_pam50(expr)
  calibration <- if (!is.null(opt$calibration)) as.integer(opt$calibration) else -1L
  collapse <- if (!is.null(opt$collapse)) opt$collapse else "mean"
  res <- predict_pam50(expr, calibration = calibration, collapse = collapse)
  utils::write.table(res, opt$out, sep = "\t", col.names = NA, quote = TRUE)
  message(sprintf("classify: wrote %d calls -> %s", nrow(res), opt$out))
}

.cli_zhao <- function(args) {
  opt <- .pam50_parse(args)
  .require_opts(opt, c("expr", "samples", "subtype-col", "out"), "zhao")
  expr <- .read_matrix(opt$expr)
  id_col <- if (!is.null(opt$`id-col`)) opt$`id-col` else "sample_id"
  ss <- utils::read.table(opt$samples, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  subt <- stats::setNames(as.character(ss[[opt$`subtype-col`]]), ss[[id_col]])
  subt <- subt[colnames(expr)]
  z <- zhao_center(expr, subt)
  utils::write.table(data.frame(Gene = rownames(z), z, check.names = FALSE),
                     opt$out, sep = "\t", row.names = FALSE, quote = FALSE)
  message(sprintf("zhao: wrote %d x %d -> %s", nrow(z), ncol(z), opt$out))
}

.cli_run <- function(args) {
  opt <- .pam50_parse(args)
  .require_opts(opt, "config", "run")
  pam50_run(opt$config)
}

.cli_add <- function(args) {
  opt <- .pam50_parse(args)
  .require_opts(opt, "config", "add")
  cfg <- yaml::read_yaml(opt$config)
  id_col <- cfg$sample_id_column %||% "sample_id"
  # expression: pass an .h5ad path through (resolved downstream), else read the TSV.
  read_expr <- function(p) if (grepl("\\.h5ad$", p, ignore.case = TRUE)) p else .read_matrix(p)
  read_sheet <- function(p) if (is.null(p)) NULL else {
    utils::read.table(p, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  }
  ref_calls <- NULL
  if (!is.null(cfg$reference$calls)) {
    ref_calls <- utils::read.table(cfg$reference$calls, sep = "\t", header = TRUE,
                                   row.names = 1, check.names = FALSE)
  }
  out <- pam50_add_samples(
    read_expr(cfg$reference$expression), read_sheet(cfg$reference$samples),
    read_expr(cfg$new$expression), read_sheet(cfg$new$samples),
    reference_calls = ref_calls, subtype_col = cfg$subtype_column,
    zhao = cfg$zhao %||% !is.null(cfg$subtype_column),
    calibration = cfg$classifier$calibration %||% -1,
    collapse = cfg$classifier$collapse_method %||% "mean",
    id_col = id_col, map_symbols = cfg$map_symbols %||% TRUE,
    claudin_low = isTRUE(cfg$claudin_low), layer = cfg$layer
  )
  outpath <- cfg$outs$output
  dir.create(dirname(outpath), recursive = TRUE, showWarnings = FALSE)
  df <- data.frame(id = rownames(out), out, check.names = FALSE, row.names = NULL)
  colnames(df)[1] <- id_col
  utils::write.table(df, outpath, sep = "\t", row.names = FALSE, quote = FALSE)
  message(sprintf("pam50 add: %d calls (%d new) -> %s", nrow(out), sum(out$added), outpath))
}

.cli_cluster <- function(args) {
  opt <- .pam50_parse(args)
  .require_opts(opt, c("expr-medctr", "outdir"), "cluster")
  m <- .read_matrix(opt$`expr-medctr`)
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
  res <- centroid_clustering(m, prefix = "intrinsic_expr_medCtr", path = opt$outdir, write_files = TRUE)
  clow <- identify_claudin_low(res, m)
  if (!is.null(opt$`out-claudin`)) writeLines(clow, opt$`out-claudin`)
  message(sprintf("cluster: %d claudin-low sample(s)%s", length(clow),
                  if (!is.null(opt$`out-claudin`)) paste0(" -> ", opt$`out-claudin`) else ""))
}

.pam50_usage <- function() {
  message("pam50 <command> [options]\n",
          "\nCommands:\n",
          "  run       --config config.yaml   (whole pipeline, config-driven)\n",
          "  add       --config add.yaml      (incremental: new samples vs a fixed reference)\n",
          "  classify  --expr E.tsv --out S.tsv [--calibration -1] [--collapse mean] [--map-symbols]\n",
          "  zhao      --expr E.tsv --samples M.tsv --subtype-col COL [--id-col sample_id] --out Z.tsv\n",
          "  cluster   --expr-medctr M.tsv --outdir DIR [--out-claudin C.txt]\n")
}

#' pam50 command-line entry point
#'
#' Dispatches `run` / `add` / `classify` / `zhao` / `cluster` subcommands to the
#' corresponding package functions. `run` executes the whole pipeline from a YAML
#' config ([pam50_run()]); `add` does incremental calling of new samples against a
#' fixed reference ([pam50_add_samples()]). Both are the commands DVC stages invoke.
#' Also callable directly: `pam50::cli(c("run", "--config", "config.yaml"))`.
#'
#' @param args Character vector of command-line arguments (default: the process
#'   arguments after `--args`).
#' @return Invisibly `0` on success; stops on error.
#' @export
cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) < 1L || args[[1]] %in% c("-h", "--help", "help")) {
    .pam50_usage(); return(invisible(0L))
  }
  cmd <- args[[1]]; rest <- args[-1]
  switch(cmd,
    run = .cli_run(rest),
    add = .cli_add(rest),
    classify = .cli_classify(rest),
    zhao = .cli_zhao(rest),
    cluster = .cli_cluster(rest),
    { .pam50_usage(); stop(sprintf("unknown command: %s", cmd), call. = FALSE) }
  )
  invisible(0L)
}

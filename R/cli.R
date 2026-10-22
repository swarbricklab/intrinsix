# Command-line interface. Exposes the same intrinsix functions used from R, so each
# operation is reachable from both a DVC stage / the shell and an R session (the
# spotty "one function, two front-ends" idea, in R form). Argument parsing uses the
# lightweight pure-R 'argparser' package: a git-style subcommand dispatch
# (run/add/classify/zhao/cluster), each subcommand with its own parser so it gets
# auto `--help`, type coercion, and typo/missing-argument rejection for free.
#
#   intrinsix run       --config config.yaml
#   intrinsix add       --config add.yaml
#   intrinsix classify  --expr E.tsv --out S.tsv [--calibration -1] [--collapse mean] [--map-symbols]
#   intrinsix zhao      --expr E.tsv --samples M.tsv --subtype-col clinical_subtype [--id-col sample_id] --out Z.tsv
#   intrinsix cluster   --expr-medctr M.tsv --outdir DIR [--out-claudin C.txt]

.read_matrix <- function(path, sep = "\t") {
  as.matrix(utils::read.table(path, sep = sep, header = TRUE, row.names = 1, check.names = FALSE))
}

# argparser leaves optional arguments without a default as NA; enforce the ones we
# require. `needed` are flag names without the leading "--" (argparser stores them
# with "-" replaced by "_").
.need <- function(argv, needed, cmd) {
  miss <- needed[vapply(gsub("-", "_", needed), function(k) {
    v <- argv[[k]]; is.null(v) || (length(v) == 1L && is.na(v))
  }, logical(1))]
  if (length(miss)) {
    stop(sprintf("intrinsix %s: missing required --%s", cmd,
                 paste(miss, collapse = ", --")), call. = FALSE)
  }
}

.cli_classify <- function(args) {
  p <- argparser::arg_parser("PAM50 classify an expression matrix", name = "intrinsix classify")
  p <- argparser::add_argument(p, "--expr", type = "character",
                               help = "expression matrix TSV (pcrID rows, or HGNC with --map-symbols)")
  p <- argparser::add_argument(p, "--out", type = "character", help = "output scores TSV")
  p <- argparser::add_argument(p, "--calibration", type = "integer", default = -1L,
                               help = "calibration parameter")
  p <- argparser::add_argument(p, "--collapse", default = "mean",
                               help = "duplicate-gene collapse method: mean|iqr")
  p <- argparser::add_argument(p, "--map-symbols", flag = TRUE,
                               help = "map HGNC symbols to pcrID and subset to the 50 PAM50 genes")
  argv <- argparser::parse_args(p, args)
  .need(argv, c("expr", "out"), "classify")
  expr <- .read_matrix(argv$expr)
  if (isTRUE(argv$map_symbols)) expr <- map_to_pam50(expr)
  res <- predict_pam50(expr, calibration = argv$calibration, collapse = argv$collapse)
  utils::write.table(res, argv$out, sep = "\t", col.names = NA, quote = TRUE)
  message(sprintf("classify: wrote %d calls -> %s", nrow(res), argv$out))
}

.cli_zhao <- function(args) {
  p <- argparser::arg_parser("Zhao4 subtype-specific quantile centering", name = "intrinsix zhao")
  p <- argparser::add_argument(p, "--expr", type = "character", help = "PAM50 expression matrix TSV")
  p <- argparser::add_argument(p, "--samples", type = "character", help = "sample sheet TSV")
  p <- argparser::add_argument(p, "--subtype-col", type = "character",
                               help = "clinical (ER x HER2) subtype column in the sample sheet")
  p <- argparser::add_argument(p, "--id-col", default = "sample_id",
                               help = "sample ID column in the sample sheet")
  p <- argparser::add_argument(p, "--out", type = "character", help = "output normalised matrix TSV")
  argv <- argparser::parse_args(p, args)
  .need(argv, c("expr", "samples", "subtype-col", "out"), "zhao")
  expr <- .read_matrix(argv$expr)
  ss <- utils::read.table(argv$samples, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  subt <- stats::setNames(as.character(ss[[argv$subtype_col]]), ss[[argv$id_col]])
  subt <- subt[colnames(expr)]
  z <- zhao_center(expr, subt)
  utils::write.table(data.frame(Gene = rownames(z), z, check.names = FALSE),
                     argv$out, sep = "\t", row.names = FALSE, quote = FALSE)
  message(sprintf("zhao: wrote %d x %d -> %s", nrow(z), ncol(z), argv$out))
}

.cli_run <- function(args) {
  p <- argparser::arg_parser("Run the whole PAM50 pipeline from a config", name = "intrinsix run")
  p <- argparser::add_argument(p, "--config", type = "character", help = "YAML config path")
  argv <- argparser::parse_args(p, args)
  .need(argv, "config", "run")
  pam50_run(argv$config)
}

.cli_add <- function(args) {
  p <- argparser::arg_parser("Incremental calling: new samples vs a fixed reference",
                             name = "intrinsix add")
  p <- argparser::add_argument(p, "--config", type = "character", help = "YAML config path")
  argv <- argparser::parse_args(p, args)
  .need(argv, "config", "add")
  cfg <- yaml::read_yaml(argv$config)
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
    claudin_low = isTRUE(cfg$claudin_low), layer = cfg$layer, gene_col = cfg$gene_col,
    pam50_id_col = cfg$pam50_id_col
  )
  outpath <- cfg$outs$output
  dir.create(dirname(outpath), recursive = TRUE, showWarnings = FALSE)
  df <- data.frame(id = rownames(out), out, check.names = FALSE, row.names = NULL)
  colnames(df)[1] <- id_col
  utils::write.table(df, outpath, sep = "\t", row.names = FALSE, quote = FALSE)
  message(sprintf("intrinsix add: %d calls (%d new) -> %s", nrow(out), sum(out$added), outpath))
}

.cli_cluster <- function(args) {
  p <- argparser::arg_parser("Intrinsic clustering + claudin-low identification",
                             name = "intrinsix cluster")
  p <- argparser::add_argument(p, "--expr-medctr", type = "character",
                               help = "median-centred intrinsic expression matrix TSV")
  p <- argparser::add_argument(p, "--outdir", type = "character", help = "output directory")
  p <- argparser::add_argument(p, "--out-claudin", type = "character",
                               help = "optional output path for the claudin-low sample list")
  argv <- argparser::parse_args(p, args)
  .need(argv, c("expr-medctr", "outdir"), "cluster")
  m <- .read_matrix(argv$expr_medctr)
  dir.create(argv$outdir, recursive = TRUE, showWarnings = FALSE)
  res <- centroid_clustering(m, prefix = "intrinsic_expr_medCtr", path = argv$outdir, write_files = TRUE)
  clow <- identify_claudin_low(res, m)
  write_clow <- !is.null(argv$out_claudin) && !is.na(argv$out_claudin)
  if (write_clow) writeLines(clow, argv$out_claudin)
  message(sprintf("cluster: %d claudin-low sample(s)%s", length(clow),
                  if (write_clow) paste0(" -> ", argv$out_claudin) else ""))
}

.pam50_usage <- function() {
  message("intrinsix <command> [options]   (run `intrinsix <command> --help` for details)\n",
          "\nCommands:\n",
          "  run       --config config.yaml   (whole pipeline, config-driven)\n",
          "  add       --config add.yaml      (incremental: new samples vs a fixed reference)\n",
          "  classify  --expr E.tsv --out S.tsv [--calibration -1] [--collapse mean] [--map-symbols]\n",
          "  zhao      --expr E.tsv --samples M.tsv --subtype-col COL [--id-col sample_id] --out Z.tsv\n",
          "  cluster   --expr-medctr M.tsv --outdir DIR [--out-claudin C.txt]\n")
}

#' intrinsix command-line entry point
#'
#' Dispatches `run` / `add` / `classify` / `zhao` / `cluster` subcommands to the
#' corresponding package functions. `run` executes the whole pipeline from a YAML
#' config ([pam50_run()]); `add` does incremental calling of new samples against a
#' fixed reference ([pam50_add_samples()]). Both are the commands DVC stages invoke.
#' Each subcommand has its own `--help` (via the `argparser` package). Also callable
#' directly: `intrinsix::cli(c("run", "--config", "config.yaml"))`.
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

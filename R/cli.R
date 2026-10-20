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
          "  classify  --expr E.tsv --out S.tsv [--calibration -1] [--collapse mean] [--map-symbols]\n",
          "  zhao      --expr E.tsv --samples M.tsv --subtype-col COL [--id-col sample_id] --out Z.tsv\n",
          "  cluster   --expr-medctr M.tsv --outdir DIR [--out-claudin C.txt]\n")
}

#' pam50 command-line entry point
#'
#' Dispatches `classify` / `zhao` / `cluster` subcommands to the corresponding
#' package functions. Invoked by the installed `pam50` executable
#' (`inst/exec/pam50`), or directly: `pam50::cli(c("classify", "--expr", ...))`.
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
    classify = .cli_classify(rest),
    zhao = .cli_zhao(rest),
    cluster = .cli_cluster(rest),
    { .pam50_usage(); stop(sprintf("unknown command: %s", cmd), call. = FALSE) }
  )
  invisible(0L)
}

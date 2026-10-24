# Test the config-driven whole-pipeline runner (pam50_run), which replaces the
# former Snakemake workflow.

# Build a tempdir with an expression matrix + sample sheet and return a base
# config (no Zhao, no clustering) that writes to <td>/out. `extra_samples` adds
# columns to the sample sheet (to exercise metadata carry-through / clashes).
.make_run_case <- function(td, extra_samples = NULL) {
  ann <- pam50_annotation()
  set.seed(3)
  genes <- c(ann$GeneName, "XIST", "GAPDH")
  m <- matrix(round(rnorm(length(genes) * 6), 3), length(genes), 6,
              dimnames = list(genes, paste0("S", 1:6)))
  write.table(data.frame(gene = rownames(m), m, check.names = FALSE),
              file.path(td, "expr.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  ss <- data.frame(sample_id = paste0("S", 1:6), stringsAsFactors = FALSE)
  if (!is.null(extra_samples)) ss <- cbind(ss, extra_samples)
  write.table(ss, file.path(td, "samples.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  list(
    inputs = list(expression = file.path(td, "expr.tsv"),
                  samples = file.path(td, "samples.tsv"),
                  sample_id_column = "sample_id"),
    filter = list(column = NULL, include = list()),
    deduplicate = list(by = NULL),
    zhao_normalization = list(enabled = FALSE, subtype_column = NULL),
    classifier = list(calibration = NA, collapse_method = "mean",
                      has_clinical = FALSE, prefix = "PAM50scores"),
    intrinsic_clustering = list(enabled = FALSE),
    outs = list(output_dir = file.path(td, "out"))
  )
}

test_that("pam50_run writes a single clean call table by default", {
  td <- withr::local_tempdir()
  cfg <- .make_run_case(td)

  res <- suppressWarnings(pam50_run(cfg))
  expect_s3_class(res, "pam50_result")

  calls <- read.table(file.path(td, "out", "pam50_calls.tsv"),
                      sep = "\t", header = TRUE, check.names = FALSE)
  expect_equal(nrow(calls), 6L)

  # Clean, unique, snake_case columns — not the Perou-native names.
  expect_true("call" %in% colnames(calls))
  expect_false("Call" %in% colnames(calls))
  expect_false(any(duplicated(colnames(calls))))          # no duplicate 'Her2'
  expect_true(all(c("score_her2", "her2_signature", "er_signature") %in% colnames(calls)))

  # The legacy Perou file is NOT written unless asked for.
  expect_false(file.exists(file.path(td, "out", "PAM50scores_pam50scores.txt")))
})

test_that("legacy_perou_scores: true additionally emits the Perou-format file", {
  td <- withr::local_tempdir()
  cfg <- .make_run_case(td)
  cfg$classifier$legacy_perou_scores <- TRUE

  suppressWarnings(pam50_run(cfg))
  legacy <- file.path(td, "out", "PAM50scores_pam50scores.txt")
  expect_true(file.exists(legacy))
  # Legacy layout: unnamed first column (blank top-left cell) and quoted cells.
  hdr <- readLines(legacy, n = 1L)
  expect_match(hdr, '^""\t')
})

test_that("sample metadata is carried through; metadata_columns is an allowlist", {
  td <- withr::local_tempdir()
  cfg <- .make_run_case(td, extra_samples = data.frame(
    tissue_id = paste0("T", 1:6), batch = rep(c("a", "b"), 3), stringsAsFactors = FALSE))

  # Default: all metadata columns carried.
  suppressWarnings(pam50_run(cfg))
  calls <- read.table(file.path(td, "out", "pam50_calls.tsv"),
                      sep = "\t", header = TRUE, check.names = FALSE)
  expect_true(all(c("tissue_id", "batch") %in% colnames(calls)))

  # Allowlist: keep only tissue_id.
  cfg$inputs$metadata_columns <- list("tissue_id")
  suppressWarnings(pam50_run(cfg))
  calls <- read.table(file.path(td, "out", "pam50_calls.tsv"),
                      sep = "\t", header = TRUE, check.names = FALSE)
  expect_true("tissue_id" %in% colnames(calls))
  expect_false("batch" %in% colnames(calls))
})

test_that("a metadata column clashing with a classifier column is an error", {
  td <- withr::local_tempdir()
  cfg <- .make_run_case(td, extra_samples = data.frame(
    call = paste0("x", 1:6), stringsAsFactors = FALSE))     # 'call' shadows the call column
  expect_error(suppressWarnings(pam50_run(cfg)), "collide")
})

test_that("metadata_columns naming a missing column is an error", {
  td <- withr::local_tempdir()
  cfg <- .make_run_case(td)
  cfg$inputs$metadata_columns <- list("not_a_column")
  expect_error(suppressWarnings(pam50_run(cfg)), "not in sample metadata")
})

test_that("pam50_run writes an opt-in HTML report", {
  skip_if_not_installed("rmarkdown")
  skip_if_not_installed("knitr")
  skip_if(!nzchar(rmarkdown::pandoc_available()), "pandoc not available")

  td <- withr::local_tempdir()
  ann <- pam50_annotation()
  set.seed(3)
  genes <- c(ann$GeneName, "XIST", "GAPDH")
  m <- matrix(round(rnorm(length(genes) * 6), 3), length(genes), 6,
              dimnames = list(genes, paste0("S", 1:6)))
  write.table(data.frame(gene = rownames(m), m, check.names = FALSE),
              file.path(td, "expr.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(data.frame(sample_id = paste0("S", 1:6)),
              file.path(td, "samples.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

  cfg <- list(
    inputs = list(expression = file.path(td, "expr.tsv"),
                  samples = file.path(td, "samples.tsv"),
                  sample_id_column = "sample_id"),
    filter = list(column = NULL, include = list()),
    deduplicate = list(by = NULL),
    zhao_normalization = list(enabled = FALSE, subtype_column = NULL),
    classifier = list(calibration = NA, collapse_method = "mean",
                      has_clinical = FALSE, prefix = "PAM50scores"),
    intrinsic_clustering = list(enabled = FALSE),
    outs = list(output_dir = file.path(td, "out"))
  )

  suppressWarnings(pam50_run(cfg, report = TRUE))
  rp <- file.path(td, "out", "intrinsix_report.html")
  expect_true(file.exists(rp))
  expect_gt(file.info(rp)$size, 1000)
  html <- paste(readLines(rp, warn = FALSE), collapse = "\n")
  expect_match(html, "PAM50", fixed = TRUE)
  expect_match(html, "Results table", fixed = TRUE)
})

test_that("pam50_run validates required config keys", {
  expect_error(pam50_run(list(inputs = list())), "missing")
})

test_that("the bundled example config parses", {
  cfg <- yaml::read_yaml(system.file("extdata", "config.example.yaml", package = "intrinsix"))
  expect_true(!is.null(cfg$inputs$expression))
  expect_true(!is.null(cfg$outs$output_dir))
})

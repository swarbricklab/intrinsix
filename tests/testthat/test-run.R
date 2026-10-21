# Test the config-driven whole-pipeline runner (pam50_run), which replaces the
# former Snakemake workflow.

test_that("pam50_run executes a config and writes the call table", {
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

  res <- suppressWarnings(pam50_run(cfg))
  expect_s3_class(res, "pam50_result")

  calls <- read.table(file.path(td, "out", "pam50_calls.tsv"),
                      sep = "\t", header = TRUE, check.names = FALSE)
  expect_equal(nrow(calls), 6L)
  expect_true("Call" %in% colnames(calls))
  expect_true(file.exists(file.path(td, "out", "PAM50scores_pam50scores.txt")))
})

test_that("pam50_run validates required config keys", {
  expect_error(pam50_run(list(inputs = list())), "missing")
})

test_that("the bundled example config parses", {
  cfg <- yaml::read_yaml(system.file("extdata", "config.example.yaml", package = "intrinsix"))
  expect_true(!is.null(cfg$inputs$expression))
  expect_true(!is.null(cfg$outs$output_dir))
})

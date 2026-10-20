# Tests for the command-line interface. cli() reaches the same package functions
# the Snakemake rules call, so a classify via CLI must match predict_pam50 directly.

test_that("cli classify matches predict_pam50 and writes a scores table", {
  expr_path <- test_path("fixtures/example_expr_200.tsv")
  out <- tempfile(fileext = ".tsv")

  suppressWarnings(cli(c("classify", "--expr", expr_path, "--out", out, "--calibration", "NA")))
  # NB: "NA" parses to as.integer("NA") = NA -> median-centring path (example driver)

  via_cli <- read.table(out, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  expr <- as.matrix(read.table(expr_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE))
  direct <- suppressWarnings(predict_pam50(expr, calibration = NA, collapse = "mean"))

  expect_identical(as.character(via_cli[["Call"]]), as.character(direct[["Call"]]))
  expect_equal(nrow(via_cli), ncol(expr))
})

test_that("cli zhao reproduces zhao_center", {
  expr_path <- test_path("fixtures/synthetic_expr.tsv")
  samples_path <- test_path("fixtures/synthetic_samples.tsv")
  out <- tempfile(fileext = ".tsv")

  cli(c("zhao", "--expr", expr_path, "--samples", samples_path,
        "--subtype-col", "clinical_subtype", "--id-col", "sample", "--out", out))

  via_cli <- as.matrix(read.table(out, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE))
  exp <- as.matrix(read.table(test_path("fixtures/expected_zhao.tsv"),
                              sep = "\t", header = TRUE, row.names = 1, check.names = FALSE))
  expect_equal(via_cli[rownames(exp), colnames(exp)], exp, tolerance = 1e-10)
})

test_that("cli errors on missing required options and unknown command", {
  expect_error(cli(c("classify", "--expr", "x.tsv")), "missing required")
  expect_error(cli(c("frobnicate")), "unknown command")
})

test_that("cli help returns 0", {
  expect_equal(suppressMessages(cli(c("--help"))), 0L)
})

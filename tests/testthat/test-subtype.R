# Tests for the high-level pam50_subtype() interface + the pam50_result S3 class.

test_that("pam50_subtype reproduces predict_pam50 / the golden calls", {
  expr <- as.matrix(read.table(test_path("fixtures/example_expr_200.tsv"),
                               sep = "\t", header = TRUE, row.names = 1, check.names = FALSE))
  exp <- read.table(test_path("fixtures/expected_pam50scores_200.tsv"),
                    sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

  res <- suppressWarnings(pam50_subtype(expr, map_symbols = FALSE, zhao = FALSE, calibration = NA))

  expect_s3_class(res, "pam50_result")
  expect_identical(as.character(res$calls[rownames(exp), "Call"]),
                   as.character(exp[["Call"]]))
})

test_that("pam50_result S3 methods behave", {
  expr <- as.matrix(read.table(test_path("fixtures/example_expr_200.tsv"),
                               sep = "\t", header = TRUE, row.names = 1, check.names = FALSE))
  res <- suppressWarnings(pam50_subtype(expr, map_symbols = FALSE, zhao = FALSE, calibration = NA))

  expect_output(print(res), "pam50_result")
  s <- summary(res)
  expect_equal(sum(s$subtype), 200L)
  expect_equal(nrow(as.data.frame(res)), 200L)

  pf <- tempfile(fileext = ".pdf")
  pdf(pf); on.exit(unlink(pf))
  expect_silent(plot(res))
  dev.off()
  expect_gt(file.size(pf), 1000)
})

test_that("pam50_subtype runs the Zhao path", {
  se <- as.matrix(read.table(test_path("fixtures/synthetic_expr.tsv"),
                             sep = "\t", header = TRUE, row.names = 1, check.names = FALSE))
  sm <- read.table(test_path("fixtures/synthetic_samples.tsv"),
                   sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  res <- suppressWarnings(pam50_subtype(se, samples = sm, subtype_col = "clinical_subtype",
                                        id_col = "sample", map_symbols = FALSE))
  expect_true(res$params$zhao)
  expect_equal(nrow(res$calls), 6L)
})

test_that("pam50_subtype errors when zhao is requested without samples", {
  se <- as.matrix(read.table(test_path("fixtures/synthetic_expr.tsv"),
                             sep = "\t", header = TRUE, row.names = 1, check.names = FALSE))
  expect_error(pam50_subtype(se, zhao = TRUE, map_symbols = FALSE), "requires")
})

# Golden regression test for the PAM50 classifier. The fixture is the PUBLIC
# upstream 200-subject example (perou-lab/PAM50_classifier/PAM50_example), run
# with the example driver's settings (calibrationParameters = NA). The expected
# output was frozen from this package after it was verified byte-identical to the
# production bcarc_wts pipeline on real (uncommittable) data.

test_that("predict_pam50 reproduces the frozen calls on the public example", {
  expr <- as.matrix(read.table(
    test_path("fixtures/example_expr_200.tsv"),
    sep = "\t", header = TRUE, row.names = 1, check.names = FALSE
  ))
  exp <- read.table(
    test_path("fixtures/expected_pam50scores_200.tsv"),
    sep = "\t", header = TRUE, row.names = 1, check.names = FALSE
  )

  got <- suppressWarnings(predict_pam50(expr, calibration = NA, collapse = "mean"))
  got <- got[rownames(exp), ]

  # Headline: identical subtype calls.
  expect_identical(got[["Call"]], exp[["Call"]])

  # Numeric columns identical to float tolerance.
  numeric_cols <- c("Basal", "Her2", "LumA", "LumB", "Normal",
                    "Confidence", "ROR-S (Subtype Only)", "Proliferation Score",
                    "ROR-P (Subtype + Proliferation)", "ER", "Her2")
  for (col in numeric_cols) {
    expect_equal(as.numeric(got[[col]]), as.numeric(exp[[col]]),
                 tolerance = 1e-8, info = col)
  }

  # Risk-group labels identical.
  for (col in c("ROR-S Group (Subtype Only)", "ROR-P Group (Subtype + Proliferation)")) {
    expect_identical(got[[col]], exp[[col]], info = col)
  }
})

test_that("predict_pam50 rejects unsupported ROR-C mode", {
  expr <- as.matrix(read.table(
    test_path("fixtures/example_expr_200.tsv"),
    sep = "\t", header = TRUE, row.names = 1, check.names = FALSE
  ))
  expect_error(predict_pam50(expr, has_clinical = TRUE), "not yet supported")
})

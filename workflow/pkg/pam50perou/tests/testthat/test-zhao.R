# Regression test for Zhao4 subtype-specific quantile centering. Uses a small
# deterministic synthetic fixture (zhao_center is a pure function of expression +
# subtypes + SIGMA). The port was separately verified to reproduce the production
# bcarc_wts pam50_zhao_normalized.tsv to ~1e-14 on real data.

test_that("zhao_center reproduces the frozen synthetic output", {
  m <- as.matrix(read.table(
    test_path("fixtures/synthetic_expr.tsv"),
    sep = "\t", header = TRUE, row.names = 1, check.names = FALSE
  ))
  samples <- read.table(
    test_path("fixtures/synthetic_samples.tsv"),
    sep = "\t", header = TRUE, stringsAsFactors = FALSE
  )
  subt <- stats::setNames(samples$clinical_subtype, samples$sample)

  exp <- as.matrix(read.table(
    test_path("fixtures/expected_zhao.tsv"),
    sep = "\t", header = TRUE, row.names = 1, check.names = FALSE
  ))

  got <- zhao_center(m, subt)
  got <- got[rownames(exp), colnames(exp)]

  expect_equal(got, exp, tolerance = 1e-10)
})

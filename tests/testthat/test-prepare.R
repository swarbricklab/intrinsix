# Unit tests for the input-preparation helpers.

test_that("map_to_pam50 relabels HGNC symbols to pcrID and subsets to 50 genes", {
  ann <- pam50_annotation()
  # Build a matrix labelled by HGNC GeneName, plus a non-PAM50 decoy row.
  mat <- matrix(
    rnorm((nrow(ann) + 1) * 3),
    nrow = nrow(ann) + 1,
    dimnames = list(c(ann$GeneName, "NOT_A_PAM50_GENE"), c("A", "B", "C"))
  )

  got <- suppressWarnings(map_to_pam50(mat, annotation = ann))

  expect_equal(nrow(got), 50L)
  expect_setequal(rownames(got), ann$pcrID)
  expect_false("NOT_A_PAM50_GENE" %in% rownames(got))
  expect_equal(ncol(got), 3L)
})

test_that("map_to_pam50 warns when PAM50 genes are missing", {
  ann <- pam50_annotation()
  mat <- matrix(
    rnorm(5 * 2), nrow = 5,
    dimnames = list(ann$GeneName[1:5], c("A", "B"))
  )
  expect_warning(map_to_pam50(mat, annotation = ann), "only 5/50")
})

test_that("median_center subtracts row medians", {
  x <- matrix(c(1, 2, 3, 10, 20, 30), nrow = 2, byrow = TRUE)
  mc <- median_center(x)
  expect_equal(apply(mc, 1, median), c(0, 0))
})

test_that("reference-data accessors load the bundled files", {
  expect_true(all(c("pcrID", "GeneName") %in% names(pam50_annotation())))
  expect_equal(ncol(pam50_centroids()), 5L)
  expect_equal(nrow(pam50_centroids()), 50L)
})

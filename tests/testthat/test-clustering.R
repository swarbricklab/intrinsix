# Tests for the Cluster 3.0 replication + claudin-low identification. The port was
# separately verified against the production bcarc_wts run on real data: the array
# and gene merge trees were byte-identical in topology (0/224 and 0/1430 nodes
# differing; similarity values matched to ~1e-15), and the 7 claudin-low calls
# were reproduced exactly. These tests use a small deterministic synthetic fixture
# with a planted claudin-low cluster (CL1-CL3).

test_that("centroid_clustering reproduces the frozen tree topology", {
  m <- as.matrix(read.table(
    test_path("fixtures/clustering_expr.tsv"),
    sep = "\t", header = TRUE, row.names = 1, check.names = FALSE
  ))
  res <- centroid_clustering(m, write_files = FALSE)

  exp_atr <- read.table(test_path("fixtures/expected_clustering_atr.tsv"),
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  exp_gtr <- read.table(test_path("fixtures/expected_clustering_gtr.tsv"),
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)

  expect_identical(res$atr[, c("NODE", "LEFT", "RIGHT")], exp_atr)
  expect_identical(res$gtr[, c("NODE", "LEFT", "RIGHT")], exp_gtr)
})

test_that("centroid_clustering is deterministic", {
  m <- as.matrix(read.table(
    test_path("fixtures/clustering_expr.tsv"),
    sep = "\t", header = TRUE, row.names = 1, check.names = FALSE
  ))
  a <- centroid_clustering(m, write_files = FALSE)
  b <- centroid_clustering(m, write_files = FALSE)
  expect_equal(a$atr, b$atr)
  expect_equal(a$gtr, b$gtr)
  expect_equal(a$cdt, b$cdt)
})

test_that("identify_claudin_low recovers the planted claudin-low samples", {
  m <- as.matrix(read.table(
    test_path("fixtures/clustering_expr.tsv"),
    sep = "\t", header = TRUE, row.names = 1, check.names = FALSE
  ))
  res <- centroid_clustering(m, write_files = FALSE)
  clow <- identify_claudin_low(res, m)

  expected <- readLines(test_path("fixtures/expected_claudin_low.txt"))
  expect_identical(clow, expected)
})

test_that("identify_claudin_low errors when claudin genes are absent", {
  m <- as.matrix(read.table(
    test_path("fixtures/clustering_expr.tsv"),
    sep = "\t", header = TRUE, row.names = 1, check.names = FALSE
  ))
  m2 <- m[!rownames(m) %in% c("CLDN3", "CLDN4", "CLDN7"), ]
  res <- centroid_clustering(m2, write_files = FALSE)
  expect_error(identify_claudin_low(res, m2), "claudin genes")
})

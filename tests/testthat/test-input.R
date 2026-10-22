# The input layer: a SummarizedExperiment / AnnData / .h5ad must give exactly the
# same calls as the equivalent plain matrix. Skipped where the optional reader
# packages aren't installed (they're Suggests).

.example_matrix <- function() {
  as.matrix(read.table(test_path("fixtures/example_expr_200.tsv"),
                       sep = "\t", header = TRUE, row.names = 1, check.names = FALSE))
}

.matrix_calls <- function(m) {
  as.character(suppressWarnings(
    pam50_subtype(m, map_symbols = FALSE, zhao = FALSE, calibration = NA))$calls$Call)
}

test_that("SummarizedExperiment input matches the matrix path", {
  skip_if_not_installed("SummarizedExperiment")
  m <- .example_matrix()
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(log2_uqn = m),
    colData = S4Vectors::DataFrame(sample_id = colnames(m), row.names = colnames(m)))
  r <- suppressWarnings(pam50_subtype(se, layer = "log2_uqn", map_symbols = FALSE,
                                      zhao = FALSE, calibration = NA))
  expect_identical(as.character(r$calls$Call), .matrix_calls(m))
})

test_that(".pcrid_from_entrez round-trips the bundled annotation", {
  ann <- pam50_annotation()
  expect_identical(unname(intrinsix:::.pcrid_from_entrez(ann$EntrezGene)),
                   as.character(ann$pcrID))
})

test_that("entrez_col selects the same PAM50 genes as pam50_id_col", {
  skip_if_not_installed("SummarizedExperiment")
  ann <- pam50_annotation()
  set.seed(7)
  genes  <- c(as.character(ann$pcrID), paste0("DECOY", 1:5))
  entrez <- c(as.character(ann$EntrezGene), rep(NA_character_, 5))
  pcr    <- c(as.character(ann$pcrID), rep("", 5))
  m <- matrix(rnorm(length(genes) * 8), length(genes), 8,
              dimnames = list(genes, paste0("S", 1:8)))
  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(x = m),
    rowData = S4Vectors::DataFrame(entrez_id = entrez, pam50_symbol = pcr),
    colData = S4Vectors::DataFrame(sample_id = colnames(m), row.names = colnames(m)))
  a <- suppressWarnings(pam50_subtype(se, entrez_col = "entrez_id",
                                      map_symbols = FALSE, zhao = FALSE, calibration = NA))
  b <- suppressWarnings(pam50_subtype(se, pam50_id_col = "pam50_symbol",
                                      map_symbols = FALSE, zhao = FALSE, calibration = NA))
  expect_identical(as.character(a$calls$Call), as.character(b$calls$Call))
})

test_that("h5ad / AnnData input matches the matrix path", {
  skip_if_not_installed("anndataR")
  m <- .example_matrix()
  ad <- anndataR::AnnData(
    X = t(m),
    obs = data.frame(row.names = colnames(m)),
    var = data.frame(gene_name = rownames(m), row.names = rownames(m)))
  h5 <- withr::local_tempfile(fileext = ".h5ad")
  anndataR::write_h5ad(ad, h5)
  r <- suppressWarnings(pam50_subtype(h5, map_symbols = FALSE, zhao = FALSE, calibration = NA))
  expect_identical(as.character(r$calls[colnames(m), "Call"]), .matrix_calls(m))
})

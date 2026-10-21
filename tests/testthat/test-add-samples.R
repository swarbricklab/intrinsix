# Tests for incremental calling (pam50_add_samples): reference calls stay fixed,
# new samples are recalled one at a time against the reference.

.load_syn <- function() {
  m <- as.matrix(read.table(test_path("fixtures/synthetic_expr.tsv"),
                            sep = "\t", header = TRUE, row.names = 1, check.names = FALSE))
  sm <- read.table(test_path("fixtures/synthetic_samples.tsv"),
                   sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  list(m = m, sm = sm)
}

test_that("reference calls stay fixed and new samples are recalled singly", {
  d <- .load_syn()
  ref_ids <- c("S1", "S2", "S3", "S4"); new_ids <- c("S5", "S6")
  ref_expr <- d$m[, ref_ids]; new_expr <- d$m[, new_ids, drop = FALSE]
  ref_sm <- d$sm[d$sm$sample %in% ref_ids, ]; new_sm <- d$sm[d$sm$sample %in% new_ids, ]

  out <- suppressWarnings(pam50_add_samples(
    ref_expr, ref_sm, new_expr, new_sm,
    subtype_col = "clinical_subtype", id_col = "sample", map_symbols = FALSE))

  expect_equal(nrow(out), 6L)
  expect_equal(out$added, c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE))

  # reference calls identical to a reference-only run
  ref_only <- suppressWarnings(pam50_subtype(
    ref_expr, samples = ref_sm, subtype_col = "clinical_subtype",
    id_col = "sample", map_symbols = FALSE))
  expect_identical(as.character(out[ref_ids, "Call"]),
                   as.character(ref_only$calls[ref_ids, "Call"]))

  # each new sample equals its own (reference + 1) augmented run
  for (nid in new_ids) {
    aug <- suppressWarnings(pam50_subtype(
      cbind(ref_expr, new_expr[, nid, drop = FALSE]),
      samples = rbind(ref_sm, new_sm[new_sm$sample == nid, ]),
      subtype_col = "clinical_subtype", id_col = "sample", map_symbols = FALSE))
    expect_identical(as.character(out[nid, "Call"]),
                     as.character(aug$calls[nid, "Call"]))
  }
})

test_that("claudin_low = TRUE recalls CLOW incrementally with reference fixed", {
  pm_genes <- rownames(pam50_centroids())
  intr <- intrinsic_genes()$Name
  intr_sub <- unique(c("CLDN3", "CLDN4", "CLDN7",
                       setdiff(intr, c("CLDN3", "CLDN4", "CLDN7"))[1:22]))
  samples <- paste0("S", 1:10)
  clow <- c("S1", "S2", "S3")
  set.seed(11)
  pm <- matrix(round(rnorm(length(pm_genes) * 10), 3), length(pm_genes), 10,
               dimnames = list(pm_genes, samples))
  im <- matrix(round(rnorm(length(intr_sub) * 10, 0, 0.5), 3), length(intr_sub), 10,
               dimnames = list(intr_sub, samples))
  for (g in c("CLDN3", "CLDN4", "CLDN7")) {
    im[g, ] <- ifelse(samples %in% clow, -2, 1) + rnorm(10, 0, 0.02)  # correlated low block
  }
  expr <- rbind(pm, im[!rownames(im) %in% rownames(pm), , drop = FALSE])
  sm <- data.frame(sample = samples,
                   clinical_subtype = rep(c("ERpos_HER2neg", "ERneg_HER2neg"), 5),
                   stringsAsFactors = FALSE)

  ref_ids <- paste0("S", 1:7); new_ids <- paste0("S", 8:10)
  out <- suppressWarnings(pam50_add_samples(
    expr[, ref_ids], sm[sm$sample %in% ref_ids, ],
    expr[, new_ids], sm[sm$sample %in% new_ids, ],
    subtype_col = "clinical_subtype", id_col = "sample",
    map_symbols = FALSE, claudin_low = TRUE))

  expect_equal(nrow(out), 10L)
  expect_true("Call_CLOW" %in% colnames(out))
  expect_equal(out$added, c(rep(FALSE, 7), rep(TRUE, 3)))

  # reference CLOW assignments unchanged vs a reference-only clustered run
  ref_only <- suppressWarnings(pam50_subtype(
    expr[, ref_ids], samples = sm[sm$sample %in% ref_ids, ],
    subtype_col = "clinical_subtype", id_col = "sample",
    map_symbols = FALSE, cluster = TRUE))
  expect_identical(as.character(out[ref_ids, "Call_CLOW"]),
                   as.character(ref_only$calls[ref_ids, "Call_CLOW"]))
})

test_that("pam50_add_samples rejects new/reference id collisions", {
  d <- .load_syn()
  expect_error(
    suppressWarnings(pam50_add_samples(
      d$m[, c("S1", "S2")], d$sm[d$sm$sample %in% c("S1", "S2"), ],
      d$m[, "S1", drop = FALSE], d$sm[d$sm$sample == "S1", ],
      subtype_col = "clinical_subtype", id_col = "sample", map_symbols = FALSE)),
    "collide")
})

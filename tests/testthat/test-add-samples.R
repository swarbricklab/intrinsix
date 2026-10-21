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

test_that("pam50_add_samples rejects new/reference id collisions", {
  d <- .load_syn()
  expect_error(
    suppressWarnings(pam50_add_samples(
      d$m[, c("S1", "S2")], d$sm[d$sm$sample %in% c("S1", "S2"), ],
      d$m[, "S1", drop = FALSE], d$sm[d$sm$sample == "S1", ],
      subtype_col = "clinical_subtype", id_col = "sample", map_symbols = FALSE)),
    "collide")
})

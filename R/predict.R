# PAM50 subtype prediction. Faithful transcription of the Perou Lab driver
# (vendor/perou-lab/subtypePrediction_distributed.R, commit 45121d8), refactored
# to take explicit arguments instead of reading global variables. The scoring
# math is byte-for-byte the same; only the I/O and plotting were lifted out.
#
# .classify_core() does the work and returns the scores table PLUS the
# intermediates needed for plotting (distances, standardised test matrix,
# proliferation score, ...). predict_pam50() is the public wrapper that returns
# just the scores data frame (its output is verified byte-identical to the
# production pipeline).

# @keywords internal
# @noRd
.classify_core <- function(expr, calibration = -1, collapse = "mean",
                           has_clinical = FALSE, centroids = pam50_centroids()) {
  if (has_clinical) {
    stop("has_clinical = TRUE (ROR-C) is not yet supported; it needs a tumour ",
         "T-stage vector. Use has_clinical = FALSE.", call. = FALSE)
  }

  proliferationGenes <- c("CCNB1", "UBE2C", "BIRC5", "KNTC2", "CDC20", "PTTG1",
                          "RRM2", "MKI67", "TYMS", "CEP55", "CDCA1")

  glthreshold <- -0.15; ghthreshold <- 0.1
  gplthreshold <- -0.25; gphthreshold <- 0.1

  # test matrix (readarray-equivalent: numeric + collapse dups)
  xd <- as.matrix(expr)
  storage.mode(xd) <- "double"
  xd <- collapseIDs(xd, rownames(xd), collapse)

  # normalisation (mirrors the driver's calibration block)
  if (is.na(calibration)) {
    xd <- medianCtr(xd)
  } else if (calibration != -1) {
    medians <- readarray(
      system.file("extdata", "mediansPerDataset_v2.txt", package = "pam50"),
      hr = 1
    )
    tm <- overlapSets(medians$xd, xd)
    xd <- tm$y - tm$x[, calibration]
  }
  # calibration == -1: no centring here (input already normalised), only std below.

  xd <- standardize(xd)

  erScore <- as.vector(t(xd["ESR1", ]))
  her2Score <- as.vector(t(xd["ERBB2", ]))

  this.proliferationGenes <- dimnames(xd)[[1]] %in% proliferationGenes
  prolifScore <- apply(xd[this.proliferationGenes, ], 2, mean, na.rm = TRUE)

  out <- sspPredict(centroids, classes = "", xd, std = FALSE,
                    distm = "spearman", centroids = TRUE)
  out$distances <- -1 * out$distances

  call.conf <- c()
  for (j in 1:length(out$predictions)) {
    call.conf[j] <- 1 - stats::cor.test(
      out$testData[, j],
      out$centroids[, which(colnames(centroids) == out$predictions[j])],
      method = "spearman"
    )$p.value
  }
  call.conf <- round(call.conf, 2)

  genomic <- 0.04210193 * out$distances[, 1] + 0.12466938 * out$distances[, 2] +
    -0.35235561 * out$distances[, 3] + 0.14213283 * out$distances[, 4]
  genomicWprolif <- -0.0009299747 * out$distances[, 1] + 0.0692289192 * out$distances[, 2] +
    -0.0951505484 * out$distances[, 3] + 0.0493487685 * out$distances[, 4] +
    0.3385116381 * prolifScore

  griskgroups <- genomic
  griskgroups[genomic > ghthreshold] <- "high"
  griskgroups[genomic > glthreshold & genomic < ghthreshold] <- "med"
  griskgroups[genomic < glthreshold] <- "low"
  gpriskgroups <- genomicWprolif
  gpriskgroups[genomicWprolif > gphthreshold] <- "high"
  gpriskgroups[genomicWprolif > gplthreshold & genomicWprolif < gphthreshold] <- "med"
  gpriskgroups[genomicWprolif < gplthreshold] <- "low"

  genomic <- 100 * (genomic + 0.35) / 0.85
  genomicWprolif <- 100 * (genomicWprolif + 0.35) / 0.85

  outtable <- cbind(out$distances, out$predictions, call.conf, genomic,
                    griskgroups, prolifScore, genomicWprolif, gpriskgroups,
                    erScore, her2Score)
  dimnames(outtable)[[2]] <- c(
    "Basal", "Her2", "LumA", "LumB", "Normal", "Call", "Confidence",
    "ROR-S (Subtype Only)", "ROR-S Group (Subtype Only)", "Proliferation Score",
    "ROR-P (Subtype + Proliferation)", "ROR-P Group (Subtype + Proliferation)",
    "ER", "Her2"
  )
  rownames(outtable) <- colnames(out$testData)
  scores <- as.data.frame(outtable, stringsAsFactors = FALSE, check.names = FALSE)

  list(
    scores = scores,
    testData = out$testData,
    centroids = out$centroids,
    distances = out$distances,
    predictions = out$predictions,
    prolifScore = prolifScore,
    erScore = erScore,
    her2Score = her2Score,
    call_conf = call.conf
  )
}

#' PAM50 intrinsic subtype prediction
#'
#' Classifies samples into the five PAM50 intrinsic subtypes (Basal, Her2, LumA,
#' LumB, Normal) by Spearman rank correlation to the published Perou centroids,
#' and computes ROR-S and ROR-P risk scores. This is the exact algorithm used by
#' the Swarbrick Lab PAM50 pipeline, packaged (not reimplemented).
#'
#' For a one-call interface from a raw expression matrix (gene-symbol mapping +
#' optional Zhao normalisation + prediction), see [pam50_subtype()].
#'
#' @param expr Numeric matrix, genes x samples. Row names are the classifier's
#'   pcrID gene symbols (see [map_to_pam50()]); column names are sample IDs.
#'   Values should already be normalised (e.g. Zhao4-centred, or log2-scale for
#'   the median-centring path).
#' @param calibration Perou `calibrationParameters`. `-1` (default) applies no
#'   calibration and no median-centring (used as-is, then column-standardised);
#'   `NA` median-centres; a positive integer selects a calibration column from the
#'   bundled `mediansPerDataset_v2.txt`.
#' @param collapse Duplicate-gene collapse method, `"mean"` (default) or `"iqr"`.
#' @param has_clinical Logical. `FALSE` (default) omits ROR-C/ROR-PC (which need
#'   tumour T-stage). `TRUE` is not yet supported here.
#' @param centroids Centroid table (genes x 5 subtypes); defaults to
#'   [pam50_centroids()].
#' @return A data frame, one row per sample (row names = sample IDs), with
#'   columns Basal/Her2/LumA/LumB/Normal (centroid distances), Call, Confidence,
#'   ROR-S + group, Proliferation Score, ROR-P + group, ER, Her2.
#' @seealso [pam50_subtype()] for the high-level interface.
#' @export
predict_pam50 <- function(expr, calibration = -1, collapse = "mean",
                          has_clinical = FALSE, centroids = pam50_centroids()) {
  .classify_core(expr, calibration = calibration, collapse = collapse,
                 has_clinical = has_clinical, centroids = centroids)$scores
}

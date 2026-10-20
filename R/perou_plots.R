# Visualisations ported from the Perou Lab scripts
# (vendor/perou-lab/subtypePrediction_functions.R) and from Alec Lobanov's
# intrinsic-clustering script. The internal .myplot/.pcaEA helpers are kept close
# to the originals; the public plot_* functions wrap them to take a pam50_result.

# --- internal: per-subtype correlation + proliferation boxplots (myplot) ---
# @noRd
.myplot <- function(y, short, pro) {
  op <- graphics::par(mfrow = c(3, 2), mar = c(5, 3, 2, 2), las = 3)
  on.exit(graphics::par(op))
  y$prediction <- factor(y$prediction, levels = c("Basal", "Her2", "LumA", "LumB", "Normal"))
  labs <- c("Basal", "Her2", "LumA", "LumB", "Normal")
  for (k in 1:5) {
    graphics::boxplot(y$distances[, k] ~ y$prediction, border = 8,
                      ylab = paste(labs[k], "Correlation"),
                      main = paste(short, ": ", labs[k], sep = ""))
    graphics::stripchart(y$distances[, k] ~ y$prediction, vertical = TRUE,
                         method = "jitter", pch = 3, add = TRUE)
  }
  graphics::boxplot(pro ~ y$prediction, border = 8, ylab = "Proliferation Index",
                    main = paste(short, ": Proliferation Index", sep = ""))
  graphics::stripchart(pro ~ y$prediction, vertical = TRUE, method = "jitter",
                       pch = 3, add = TRUE)
  invisible(NULL)
}

#' Plot PAM50 score distributions
#'
#' Six-panel figure (ported from the Perou Lab `myplot`): for each intrinsic
#' subtype, a boxplot of its centroid correlation across samples grouped by the
#' predicted call, plus the proliferation index. Useful for eyeballing call
#' confidence and separation.
#'
#' @param x A [pam50_result] from [pam50_subtype()].
#' @param short Title prefix (default `"PAM50"`).
#' @param ... Ignored.
#' @return `x`, invisibly (called for the plot side effect).
#' @export
plot_pam50_scores <- function(x, short = "PAM50", ...) {
  if (!inherits(x, "pam50_result")) {
    stop("plot_pam50_scores() expects a pam50_result (from pam50_subtype()).", call. = FALSE)
  }
  cl <- x$classify
  .myplot(list(distances = cl$distances, prediction = cl$predictions), short, cl$prolifScore)
  invisible(x)
}

# --- internal: PCA scatter coloured by class (pcaEA) ---
# @noRd
.pcaEA <- function(x, classes, size = 1, showLegend = TRUE, legendloc = "topright",
                   mainStr = "", startPC = 1, stopPC = 2, showNames = TRUE,
                   showClasses = FALSE, axisExpansion = 0, groupColors = NA) {
  x <- apply(x, 2, as.numeric)
  data.pca <- stats::prcomp(as.matrix(x))
  tmp <- data.pca$sdev[1:10]^2 / sum(data.pca$sdev^2)
  gr.labels <- as.vector(t(classes))
  gr.labels.fac <- factor(as.vector(t(classes)), exclude = "")
  nlabels <- nlevels(gr.labels.fac)
  legendLabels <- vector(); legendColors <- vector()
  for (k in 1:nlabels) {
    group <- levels(gr.labels.fac)[k]
    legendLabels[k] <- group
    if (length(groupColors) > 1) {
      gr.labels[gr.labels.fac == group] <- groupColors[k]; legendColors[k] <- groupColors[k]
    } else {
      gr.labels[gr.labels.fac == group] <- k; legendColors[k] <- k
    }
  }
  if (length(groupColors) == 1) gr.labels <- as.numeric(gr.labels)
  i <- startPC; j <- stopPC
  graphics::par(lab = c(3, 4, 3), mgp = c(.3, .5, .0), mai = c(.5, .5, .5, .5),
                xaxt = "n", yaxt = "n")
  strX <- paste("PC", i, paste("(", round(tmp[i], 4) * 100, "%)", sep = ""), sep = " ")
  strY <- paste("PC", j, paste("(", round(tmp[j], 4) * 100, "%)", sep = ""), sep = " ")
  graphics::plot(data.pca$rotation[, i], data.pca$rotation[, j], xlab = strX, ylab = strY,
                 main = mainStr, col = gr.labels, cex = size, pch = "")
  graphics::points(data.pca$rotation[, i], data.pca$rotation[, j], col = gr.labels,
                   cex = size * 1.5, pch = 19)
  if (showLegend) {
    graphics::legend(legendloc, legend = legendLabels, col = legendColors, pch = 19,
                     x.intersp = .3, yjust = .5, bty = "n", cex = size)
  }
}

#' Plot a PAM50 training-vs-test PCA
#'
#' Three-panel PCA (ported from the Perou Lab driver): the published training set
#' plus the classified test samples, coloured by cohort, by training subtype, and
#' by predicted subtype. Shows where the test cohort sits relative to the
#' reference. Uses the bundled training set.
#'
#' @param x A [pam50_result] from [pam50_subtype()].
#' @param short Label for the test cohort (default `"PAM50"`).
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
plot_pam50_pca <- function(x, short = "PAM50", ...) {
  if (!inherits(x, "pam50_result")) {
    stop("plot_pam50_pca() expects a pam50_result.", call. = FALSE)
  }
  subtypeColors <- c("red", "hotpink", "darkblue", "skyblue", "green")
  xtrain <- readarray(training_set(), hr = 2)
  xtrain$xd <- standardize(medianCtr(xtrain$xd))
  y_xd <- x$classify$xd
  preds <- x$classify$predictions

  tm <- overlapSets(xtrain$xd, y_xd)
  keep <- !is.na(xtrain$classes$subtype)
  tm$x <- tm$x[, keep]
  if (requireNamespace("impute", quietly = TRUE)) {
    yimp <- impute::impute.knn(as.matrix(tm$y))$data
  } else {
    yimp <- as.matrix(tm$y)
  }
  M <- cbind(tm$x, yimp)
  nTrain <- sum(keep)
  nTest <- ncol(y_xd)
  cohort <- c(rep("train", nTrain), rep(short, nTest))
  trainSub <- c(as.character(xtrain$classes$subtype[keep]), rep(NA, ncol(M) - nTrain))
  testSub <- c(rep(NA, ncol(M) - length(preds)), preds)
  M <- scale(M, center = FALSE)

  op <- graphics::par(mfrow = c(1, 3)); on.exit(graphics::par(op))
  .pcaEA(M, cohort, mainStr = "Training + test", showNames = FALSE, showClasses = FALSE)
  .pcaEA(M[, !is.na(trainSub)], trainSub[!is.na(trainSub)], mainStr = "Training subtypes",
         showNames = FALSE, showClasses = FALSE, groupColors = subtypeColors)
  .pcaEA(M[, !is.na(testSub)], testSub[!is.na(testSub)], mainStr = "Test calls",
         showNames = FALSE, showClasses = FALSE, groupColors = subtypeColors)
  invisible(x)
}

#' Intrinsic-subtype heatmap
#'
#' Draws the median-centred intrinsic expression heatmap (ComplexHeatmap), with
#' optional per-sample annotations and Cluster 3.0 dendrograms. This is the same
#' function the workflow's intrinsic-clustering step uses.
#'
#' @param expr_medctr Median-centred intrinsic expression matrix (genes x samples).
#' @param annotations Optional data frame of per-sample annotation columns (rows
#'   aligned to `colnames(expr_medctr)`), e.g. subtype / claudin-low / silhouette.
#' @param cluster_dir Optional directory holding `intrinsic_expr_medCtr.atr`/`.gtr`
#'   (from [centroid_clustering()]); used for the dendrograms. If absent,
#'   ComplexHeatmap clusters internally.
#' @param title Optional plot title.
#' @param file Optional PDF path; if given, the heatmap is written there and `x`
#'   is returned invisibly, else the Heatmap object is returned.
#' @return A ComplexHeatmap object (or invisibly, when `file` is written).
#' @export
pam50_heatmap <- function(expr_medctr, annotations = NULL, cluster_dir = NULL,
                          title = NULL, file = NULL) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) ||
      !requireNamespace("circlize", quietly = TRUE)) {
    stop("pam50_heatmap() needs the ComplexHeatmap and circlize packages.", call. = FALSE)
  }
  col_dend <- TRUE; row_dend <- TRUE
  if (!is.null(cluster_dir) && requireNamespace("ctc", quietly = TRUE)) {
    col_dend <- tryCatch(ctc::xcluster2r(file = file.path(cluster_dir, "intrinsic_expr_medCtr.atr"),
                         distance = "pearson", labels = FALSE, fast = FALSE, clean = TRUE, dec = "."),
                         error = function(e) TRUE)
    row_dend <- tryCatch(ctc::xcluster2r(file = file.path(cluster_dir, "intrinsic_expr_medCtr.gtr"),
                         distance = "pearson", labels = FALSE, fast = FALSE, clean = TRUE, dec = "."),
                         error = function(e) TRUE)
  }
  top <- NULL
  if (!is.null(annotations)) {
    top <- ComplexHeatmap::HeatmapAnnotation(df = as.data.frame(annotations, check.names = FALSE))
  }
  hm <- ComplexHeatmap::Heatmap(
    expr_medctr, name = "Gene exp.",
    show_column_names = FALSE, show_row_names = FALSE,
    top_annotation = top, row_title = "genes", column_title = "samples",
    cluster_rows = row_dend, cluster_columns = col_dend,
    col = circlize::colorRamp2(c(-4, 0, 4), c("blue", "black", "yellow")),
    use_raster = TRUE
  )
  if (!is.null(file)) {
    grDevices::pdf(file, width = 11, height = 8)
    ComplexHeatmap::draw(hm, column_title = title)
    grDevices::dev.off()
    return(invisible(hm))
  }
  hm
}


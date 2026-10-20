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

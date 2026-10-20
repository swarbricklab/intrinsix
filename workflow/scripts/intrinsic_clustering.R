# Snakemake shim: intrinsic clustering, claudin-low identification, silhouette
# widths, and the intrinsic heatmap. Clustering + claudin-low are in pam50;
# the heatmap/silhouette plotting (ComplexHeatmap / ctc / cluster) stays here.
# Unlike the legacy pipeline this writes a NEW call table rather than overwriting
# an intermediate in place.

suppressPackageStartupMessages({
  library(pam50)
  library(cluster)
})

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message"); sink(log, type = "output")

cfg <- snakemake@config
id_col <- cfg$inputs$sample_id_column
subtype_col <- cfg$zhao_normalization$subtype_column

# --- inputs ---
scores <- read.table(snakemake@input[["scores"]],
                     sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
samples <- read.table(snakemake@input[["samples_filtered"]],
                      sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
rownames(samples) <- samples[[id_col]]
int_expr <- as.matrix(read.table(snakemake@input[["int_expr"]],
                                 sep = "\t", header = TRUE, row.names = 1, check.names = FALSE))
int_medctr <- as.matrix(read.table(snakemake@input[["int_expr_medctr"]],
                                   sep = "\t", header = TRUE, row.names = 1, check.names = FALSE))

heatmap_out <- snakemake@output[["heatmap"]]
hdir <- dirname(heatmap_out)
dir.create(hdir, recursive = TRUE, showWarnings = FALSE)

# --- clustering + claudin-low (pam50) ---
res <- centroid_clustering(int_medctr, prefix = "intrinsic_expr_medCtr",
                           path = hdir, write_files = TRUE)
clow <- identify_claudin_low(res, int_medctr)
message(sprintf("claudin-low: %d samples", length(clow)))

# --- assemble call table ---
samp_ids <- colnames(int_medctr)
meta <- data.frame(sample = samp_ids, stringsAsFactors = FALSE)
if (!is.null(subtype_col) && !is.na(subtype_col) && subtype_col %in% colnames(samples)) {
  meta$clinical_subtype <- samples[samp_ids, subtype_col]
}
meta$PAM50_Call <- scores[samp_ids, "Call"]
meta$PAM50_Call_CLOW <- ifelse(samp_ids %in% clow, "CLOW", meta$PAM50_Call)

# --- silhouette widths (clusters = CLOW-adjusted calls; 1 - Pearson distance) ---
clusters <- as.numeric(as.factor(meta$PAM50_Call_CLOW))
names(clusters) <- samp_ids
dist_matrix <- as.dist(1 - cor(int_expr, method = "pearson"))
sil <- silhouette(clusters, dist_matrix)
meta$silhouette_width <- sil[, 3]
meta$sil_width_bin <- ifelse(meta$silhouette_width < 0, "Negative", "Positive")

write.table(meta, snakemake@output[["calls_clow"]],
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- silhouette plot ---
pdf(snakemake@output[["silhouette"]], width = 9, height = 8)
plot(sil, border = NA, col = c("red", "gold", "pink", "blue", "lightblue", "green")[seq_len(max(clusters))])
abline(v = 0)
dev.off()

# --- intrinsic heatmap ---
ok_hm <- requireNamespace("ComplexHeatmap", quietly = TRUE) &&
  requireNamespace("circlize", quietly = TRUE)
if (ok_hm) {
  # Column/row dendrograms from the Cluster 3.0 .atr/.gtr (ctc), else let
  # ComplexHeatmap cluster internally (visual only; does not affect calls).
  col_dend <- TRUE; row_dend <- TRUE
  if (requireNamespace("ctc", quietly = TRUE)) {
    col_dend <- tryCatch(
      ctc::xcluster2r(file = file.path(hdir, "intrinsic_expr_medCtr.atr"),
                      distance = "pearson", labels = FALSE, fast = FALSE, clean = TRUE, dec = "."),
      error = function(e) { message("xcluster2r(atr) failed; internal clustering"); TRUE })
    row_dend <- tryCatch(
      ctc::xcluster2r(file = file.path(hdir, "intrinsic_expr_medCtr.gtr"),
                      distance = "pearson", labels = FALSE, fast = FALSE, clean = TRUE, dec = "."),
      error = function(e) { message("xcluster2r(gtr) failed; internal clustering"); TRUE })
  }

  ann_list <- list("PAM50_zhao4+CLOW" = meta$PAM50_Call_CLOW,
                   "Silhouette_width" = meta$sil_width_bin)
  if (!is.null(meta$clinical_subtype)) {
    ann_list <- c(list("Clinical_subtype" = meta$clinical_subtype), ann_list)
  }
  column_ha <- ComplexHeatmap::HeatmapAnnotation(df = as.data.frame(ann_list, check.names = FALSE))

  hm <- ComplexHeatmap::Heatmap(
    int_medctr, name = "Gene exp.",
    show_column_names = FALSE, show_row_names = FALSE,
    top_annotation = column_ha,
    row_title = "genes", column_title = "samples",
    cluster_rows = row_dend, cluster_columns = col_dend,
    col = circlize::colorRamp2(c(-4, 0, 4), c("blue", "black", "yellow")),
    use_raster = TRUE
  )
  pdf(heatmap_out, width = 11, height = 8)
  ComplexHeatmap::draw(hm, column_title = sprintf("Intrinsic subtype heatmap (n=%d)", nrow(meta)))
  dev.off()
} else {
  message("ComplexHeatmap/circlize unavailable; writing placeholder heatmap")
  pdf(heatmap_out); plot.new(); text(0.5, 0.5, "heatmap unavailable"); dev.off()
}
message("intrinsic_clustering: done")

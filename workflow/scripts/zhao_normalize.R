# Snakemake shim: Zhao4 subtype-specific quantile centering of the PAM50 matrix.

suppressPackageStartupMessages(library(pam50perou))

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message"); sink(log, type = "output")

cfg <- snakemake@config
id_col <- cfg$inputs$sample_id_column
subtype_col <- cfg$zhao_normalization$subtype_column
if (is.null(subtype_col) || is.na(subtype_col) || !nzchar(subtype_col)) {
  stop("zhao_normalization.enabled is true but subtype_column is not set")
}

expr <- as.matrix(read.table(snakemake@input[["pam50_expr"]],
                             sep = "\t", header = TRUE, row.names = 1, check.names = FALSE))
samples <- read.table(snakemake@input[["samples_filtered"]],
                      sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
if (!subtype_col %in% colnames(samples)) {
  stop(sprintf("subtype_column '%s' not in sample sheet", subtype_col))
}

subtypes <- stats::setNames(as.character(samples[[subtype_col]]), samples[[id_col]])
subtypes <- subtypes[colnames(expr)]
message(sprintf("Zhao groups: %s",
                paste(names(table(subtypes)), table(subtypes), sep = "=", collapse = ", ")))

z <- zhao_center(expr, subtypes)

out <- snakemake@output[[1]]
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
write.table(data.frame(Gene = rownames(z), z, check.names = FALSE),
            out, sep = "\t", row.names = FALSE, quote = FALSE)
message("zhao_normalize: done")

# Snakemake shim: PAM50 classification. Writes the scores table in the same
# format as the reference pipeline (*_pam50scores.txt: tab-separated, empty
# top-left header cell, sample IDs as row names).

suppressPackageStartupMessages(library(pam50))

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message"); sink(log, type = "output")

cfg <- snakemake@config$classifier

expr <- as.matrix(read.table(snakemake@input[["expr"]],
                             sep = "\t", header = TRUE, row.names = 1, check.names = FALSE))

scores <- predict_pam50(
  expr,
  calibration  = cfg$calibration,
  collapse     = cfg$collapse_method,
  has_clinical = cfg$has_clinical
)

out <- snakemake@output[["scores"]]
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
# col.names = NA reproduces the reference layout (leading empty header cell).
write.table(scores, out, sep = "\t", col.names = NA, quote = TRUE)
message(sprintf("classify: %d samples; calls: %s", nrow(scores),
                paste(names(table(scores$Call)), table(scores$Call),
                      sep = "=", collapse = " ")))

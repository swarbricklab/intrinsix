# Snakemake shim: enrich the raw PAM50 scores with the sample-sheet metadata and
# write the final, keyed call table.

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message"); sink(log, type = "output")

cfg <- snakemake@config
id_col <- cfg$inputs$sample_id_column

scores <- read.table(snakemake@input[["scores"]],
                     sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
samples <- read.table(snakemake@input[["samples_filtered"]],
                      sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# Join sample-sheet columns (minus the id) onto the scores by sample ID.
meta_cols <- setdiff(colnames(samples), id_col)
meta <- samples[match(rownames(scores), samples[[id_col]]), meta_cols, drop = FALSE]

out_df <- data.frame(
  sample = rownames(scores),
  meta,
  scores,
  check.names = FALSE,
  row.names = NULL
)
colnames(out_df)[1] <- id_col

out <- snakemake@output[[1]]
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
write.table(out_df, out, sep = "\t", row.names = FALSE, quote = FALSE)
message(sprintf("finalize: wrote %d calls with %d metadata columns",
                nrow(out_df), length(meta_cols)))

# intrinsix (development version)

## Cleaner, well-formed `pam50_run` output

* `pam50_run()` / `intrinsix run` now writes a **single** well-formed call table,
  `pam50_calls.tsv`: one row per sample, a named id column, the carried-through
  sample metadata, then the classifier columns under clean, unique, snake_case
  names (`call`, `confidence`, `score_basal`..`score_normal`, `ror_s`,
  `ror_s_group`, `proliferation_score`, `ror_p`, `ror_p_group`, `er_signature`,
  `her2_signature`; plus `call_clow`, `silhouette_width`, `sil_width_bin` with
  clustering). This removes the duplicate `Her2` column of the old raw scores
  file and folds in the former `pam50_calls_clow.tsv`. The classifier is
  unchanged — calls and values are identical; only the serialization is cleaner.
  The in-memory `pam50_result` keeps the raw Perou-native column names.
* The legacy Perou-format `"<prefix>_pam50scores.txt"` (quoted cells, unnamed
  first column, duplicate `Her2`) is **no longer written by default**. Set
  `classifier$legacy_perou_scores: true` to emit it for downstream Perou Lab
  scripts that expect that exact format.
* Sample metadata is carried from the input's `obs` / `colData` (or the
  `inputs$samples` sheet) into `pam50_calls.tsv`. New `inputs$metadata_columns`
  allowlist keeps only the named columns (default: all). A carried metadata
  column whose name collides with a classifier output column is now an **error**
  (rename it upstream, or drop it via `inputs$metadata_columns`).

# intrinsix 0.1.0

First packaged release — an installable R package (`intrinsix`) with a config-driven
`intrinsix run` pipeline.

## Inputs

* Accepts a matrix, a `SummarizedExperiment` / `SingleCellExperiment`, an
  `anndataR` `AnnData`, or an `.h5ad` path (via `anndataR`, R >= 4.5); the sample
  sheet is taken from `colData` / `obs`, and `layer=` picks the assay/layer.
  `SummarizedExperiment` / `anndataR` are optional (`Suggests`).

## Algorithm (verified byte-identical to the production bcarc_wts pipeline)

* `predict_pam50()` — PAM50 classification (Spearman correlation to the Perou
  centroids) with ROR-S / ROR-P risk scores.
* `zhao_center()` — Zhao4 subtype-specific quantile centering.
* `map_to_pam50()`, `median_center()` — input preparation.
* `centroid_clustering()`, `identify_claudin_low()` — intrinsic Cluster 3.0
  replication and claudin-low identification.
* Bundled reference data: `pam50_centroids()`, `pam50_annotation()`,
  `sigma_table()`, `intrinsic_genes()`, `training_set()`.

## Interactive interface

* `pam50_subtype()` — one call from an expression matrix to a `pam50_result`
  (map → optional Zhao → classify → optional intrinsic clustering).
* `pam50_result` S3 class with `print()`, `summary()`, `plot()`, `as.data.frame()`.
* `pam50_add_samples()` — incremental calling: add new samples to an
  already-called cohort one at a time (Perou-recommended), keeping prior calls
  fixed.

## Visualisations (ported from the Perou Lab / Alec Lobanov scripts)

* `plot_pam50_scores()` — per-subtype correlation + proliferation boxplots.
* `plot_pam50_pca()` — training-vs-test PCA.
* `pam50_heatmap()` — intrinsic-subtype heatmap (shared with the workflow).

## Command line

* `intrinsix` CLI (`classify` / `zhao` / `cluster`) over the same functions.

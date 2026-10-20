# pam50 0.1.0

First packaged release. The repo is now both an installable R package (`pam50`)
and a Snakemake workflow.

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

## Visualisations (ported from the Perou Lab / Alec Lobanov scripts)

* `plot_pam50_scores()` — per-subtype correlation + proliferation boxplots.
* `plot_pam50_pca()` — training-vs-test PCA.
* `pam50_heatmap()` — intrinsic-subtype heatmap (shared with the workflow).

## Command line

* `pam50` CLI (`classify` / `zhao` / `cluster`) over the same functions.

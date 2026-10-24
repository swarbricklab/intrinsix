# intrinsix

[![CI](https://github.com/swarbricklab/intrinsix/actions/workflows/ci.yml/badge.svg)](https://github.com/swarbricklab/intrinsix/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

PAM50 intrinsic molecular subtyping of breast-cancer expression data — an **R
package** (`intrinsix`) for interactive / scripted use, with a config-driven
`intrinsix run` pipeline that plugs into a DVC super-pipeline.

It packages the exact Perou Lab PAM50 classifier (published centroids, Spearman
rank correlation) with optional Zhao4 subtype-specific normalisation, ROR-S/ROR-P
risk scores, intrinsic clustering and claudin-low identification. The algorithm
and all reference data are bundled in the package, so there is no git-submodule
checkout or external data pull at runtime.

Extracted from the previously duplicated `modules/PAM50/` stages in `bcarc_wts`
and `brca_mega_atlas-data`, which it is intended to replace. Verified
byte-identical to that pipeline at unit, golden, and full end-to-end levels.

## Provenance & genealogy

**This package is a faithful repackaging of Perou Lab research code — not a
reimplementation.** The scoring mathematics is unchanged (byte-for-byte), and most
of the R under `R/` is ported verbatim or transcribed from existing scripts. What
is new here is packaging: an installable/tested R package, a config-driven runner,
annotated-object (SE/h5ad) inputs, incremental calling, and documentation. The
scientific credit for the method and its implementation belongs upstream.

The PAM50 subtype predictor is from **Parker & Perou et al.**, *J Clin Oncol*
2009. The Garvan/Swarbrick driver, Zhao4 normalisation and intrinsic
clustering / claudin-low scripts are by **Alec Lobanov** (Perou Lab), with the
PAM50 assay protocol from **Arantxa Fernandez-Martinez** (Perou Lab).

The upstream sources are kept under [`vendor/`](vendor/) for audit (never sourced
at runtime, and `.Rbuildignore`d so not part of the built package): the two
perou-lab repositories are pinned **git submodules** (exact, tamper-evident
upstream commits), and Alec's scripts — the lab's own, with no standalone public
upstream — are kept as unmodified copies. To populate the submodules in a fresh
checkout (optional; only to inspect the upstream code):

```bash
git submodule update --init vendor/perou-lab/*
```

Where each piece of `R/` comes from:

| `intrinsix` source | Derived from | Origin |
|---|---|---|
| `R/perou_internal.R` | `subtypePrediction_functions.R` (**verbatim**) | Parker & Perou (Perou Lab) |
| `R/predict.R` | `subtypePrediction_distributed.R` (**transcribed**) | Parker & Perou (Perou Lab) |
| `R/perou_plots.R` | `subtypePrediction_functions.R` + intrinsic-clustering script | Parker & Perou; Lobanov |
| `R/zhao.R` | `PAM50_zhao_normalization.R` | Lobanov (Perou Lab) |
| `R/clustering.R` | `clustering_scripts.R`, `cluster3_replicate.R`, `identify_claudin_low.R` | Lobanov (Perou Lab) |
| `inst/extdata/` (centroids, SIGMA, training set, annotation, intrinsic gene list) | PAM50 published reference data | Perou Lab |
| `R/prepare.R` | logic of `prepare_pam50_input.py` | Swarbrick Lab |
| `R/data.R`, `input.R`, `subtype.R`, `run.R`, `add_samples.R`, `cli.R` | new packaging / interface code | Swarbrick Lab (Reeves) |

See [Attribution & citation](#attribution--citation) for how to cite.

## Input

The expression can be a **matrix**, a **SummarizedExperiment / SingleCellExperiment**,
an **AnnData**, or a path to an **`.h5ad`** file (see [`DESIGN.md`](DESIGN.md#input-contract)).
For a plain matrix you also supply a sample sheet; for the annotated objects the
sample sheet is taken from them (`colData` / `obs`):

- **expression** — genes × samples, HGNC symbols, log2-scale normalised values
  (pick the layer/assay with `layer =`, e.g. `"log2_uqn"`).
- **samples** — one row per sample; an ID column plus any columns used for
  filtering, subtype grouping, or passthrough metadata.

By default PAM50 is called on **all samples**; set `filter.column` to restrict
(e.g. to an eligibility flag).

> `.h5ad` support needs the `anndataR` package (R ≥ 4.5); `SummarizedExperiment`
> input needs that package. Both are optional (`Suggests`) — the matrix path needs
> neither. The conda `environment.yml` and container include them.

## Use as an R package

```r
# install.packages("remotes")
remotes::install_github("swarbricklab/intrinsix")
library(intrinsix)

res   <- pam50_subtype(expr, samples = meta, subtype_col = "clinical_subtype")
res   <- pam50_subtype(sce, layer = "log2_uqn")      # SummarizedExperiment / SCE
res   <- pam50_subtype("cohort.h5ad", layer = "log2_uqn")   # .h5ad path
calls <- predict_pam50(expr)            # low-level: expr = PAM50 pcrID rows
```

Reference data is bundled: `pam50_centroids()`, `pam50_annotation()`,
`sigma_table()`, `intrinsic_genes()`.

## Use from the command line

The same functions, via the `intrinsix` CLI (installed with the package; available as
`intrinsix` inside the container, or `Rscript -e 'intrinsix::cli()' --args ...`):

```bash
intrinsix classify --expr expr.tsv --out scores.tsv [--map-symbols]
intrinsix zhao     --expr expr.tsv --samples samples.tsv --subtype-col clinical_subtype --out zhao.tsv
intrinsix cluster  --expr-medctr intrinsic_medctr.tsv --outdir out/ --out-claudin claudin_low.txt
```

## Run the whole pipeline (config-driven)

`intrinsix run` reads a config (filter/dedup → subtype → optional clustering) and
writes the outputs — the command a DVC stage invokes. A template config ships with
the package (`system.file("extdata", "config.example.yaml", package = "intrinsix")`):

```bash
intrinsix run --config config.yaml           # or: Rscript -e 'intrinsix::pam50_run("config.yaml")'
```

It writes a single well-formed call table, `pam50_calls.tsv` (one row per sample:
a named id column, the carried-through sample metadata, then clean snake_case
classifier columns — `call`, `confidence`, `score_basal`..`score_normal`, the
`ror_*` scores, `er_signature`, `her2_signature`, and `call_clow` /
`silhouette_width` / `sil_width_bin` when clustering is on). Sample metadata is
inherited from the input's `obs`/`colData` (narrow it with
`inputs$metadata_columns`); a metadata column that collides with a classifier
column is an error. Set `classifier.legacy_perou_scores: true` to also emit the
legacy Perou-format `"<prefix>_pam50scores.txt"` for scripts that expect it.

Add `--report` for an **opt-in, self-contained HTML report** — the run's
visualisations (per-subtype correlation boxplots, training-vs-test PCA, and, when
clustering is on, the intrinsic heatmap + silhouette) plus the per-sample results
table, with explanatory prose for a reader new to the run. It is written to
`<output_dir>/intrinsix_report.html` (override with `--report-html PATH`, or set
`outs.report: true` in the config). The report needs the `rmarkdown` + `knitr`
packages (`Suggests`); default behaviour is unchanged when the flag is off.

```bash
intrinsix run --config config.yaml --report                 # -> <output_dir>/intrinsix_report.html
intrinsix run --config config.yaml --report-html report.html
# or from R: pam50_run("config.yaml", report = TRUE)
```

### Running environments

The package is pure R, so pick whatever the platform has — **no singularity
required**:

| Context | Environment | Invocation |
|---|---|---|
| Laptop / interactive | native install (`remotes` + `BiocManager`) | `intrinsix run --config config.yaml` |
| Reproducible, any OS | **conda** (`environment.yml`) | `conda run -n intrinsix intrinsix run --config config.yaml` |
| HPC, locked-down | singularity image | `module load singularity && singularity exec intrinsix.sif intrinsix run --config config.yaml` |

```bash
# conda env (portable default)
conda env create -f environment.yml && conda activate intrinsix
R -e 'remotes::install_local(".")'
```

### DVC integration

Each parent DVC pipeline keeps a `intrinsix` stage whose `cmd` is one of the above —
DVC does the stage-level caching; no second workflow engine.

## Layout

The repo root **is** the R package.

```
DESCRIPTION, NAMESPACE, R/, inst/   the intrinsix R package (algorithm, data, CLI, runner)
tests/                              package tests (testthat)
vignettes/                          'Getting started' vignette
containers/                         Dockerfile + env for the optional HPC image
environment.yml                     conda env (portable default)
inst/extdata/config.example.yaml    config template for `intrinsix run`
vendor/                             upstream provenance (perou-lab submodules + copies)
```

## Attribution & citation

This package recycles Perou Lab research code heavily (see
[Provenance & genealogy](#provenance--genealogy) for the per-file breakdown).
Please credit the originators, not just this package:

- **PAM50 method:** Parker JS, Perou CM, *et al.* "Supervised risk predictor of
  breast cancer based on intrinsic subtypes." *J Clin Oncol* 2009;27(8):1160–1167.
  [doi:10.1200/JCO.2008.18.1370](https://doi.org/10.1200/JCO.2008.18.1370)
- **Classifier R implementation:** Joel S. Parker & Charles M. Perou (Perou Lab).
- **Driver / Zhao4 / intrinsic clustering / claudin-low scripts:** Alec Lobanov
  (Perou Lab); PAM50 assay protocol: Arantxa Fernandez-Martinez (Perou Lab).
- **Packaging & refactor:** John Reeves (Swarbrick Lab, Garvan Institute).

`citation("intrinsix")` lists the software and the method paper. See also
[`CITATION.cff`](CITATION.cff), [`inst/CITATION`](inst/CITATION), and the
provenance notes under [`vendor/`](vendor/).

## Licence

Workflow code: MIT (see [`LICENSE`](LICENSE)). The vendored upstream `perou-lab`
code retains its own licensing — to be confirmed before public release.

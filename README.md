# pam50

[![CI](https://github.com/swarbricklab/pam50/actions/workflows/ci.yml/badge.svg)](https://github.com/swarbricklab/pam50/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

PAM50 intrinsic molecular subtyping of breast-cancer expression data — an **R
package** (`pam50`) for interactive / scripted use, with a config-driven
`pam50 run` pipeline that plugs into a DVC super-pipeline.

It packages the exact Perou Lab PAM50 classifier (published centroids, Spearman
rank correlation) with optional Zhao4 subtype-specific normalisation, ROR-S/ROR-P
risk scores, intrinsic clustering and claudin-low identification. The algorithm
and all reference data are bundled in the package, so there is no git-submodule
checkout or external data pull at runtime.

Extracted from the previously duplicated `modules/PAM50/` stages in `bcarc_wts`
and `brca_mega_atlas-data`, which it is intended to replace. Verified
byte-identical to that pipeline at unit, golden, and full end-to-end levels.

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
remotes::install_github("swarbricklab/pam50")
library(pam50)

res   <- pam50_subtype(expr, samples = meta, subtype_col = "clinical_subtype")
res   <- pam50_subtype(sce, layer = "log2_uqn")      # SummarizedExperiment / SCE
res   <- pam50_subtype("cohort.h5ad", layer = "log2_uqn")   # .h5ad path
calls <- predict_pam50(expr)            # low-level: expr = PAM50 pcrID rows
```

Reference data is bundled: `pam50_centroids()`, `pam50_annotation()`,
`sigma_table()`, `intrinsic_genes()`.

## Use from the command line

The same functions, via the `pam50` CLI (installed with the package; available as
`pam50` inside the container, or `Rscript -e 'pam50::cli()' --args ...`):

```bash
pam50 classify --expr expr.tsv --out scores.tsv [--map-symbols]
pam50 zhao     --expr expr.tsv --samples samples.tsv --subtype-col clinical_subtype --out zhao.tsv
pam50 cluster  --expr-medctr intrinsic_medctr.tsv --outdir out/ --out-claudin claudin_low.txt
```

## Run the whole pipeline (config-driven)

`pam50 run` reads a config (filter/dedup → subtype → optional clustering) and
writes the outputs — the command a DVC stage invokes. A template config ships with
the package (`system.file("extdata", "config.example.yaml", package = "pam50")`):

```bash
pam50 run --config config.yaml           # or: Rscript -e 'pam50::pam50_run("config.yaml")'
```

### Running environments

The package is pure R, so pick whatever the platform has — **no singularity
required**:

| Context | Environment | Invocation |
|---|---|---|
| Laptop / interactive | native install (`remotes` + `BiocManager`) | `pam50 run --config config.yaml` |
| Reproducible, any OS | **conda** (`environment.yml`) | `conda run -n pam50 pam50 run --config config.yaml` |
| HPC, locked-down | singularity image | `module load singularity && singularity exec pam50.sif pam50 run --config config.yaml` |

```bash
# conda env (portable default)
conda env create -f environment.yml && conda activate pam50
R -e 'remotes::install_local(".")'
```

### DVC integration

Each parent DVC pipeline keeps a `pam50` stage whose `cmd` is one of the above —
DVC does the stage-level caching; no second workflow engine.

## Layout

The repo root **is** the R package.

```
DESCRIPTION, NAMESPACE, R/, inst/   the pam50 R package (algorithm, data, CLI, runner)
tests/                              package tests (testthat)
vignettes/                          'Getting started' vignette
containers/                         Dockerfile + env for the optional HPC image
environment.yml                     conda env (portable default)
inst/extdata/config.example.yaml    config template for `pam50 run`
vendor/                             unmodified upstream snapshots (provenance)
```

## Attribution & citation

Core algorithm and the Garvan driver / Zhao4 / clustering scripts:
**Alec Lobanov** (Perou Lab); PAM50 protocol via Arantxa Fernandez-Martinez
(Perou Lab). PAM50: Parker et al., *J Clin Oncol* 2009. See
[`CITATION.cff`](CITATION.cff) and `vendor/perou-lab/PROVENANCE.txt`.

## Licence

Workflow code: MIT (see [`LICENSE`](LICENSE)). The vendored upstream `perou-lab`
code retains its own licensing — to be confirmed before public release.

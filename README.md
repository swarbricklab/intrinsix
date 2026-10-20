# pam50

[![CI](https://github.com/swarbricklab/pam50/actions/workflows/ci.yml/badge.svg)](https://github.com/swarbricklab/pam50/actions/workflows/ci.yml)
[![Snakemake](https://img.shields.io/badge/snakemake-%E2%89%A57.32.4-brightgreen.svg)](https://snakemake.github.io)
[![Code style: snakefmt](https://img.shields.io/badge/code%20style-snakefmt-000000.svg)](https://github.com/snakemake/snakefmt)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

PAM50 intrinsic molecular subtyping of breast-cancer expression data — as **both
an R package** (`pam50`, for interactive / scripted use) **and a portable
Snakemake workflow** (to plug into a DVC super-pipeline). One repo, one source of
truth.

It packages the exact Perou Lab PAM50 classifier (published centroids, Spearman
rank correlation) with optional Zhao4 subtype-specific normalisation, ROR-S/ROR-P
risk scores, intrinsic clustering and claudin-low identification. The algorithm
and all reference data are bundled in the package, so there is no git-submodule
checkout or external data pull at runtime.

Extracted from the previously duplicated `modules/PAM50/` stages in `bcarc_wts`
and `brca_mega_atlas-data`, which it is intended to replace. Verified
byte-identical to that pipeline at unit, golden, and full end-to-end levels.

## Input

Two tables (see [`DESIGN.md`](DESIGN.md#input-contract)):

- **expression** — genes × samples, HGNC symbols, log2-scale normalised values.
- **samples** — one row per sample; an ID column plus any columns used for
  filtering, subtype grouping, or passthrough metadata.

By default PAM50 is called on **all samples**; set `filter.column` to restrict
(e.g. to an eligibility flag). Users starting from an `.h5ad` can use the optional
Python adapter rule.

## Use as an R package

```r
# install.packages("remotes")
remotes::install_github("swarbricklab/pam50")
library(pam50)

calls <- predict_pam50(expr)            # expr: genes x samples, PAM50 pcrID rows
z     <- zhao_center(expr, subtypes)    # Zhao4 subtype-specific centering
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

## Use as a Snakemake workflow

```bash
snakemake --use-singularity --configfile config/config.yaml --cores 4
```

Each parent DVC pipeline keeps a `pam50` stage whose `cmd` shells out to
`snakemake`; the parent's preprocessing emits `expression.tsv` + `samples.tsv`.

## Layout

The repo root **is** the R package (`DESCRIPTION`, `R/`, `inst/`); the Snakemake
workflow sits alongside it (`.Rbuildignore`d from the package build).

```
DESCRIPTION, NAMESPACE, R/, inst/   the pam50 R package (algorithm, data, CLI)
tests/                              package tests (testthat)
workflow/Snakefile                  DAG
workflow/rules/pam50.smk            steps (prepare -> zhao -> classify -> finalize [-> cluster])
workflow/scripts/                   thin R wrappers (library(pam50) -> call it)
workflow/envs/                      conda fallback for platforms without singularity
containers/                         Dockerfile + env for the pam50 image
config/                             config.yaml (the input contract) + schema
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

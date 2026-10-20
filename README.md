# pam50

[![CI](https://github.com/swarbricklab/pam50/actions/workflows/ci.yml/badge.svg)](https://github.com/swarbricklab/pam50/actions/workflows/ci.yml)
[![Snakemake](https://img.shields.io/badge/snakemake-%E2%89%A57.32.4-brightgreen.svg)](https://snakemake.github.io)
[![Code style: snakefmt](https://img.shields.io/badge/code%20style-snakefmt-000000.svg)](https://github.com/snakemake/snakefmt)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A standalone, portable **Snakemake** workflow for PAM50 intrinsic molecular
subtyping of breast-cancer expression data.

> **Status: scaffold / work in progress.** The structure, input contract,
> vendored reference data and design are in place. The algorithm packaging
> (`pam50perou`), the thin rule scripts, containers and tests are the next
> milestone — see [`DESIGN.md`](DESIGN.md).

It packages the exact Perou Lab PAM50 classifier (published centroids, Spearman
rank correlation) with optional Zhao4 subtype-specific normalisation, ROR-S/ROR-P
risk scores, intrinsic clustering and claudin-low identification. The algorithm
and all reference data are bundled in the `pam50perou` R package, so there is no
git-submodule checkout or external data pull at runtime.

This workflow is extracted from the previously duplicated `modules/PAM50/` stages
in `bcarc_wts` and `brca_mega_atlas-data`, which it is intended to replace.

## Input

Two tables (see [`DESIGN.md`](DESIGN.md#input-contract)):

- **expression** — genes × samples, HGNC symbols, log2-scale normalised values.
- **samples** — one row per sample; an ID column plus any columns used for
  filtering, subtype grouping, or passthrough metadata.

By default PAM50 is called on **all samples**; set `filter.column` to restrict
(e.g. to an eligibility flag). Users starting from an `.h5ad` can use the optional
Python adapter rule.

## Quick start (target)

```bash
snakemake --use-singularity --configfile config/config.yaml --cores 4
```

## Layout

```
workflow/Snakefile          DAG
workflow/rules/pam50.smk     steps (prepare -> zhao -> classify -> finalize [-> cluster])
workflow/scripts/            thin R wrappers around pam50perou
workflow/pkg/pam50perou/     the algorithm + bundled reference data (R package)
workflow/envs/               conda fallback for platforms without singularity
config/                      config.yaml (the input contract) + schema
vendor/                      unmodified upstream snapshots (provenance)
.tests/                      integration test data + expected calls
```

## Attribution & citation

Core algorithm and the Garvan driver / Zhao4 / clustering scripts:
**Alec Lobanov** (Perou Lab); PAM50 protocol via Arantxa Fernandez-Martinez
(Perou Lab). PAM50: Parker et al., *J Clin Oncol* 2009. See
[`CITATION.cff`](CITATION.cff) and `vendor/perou-lab/PROVENANCE.txt`.

## Licence

Workflow code: MIT (see [`LICENSE`](LICENSE)). The vendored upstream `perou-lab`
code retains its own licensing — to be confirmed before public release.

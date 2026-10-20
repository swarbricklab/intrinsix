# PAM50 workflow — design

Status: **implemented + verified**. Algorithm ported into `pam50perou` (classifier,
Zhao4, intrinsic clustering, claudin-low), rule wrappers written, container built
(`ghcr.io/swarbricklab/pam50`), CI + golden tests in place, and a `pam50` CLI.
Verified byte-identical to the production bcarc_wts pipeline at unit, golden, and
full end-to-end (containerised, 225 samples) levels. Remaining: switch the two DVC
parents over; confirm the upstream perou-lab licence before going public.

A standalone, portable Snakemake workflow for PAM50 intrinsic molecular
subtyping of breast-cancer expression data, extracted from the duplicated
`modules/PAM50/` stages in `bcarc_wts` and `brca_mega_atlas-data`.

## Why this exists

PAM50 currently lives as copy-adapted modules in two DVC repos. Problems the
extraction fixes:

1. **Duplication** — one source of truth instead of two drifting copies.
2. **Script-sourcing fragility** — the Perou algorithm was pulled in as two live
   git submodules and `source()`d by path, relying on global-variable side
   effects. In `brca_mega_atlas-data` the submodules weren't even checked out, so
   the stage was silently broken. Here the algorithm becomes an **installed R
   package** (`pam50perou`) with reference data bundled as `inst/extdata`.
3. **Hidden cross-repo input coupling** — eligibility (`is_pam50_eligible`) and
   `clinical_subtype` came from the `Swarbricklab/metadata` repo via `dt`. Those
   are now **declared inputs**, not reached across repos.

Scope preserved deliberately: this packages **the exact Perou centroids + Zhao4
normalisation** the lab already uses. Swapping to `genefu`/`sspbc` would change
the calls — that is a separate, results-changing (goal-2) decision.

## Architecture

```
pam50perou/   R package: the algorithm + bundled reference data
  R/                        refactored functions (exported API below)
  inst/extdata/             centroids, annotation, calibration medians,
                            training set, SIGMA, intrinsic gene list
vendor/                     UNMODIFIED upstream snapshots (provenance only)
workflow/Snakefile          thin DAG
workflow/rules/pam50.smk    one rule per step, container: directives
workflow/scripts/           ~5-line wrappers: read -> call pam50perou -> write
config/config.yaml          the input contract (below)
```

Two containers (mixing Py/R is fine):

- `ghcr.io/swarbricklab/pam50` — R + `pam50perou` + ComplexHeatmap/ctc/impute.
  Runs every core step. The workflow is **R-only on the TSV entry path**.
- an existing GHCR Python/anndata image — used *only* by the optional `h5ad`
  input adapter, for users who start from AnnData rather than a TSV.

## Input contract

The standalone entry point is **two tables** (TSV/CSV):

- `expression` — genes × samples. Rows = HGNC gene symbols, columns = sample IDs,
  values = **log2-scale normalised** expression (log2 UQN, log2(TPM+1), …). The
  workflow maps symbols → the classifier's `pcrID` via the bundled
  `pam50_annotation.txt` and subsets the 50 PAM50 genes + the intrinsic list.
- `samples` — one row per sample; an ID column matching the expression columns,
  plus any columns used for filtering / subtype grouping / passthrough metadata.

Users who have an `.h5ad` can use the optional adapter rule
(`h5ad_to_expression`, Python container) to emit `expression`/`samples` first;
this is where the parent repos' preprocessing would plug in.

### Generalised sample filter (new vs. the hard-coded eligibility flag)

```yaml
filter:
  column: null          # column in `samples` to filter on; null => ALL samples
  include: []           # keep rows whose value is in this set (used if column set)
```

Default (`column: null`) calls PAM50 on **all samples** — the common case for
new users. The old bcarc_wts behaviour is just
`filter: {column: is_pam50_eligible, include: ["TRUE"]}`.

### Generalised replicate dedup (optional)

```yaml
deduplicate:
  by: null              # group column (e.g. tissue_id); null => no dedup
  keep: first           # lowest sample id lexically within each group
  selections: {}        # explicit overrides: "<group value>": "<sample id>"
```

### Zhao4 subtype-specific normalisation (optional)

```yaml
zhao_normalization:
  enabled: true
  subtype_column: clinical_subtype   # REQUIRED if enabled
  # values must match SIGMA columns:
  # ERpos_HER2neg, ERneg_HER2neg, HER2pos_ERpos, HER2pos_ERneg
```

When disabled (e.g. no clinical subtype available — the "all samples" case), the
classifier's own median-centering is used instead.

### Classifier + clustering

```yaml
classifier:
  calibration: -1          # Perou calibrationParameters (-1 = none)
  collapse_method: mean    # mean | iqr
  has_clinical: false      # ROR-C needs a T-stage column; keep false unless provided
  prefix: PAM50scores
intrinsic_clustering:
  enabled: true            # claudin-low identification + silhouette + heatmap
```

## Exported package API (target)

```r
pam50perou::pam50_annotation()      # bundled reference-data accessors
pam50perou::pam50_centroids()
pam50perou::sigma_table()
pam50perou::intrinsic_genes()

map_to_pam50(expr, annotation)                 # symbol -> pcrID, subset 50, dedup
zhao_center(expr, subtypes, sigma)             # Zhao4 quantile centering
predict_pam50(expr, calibration=-1,            # wraps sspPredict(); returns the
              collapse="mean", has_clinical=FALSE)   #   scores/Call table + plots
centroid_clustering(expr_medctr, outdir)       # Cluster3.0 replication
identify_claudin_low(cluster_files, expr_medctr)
```

Internal (not exported): the Perou helpers `readarray`, `standardize`,
`medianCtr`, `sspPredict`, `overlapSets`, `myHeatmap`, `myplot`, `pcaEA` from
`subtypePrediction_functions.R`, refactored to take explicit args instead of
reading globals.

## Acceptance test (goal-1 insurance)

A golden fixture (small expression matrix) must reproduce **call-identical**
output to the current `bcarc_wts` module before either parent is switched over.
That is the gate for retiring the old `modules/PAM50/`.

## Consumption by the parent DVC repos

Each parent keeps a single `pam50` DVC stage whose `cmd` shells out to
`snakemake` (the lab already nests Snakemake under DVC). The parent's
preprocessing emits `expression.tsv` + `samples.tsv`; DVC tracks deps/outs as
today.

## Open questions

- **Licensing of the vendored Perou code.** The workflow code is MIT. The
  upstream `perou-lab` scripts carry their own terms — confirm the upstream
  licence before the repo goes public. (`vendor/perou-lab/` keeps the raw
  snapshots + `PROVENANCE.txt` with exact commits.)
- **Public test data.** Need a small, redistributable expression matrix with
  known PAM50 calls for the acceptance test + CI (candidate: a public GEO/TCGA
  subset, or the upstream `PAM50_example/sampleInputFile_200subjects.txt`).
- **Package home.** Vendored inside this repo for v1 (one repo, as requested).
  Split into its own `swarbricklab/pam50perou` repo later if a second consumer
  appears.

## Attribution

Core algorithm and the Garvan driver/Zhao4/clustering scripts: **Alec Lobanov**
(Perou Lab); PAM50 protocol via Arantxa Fernandez-Martinez (Perou Lab). PAM50 is
Parker et al., JCO 2009. See `vendor/perou-lab/PROVENANCE.txt`.

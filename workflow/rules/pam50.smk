# PAM50 rules. Each step is a thin wrapper that reads TSVs, calls the pam50perou
# package, and writes TSVs. The algorithm + reference data live in the package
# (workflow/pkg/pam50perou), baked into the `pam50` container.
#
# Paths use the out() helper (os.path.join) from the Snakefile, NOT f-strings:
# f-strings in rule directives interact badly with Snakemake's wildcard handling.

PAM50_IMG = config["containers"]["pam50"]
SCORES = out(PREFIX + "_pam50scores.txt")


# Step 1 -- filter samples, map symbols -> pcrID, subset PAM50 + intrinsic genes,
# median-centre the intrinsic matrix.
rule prepare_input:
    input:
        expression=config["inputs"]["expression"],
        samples=config["inputs"]["samples"],
    output:
        pam50_expr=out("pam50_expr.tsv"),
        samples_filtered=out("samples_filtered.tsv"),
        int_expr=out("intrinsic_heatmap", "int_expr.tsv"),
        int_expr_medctr=out("intrinsic_heatmap", "int_expr_medCtr.tsv"),
    log:
        out("logs", "prepare_input.log"),
    container:
        PAM50_IMG
    # conda:  # conda fallback for platforms without singularity
    #     "../envs/pam50.yaml"
    script:
        "../scripts/prepare_input.R"


# Step 2 (optional) -- Zhao4 subtype-specific quantile centering.
rule zhao_normalize:
    input:
        pam50_expr=out("pam50_expr.tsv"),
        samples_filtered=out("samples_filtered.tsv"),
    output:
        out("pam50_zhao_normalized.tsv"),
    log:
        out("logs", "zhao_normalize.log"),
    container:
        PAM50_IMG
    # conda:  # conda fallback for platforms without singularity
    #     "../envs/pam50.yaml"
    script:
        "../scripts/zhao_normalize.R"


# The classifier reads the Zhao-normalised matrix when enabled, else the raw
# PAM50 matrix (the classifier applies its own median-centering).
def _classify_input(wildcards):
    if config["zhao_normalization"]["enabled"]:
        return out("pam50_zhao_normalized.tsv")
    return out("pam50_expr.tsv")


# Step 3 -- PAM50 classification (Spearman correlation to Perou centroids).
rule classify:
    input:
        expr=_classify_input,
    output:
        scores=SCORES,
    log:
        out("logs", "classify.log"),
    container:
        PAM50_IMG
    # conda:  # conda fallback for platforms without singularity
    #     "../envs/pam50.yaml"
    script:
        "../scripts/classify.R"


# Step 4 -- enrich the raw scores with sample metadata -> pam50_calls.tsv.
rule finalize:
    input:
        scores=SCORES,
        samples_filtered=out("samples_filtered.tsv"),
    output:
        out("pam50_calls.tsv"),
    log:
        out("logs", "finalize.log"),
    container:
        PAM50_IMG
    # conda:  # conda fallback for platforms without singularity
    #     "../envs/pam50.yaml"
    script:
        "../scripts/finalize.R"


# Step 5 (optional) -- intrinsic clustering, claudin-low, silhouette, heatmap.
rule intrinsic_clustering:
    input:
        scores=SCORES,
        samples_filtered=out("samples_filtered.tsv"),
        int_expr=out("intrinsic_heatmap", "int_expr.tsv"),
        int_expr_medctr=out("intrinsic_heatmap", "int_expr_medCtr.tsv"),
    output:
        calls_clow=out("pam50_calls_clow.tsv"),
        heatmap=out("intrinsic_heatmap", "intrinsic_medctr_heatmap.pdf"),
        silhouette=out("intrinsic_heatmap", "silhouette_plot.pdf"),
    log:
        out("logs", "intrinsic_clustering.log"),
    container:
        PAM50_IMG
    # conda:  # conda fallback for platforms without singularity
    #     "../envs/pam50.yaml"
    script:
        "../scripts/intrinsic_clustering.R"

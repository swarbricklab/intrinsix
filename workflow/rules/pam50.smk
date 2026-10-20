# PAM50 rules. Each step is a thin wrapper that reads TSVs, calls the pam50perou
# package, and writes TSVs. The algorithm + reference data live in the package
# (workflow/pkg/pam50perou), baked into the `pam50` container.
#
# NOTE: the workflow/scripts/*.R wrappers are the next milestone; see DESIGN.md.

PAM50_IMG = config["containers"]["pam50"]


# Step 1 -- filter samples, map symbols -> pcrID, subset PAM50 + intrinsic genes,
# median-centre the intrinsic matrix.
rule prepare_input:
    input:
        expression=config["inputs"]["expression"],
        samples=config["inputs"]["samples"],
    output:
        pam50_expr=f"{OUT}/pam50_expr.tsv",
        samples_filtered=f"{OUT}/samples_filtered.tsv",
        int_expr=f"{OUT}/intrinsic_heatmap/int_expr.tsv",
        int_expr_medctr=f"{OUT}/intrinsic_heatmap/int_expr_medCtr.tsv",
    log:
        f"{OUT}/logs/prepare_input.log",
    container:
        PAM50_IMG
    # conda:  # conda fallback for platforms without singularity
    #     "../envs/pam50.yaml"
    script:
        "../scripts/prepare_input.R"


# Step 2 (optional) -- Zhao4 subtype-specific quantile centering.
rule zhao_normalize:
    input:
        pam50_expr=f"{OUT}/pam50_expr.tsv",
        samples_filtered=f"{OUT}/samples_filtered.tsv",
    output:
        f"{OUT}/pam50_zhao_normalized.tsv",
    log:
        f"{OUT}/logs/zhao_normalize.log",
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
        return f"{OUT}/pam50_zhao_normalized.tsv"
    return f"{OUT}/pam50_expr.tsv"


# Step 3 -- PAM50 classification (Spearman correlation to Perou centroids).
rule classify:
    input:
        expr=_classify_input,
    output:
        scores=f"{OUT}/{config['classifier']['prefix']}_pam50scores.txt",
    log:
        f"{OUT}/logs/classify.log",
    container:
        PAM50_IMG
    # conda:  # conda fallback for platforms without singularity
    #     "../envs/pam50.yaml"
    script:
        "../scripts/classify.R"


# Step 4 -- enrich the raw scores with sample metadata -> pam50_calls.tsv.
rule finalize:
    input:
        scores=f"{OUT}/{config['classifier']['prefix']}_pam50scores.txt",
        samples_filtered=f"{OUT}/samples_filtered.tsv",
    output:
        f"{OUT}/pam50_calls.tsv",
    log:
        f"{OUT}/logs/finalize.log",
    container:
        PAM50_IMG
    # conda:  # conda fallback for platforms without singularity
    #     "../envs/pam50.yaml"
    script:
        "../scripts/finalize.R"


# Step 5 (optional) -- intrinsic clustering, claudin-low, silhouette, heatmap.
rule intrinsic_clustering:
    input:
        scores=f"{OUT}/{config['classifier']['prefix']}_pam50scores.txt",
        samples_filtered=f"{OUT}/samples_filtered.tsv",
        int_expr=f"{OUT}/intrinsic_heatmap/int_expr.tsv",
        int_expr_medctr=f"{OUT}/intrinsic_heatmap/int_expr_medCtr.tsv",
    output:
        calls_clow=f"{OUT}/pam50_calls_clow.tsv",
        heatmap=f"{OUT}/intrinsic_heatmap/intrinsic_medctr_heatmap.pdf",
        silhouette=f"{OUT}/intrinsic_heatmap/silhouette_plot.pdf",
    log:
        f"{OUT}/logs/intrinsic_clustering.log",
    container:
        PAM50_IMG
    # conda:  # conda fallback for platforms without singularity
    #     "../envs/pam50.yaml"
    script:
        "../scripts/intrinsic_clustering.R"

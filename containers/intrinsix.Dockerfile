# Container image for the intrinsix workflow: R + the PAM50 dependency stack, with the
# intrinsix package (algorithm + bundled reference data + CLI) installed on top.
#
# Build (via absconda, from the repo root so the build context includes the pkg):
#   module load absconda/0.2.8
#   absconda build \
#     --dockerfile containers/intrinsix.Dockerfile \
#     --context . \
#     --repository ghcr.io/swarbricklab/intrinsix \
#     --tag "$(date +%Y%m%d)" \
#     --push
#
# The image tag is the consuming DVC stage's choice (not pinned in-repo).

FROM mambaorg/micromamba:1.5.10

USER root
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates && \
    rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

# --- conda environment (R + deps + testthat) ---
COPY --chown=$MAMBA_USER:$MAMBA_USER containers/intrinsix-env.yaml /tmp/intrinsix-env.yaml
RUN micromamba install -y -n base -f /tmp/intrinsix-env.yaml && \
    micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Put the base env on PATH unconditionally. `singularity exec` bypasses the
# micromamba entrypoint that would normally activate the env, so without this
# Rscript/R are not found at runtime (only during build via MAMBA_DOCKERFILE_ACTIVATE).
ENV PATH=/opt/conda/bin:$PATH

# --- install the intrinsix package (pure R: fast, no compilation) ---
# The package now lives at the repo root; copy just its components into the image.
COPY --chown=$MAMBA_USER:$MAMBA_USER DESCRIPTION NAMESPACE LICENSE /tmp/pkg/
COPY --chown=$MAMBA_USER:$MAMBA_USER R /tmp/pkg/R/
COPY --chown=$MAMBA_USER:$MAMBA_USER inst /tmp/pkg/inst/
COPY --chown=$MAMBA_USER:$MAMBA_USER man /tmp/pkg/man/
RUN R CMD INSTALL --with-keep.source /tmp/pkg && \
    rm -rf /tmp/pkg /tmp/intrinsix-env.yaml

# --- expose the `intrinsix` CLI on PATH ---
RUN chmod +x /opt/conda/lib/R/library/intrinsix/exec/intrinsix && \
    ln -sf /opt/conda/lib/R/library/intrinsix/exec/intrinsix /opt/conda/bin/intrinsix

# --- sanity: package loads, reference data present, CLI runs ---
RUN Rscript -e 'library(intrinsix); stopifnot(nrow(pam50_centroids()) == 50L, ncol(sigma_table()) >= 4L)' && \
    intrinsix --help

# Default to a login shell with the env active (singularity exec picks up PATH).
ENV R_LIBS_USER=/opt/conda/lib/R/library

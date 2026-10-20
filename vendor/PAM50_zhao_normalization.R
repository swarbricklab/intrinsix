###
# Applying Subtype-specific gene centering (Zhao4 Normalization)
# 3/3/25
# Alec Lobanov
# Perou Lab
#
# Paths updated for bcarc_wts repo structure (Mar 2026)
###

# Load necessary libraries
library(dplyr)
library(data.table)
library(yaml)
set.seed(1)

config <- yaml::read_yaml("config/PAM50/config.yaml")

# Read gene expression data (PAM50 genes × samples, log2_uqn)
gene_expr <- read.table(config$intermediates$pam50_expr, sep = '\t', header = TRUE, row.names = 1)
print(dim(gene_expr))

# Load metadata: gussid → clinical_subtype
# Written by prepare_pam50_input.py; uses the same 4-category encoding as
# the SIGMA table (ERpos_HER2neg, ERneg_HER2neg, HER2pos_ERpos, HER2pos_ERneg)
metadata <- read.csv(config$intermediates$pam50_metadata)
print(head(metadata))

# Read sigma data (quantile) for each PAM50 gene
# From perou-lab/PAM50_ER_HER2_normalization submodule
sigma <- read.table(config$deps$sigma, sep = '\t', header = TRUE, row.names = 1)
sigma[1:5,1:5]

# Function for quantile centering
quantile_centering <- function(expr_matrix, gene_quantile) {
  res <- expr_matrix
  for (name in rownames(expr_matrix)) {
    q <- gene_quantile[name, ]
    q_value <- quantile(as.numeric(expr_matrix[name, ]), probs = q, na.rm = TRUE)
    res[name, ] <- expr_matrix[name, ] - q_value
  }
  return(res)
}

# Run quantile centering per clinical subtype
centered_matrix <- list()
for (subtype in unique(metadata$clinical_subtype)) {
  samples <- metadata$gussid[metadata$clinical_subtype == subtype]
  gene_expr_subtype <- gene_expr[, colnames(gene_expr) %in% samples]
  percentile_group <- sigma[,subtype,drop=F]
  centered_matrix[[subtype]] <- quantile_centering(gene_expr_subtype, percentile_group)
}

# Combine results
normalized_matrix <- do.call(cbind, centered_matrix)
colnames(normalized_matrix) <- tstrsplit(colnames(normalized_matrix),"\\.")[[2]]


# Save output to file
write.table(data.frame(Gene=rownames(normalized_matrix),normalized_matrix), 
            config$intermediates$pam50_zhao_normalized,
            sep = '\t', row.names = FALSE, quote = FALSE)

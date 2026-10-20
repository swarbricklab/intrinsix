#' Identify Claudin-Low Samples
#'
#' @param cluster_results The list output from replicate_cluster3 (contains atr, gtr, cdt).
#' @param expr_matrix The original median-centered expression matrix (genes as rows, samples as columns).
#' @param array_corr_cutoff Minimum node correlation for sample clusters. Default 0.5.
#' @param gene_corr_cutoff Minimum node correlation for gene cluster. Default 0.6.
#' @param max_gene_node_size Maximum genes allowed in the claudin gene node. Default 20.
#' @param median_expr_cutoff Claudin expression must be below this. Default -1.
identify_claudin_low <- function(cluster_results, 
                                 expr_matrix,
                                 array_corr_cutoff = 0.5,
                                 gene_corr_cutoff = 0.6,
                                 max_gene_node_size = 20,
                                 median_expr_cutoff = -1) {
  
  # --- Setup and Mapping ---
  CLAUDIN_GENES <- c("CLDN3", "CLDN4", "CLDN7")
  atr <- cluster_results$atr
  gtr <- cluster_results$gtr
  
  # Map GENE/ARRY IDs back to names based on original matrix indices
  # GENEX indices are 0-based based on original row/col order
  gene_names   <- rownames(expr_matrix)
  sample_names <- colnames(expr_matrix)
  
  gene_id_map    <- setNames(gene_names, paste0("GENE", seq_along(gene_names) - 1, "X"))
  arry_to_sample <- setNames(sample_names, paste0("ARRY", seq_along(sample_names) - 1, "X"))
  
  # --- Helper: Recursively collect leaf IDs ---
  get_leaves <- function(node_id, tree_dict, leaf_prefix) {
    if (startsWith(node_id, leaf_prefix)) return(node_id)
    children <- tree_dict[[node_id]]
    if (is.null(children)) return(character(0))
    c(get_leaves(children[1], tree_dict, leaf_prefix),
      get_leaves(children[2], tree_dict, leaf_prefix))
  }
  
  # Build lookup dictionaries
  gtr_dict <- setNames(Map(c, gtr$LEFT, gtr$RIGHT), gtr$NODE)
  atr_dict <- setNames(Map(c, atr$LEFT, atr$RIGHT), atr$NODE)
  
  # --- STEP 1: Identify the claudin gene node ---
  message("Finding claudin gene node...")
  
  claudin_gids <- names(gene_id_map)[gene_id_map %in% CLAUDIN_GENES]
  if (length(claudin_gids) < length(CLAUDIN_GENES)) {
    stop("Not all claudin genes found in matrix: ", 
         paste(setdiff(CLAUDIN_GENES, gene_names), collapse = ", "))
  }
  
  # Evaluate gene nodes meeting the correlation threshold
  qualifying_gene_nodes <- lapply(
    which(gtr$SIMILARITY > gene_corr_cutoff),
    function(i) {
      nd     <- gtr$NODE[i]
      corr   <- gtr$SIMILARITY[i]
      leaves <- get_leaves(nd, gtr_dict, "GENE")
      n_genes <- length(leaves)
      
      if (n_genes > max_gene_node_size) return(NULL)
      if (!all(claudin_gids %in% leaves)) return(NULL)
      
      list(node = nd, n_genes = n_genes, correlation = corr,
           gene_names = gene_id_map[leaves])
    }
  )
  qualifying_gene_nodes <- Filter(Negate(is.null), qualifying_gene_nodes)
  
  if (length(qualifying_gene_nodes) == 0) {
    stop("No gene node found containing all claudins with correlation > ", 
         gene_corr_cutoff, " and size <= ", max_gene_node_size)
  }
  
  # Select largest qualifying node
  n_genes_vec    <- sapply(qualifying_gene_nodes, `[[`, "n_genes")
  best_gene_node <- qualifying_gene_nodes[[which.max(n_genes_vec)]]
  
  message(sprintf("  Selected gene node: %s | correlation=%.6f | n_genes=%d",
                  best_gene_node$node, best_gene_node$correlation, best_gene_node$n_genes))
  message("  Genes: ", paste(sort(best_gene_node$gene_names), collapse = ", "))
  
  # --- STEP 2: Identify claudin-low sample clusters ---
  message("Finding claudin-low sample clusters...")
  
  claudin_low_samples <- character(0)
  node = character(0)
  similarity = character(0)
  cldn_exp = character(0)
  
  # Iterate through array nodes meeting correlation threshold
  for (i in which(atr$SIMILARITY > array_corr_cutoff)) {
    nd <- atr$NODE[i]
    sim <- atr$SIMILARITY[i]
    #print(nd)
    
    arry_leaves <- get_leaves(nd, atr_dict, "ARRY")
    samp_names  <- arry_to_sample[arry_leaves]
    samp_names  <- samp_names[!is.na(samp_names)]
    
    if (length(samp_names) == 0) next
    
    # Calculate median expression for claudins in this cluster
    claudin_expr <- expr_matrix[CLAUDIN_GENES, samp_names, drop = FALSE]
    gene_medians <- apply(claudin_expr, 1, median, na.rm = TRUE)
    
    if (all(gene_medians < median_expr_cutoff)) {
      claudin_low_samples <- union(claudin_low_samples, samp_names)
      node = nd
      similarity = sim
      cldn_exp = gene_medians
    }
  }
  
  claudin_low_samples <- sort(claudin_low_samples)
  message(sprintf("  Selected array node: %s | correlation=%.6f | n_samples=%d",
                  node, similarity, length(claudin_low_samples)))
  message("  Samples: ", paste(sort(claudin_low_samples), collapse = ", "))
  message(sprintf("Median Claudins Expression : CLDN3 = %.4f | CLDN4 = %.4f | CLDN7 = %.4f",
          cldn_exp[1],cldn_exp[2],cldn_exp[3]))
  
  return(claudin_low_samples)
}

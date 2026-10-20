#' Replicate Cluster 3.0 and return data structures with progress messages
#'
#' @param input File path (string) or a data frame/matrix.
#' @param prefix Optional prefix for file writing.
#' @param write_files Logical. If TRUE, saves .atr, .gtr, and .cdt files to disk.
centroid_clustering <- function(input, prefix = deparse(substitute(input)), write_files = TRUE, path = ".") {
  
  # --- Internal Helper: Core Clustering ---
  centroid_linkage_cluster3 <- function(obs, label = "") {
    n   <- nrow(obs)
    pfx <- if (nzchar(label)) sprintf("[%s] ", label) else ""
    
    cat(sprintf("%sBuilding %dx%d correlation matrix...\n", pfx, n, n))
    t0  <- proc.time()["elapsed"]
    sim <- cor(t(obs))
    diag(sim) <- -Inf
    cat(sprintf("%sDone in %.1fs\n", pfx, proc.time()["elapsed"] - t0))
    
    max_nodes <- 2L * n - 1L
    full_sim  <- matrix(-Inf, max_nodes, max_nodes)
    full_sim[1:n, 1:n] <- sim
    full_cent <- matrix(0, max_nodes, ncol(obs))
    full_cent[1:n, ] <- obs
    
    sizes      <- numeric(max_nodes);  sizes[1:n]      <- 1
    stored_sim <- numeric(max_nodes);  stored_sim[1:n] <- Inf
    active     <- logical(max_nodes);  active[1:n]     <- TRUE
    
    merges  <- vector("list", n - 1L)
    next_id <- n + 1L
    
    cat(sprintf("%sClustering %d items...\n", pfx, n))
    t1 <- proc.time()["elapsed"]
    
    for (step in seq_len(n - 1L)) {
      if (n > 300L && step %% 200L == 0L)
        cat(sprintf("%s  Step %d/%d  (%.1fs)\n", pfx, step, n - 1L, proc.time()["elapsed"] - t1))
      
      active_idx <- which(active)
      m          <- length(active_idx)
      sub_ut     <- full_sim[active_idx, active_idx]
      sub_ut[!upper.tri(sub_ut)] <- -Inf
      
      flat    <- which.max(sub_ut)
      col_pos <- (flat - 1L) %/% m + 1L
      row_pos <- (flat - 1L) %% m + 1L
      
      ci <- active_idx[row_pos]; cj <- active_idx[col_pos]
      cv <- full_sim[ci, cj]
      
      new_id <- next_id; next_id <- next_id + 1L
      sz_i <- sizes[ci]; sz_j <- sizes[cj]; new_sz <- sz_i + sz_j
      new_cent <- (full_cent[ci, ] * sz_i + full_cent[cj, ] * sz_j) / new_sz
      new_stored <- min(cv, stored_sim[ci], stored_sim[cj])
      
      merges[[step]] <- list(left = ci, right = cj, sim = new_stored)
      
      active[ci] <- FALSE; active[cj] <- FALSE; active[new_id] <- TRUE
      full_cent[new_id, ] <- new_cent
      sizes[new_id]       <- new_sz
      stored_sim[new_id]  <- new_stored
      
      ar <- which(active); ar <- ar[ar != new_id]
      if (length(ar) > 0L) {
        others <- full_cent[ar, , drop = FALSE]
        a <- new_cent - mean(new_cent)
        b <- others - rowMeans(others)
        da <- max(sqrt(sum(a^2)), 1e-15)
        db <- pmax(sqrt(rowSums(b^2)), 1e-15)
        r <- as.vector((b %*% a) / (da * db))
        full_sim[new_id, ar] <- r; full_sim[ar, new_id] <- r
      }
      full_sim[new_id, new_id] <- -Inf
    }
    cat(sprintf("%sDone in %.1fs\n", pfx, proc.time()["elapsed"] - t1))
    return(merges)
  }
  
  # --- Internal Helper: Treesort ---
  cluster3_sort <- function(merges, n_items) {
    n_nodes <- length(merges)
    node_order <- numeric(n_nodes); node_count <- numeric(n_nodes)
    for (i in seq_len(n_nodes)) {
      m <- merges[[i]]
      l <- m$left; r <- m$right
      lo <- if (l <= n_items) as.numeric(l) else node_order[l - n_items]
      lc <- if (l <= n_items) 1 else node_count[l - n_items]
      ro <- if (r <= n_items) as.numeric(r) else node_order[r - n_items]
      rc <- if (r <= n_items) 1 else node_count[r - n_items]
      node_order[i] <- (lo * lc + ro * rc) / (lc + rc)
      node_count[i] <- lc + rc
    }
    order_out <- integer(n_items); pos <- 1L
    stack <- list(list(is_leaf = FALSE, idx = n_nodes))
    while (length(stack) > 0L) {
      top <- stack[[length(stack)]]; stack <- stack[-length(stack)]
      if (top$is_leaf) { order_out[pos] <- top$idx; pos <- pos + 1L }
      else {
        m <- merges[[top$idx]]; l <- m$left; r <- m$right
        lo <- if (l <= n_items) as.numeric(l) else node_order[l - n_items]
        ro <- if (r <= n_items) as.numeric(r) else node_order[r - n_items]
        if (lo <= ro) { first <- l; second <- r } else { first <- r; second <- l }
        for (child in c(second, first)) {
          if (child <= n_items) stack <- c(stack, list(list(is_leaf = TRUE, idx = child)))
          else stack <- c(stack, list(list(is_leaf = FALSE, idx = child - n_items)))
        }
      }
    }
    return(order_out)
  }
  
  format_tree_df <- function(merges, n_items, type = "GENE") {
    data.frame(
      NODE = paste0("NODE", seq_along(merges), "X"),
      LEFT = sapply(merges, function(m) if(m$left <= n_items) paste0(type, m$left - 1, "X") else paste0("NODE", m$left - n_items, "X")),
      RIGHT = sapply(merges, function(m) if(m$right <= n_items) paste0(type, m$right - 1, "X") else paste0("NODE", m$right - n_items, "X")),
      SIMILARITY = sapply(merges, function(m) m$sim),
      stringsAsFactors = FALSE
    )
  }
  
  # --- Data Ingest ---
  if (is.character(input)) {
    cat(sprintf("Loading %s...\n", input))
    data <- as.matrix(read.table(input, sep="\t", header=T, row.names=1, check.names=F))
  } else {
    data <- as.matrix(input)
  }
  
  genes <- rownames(data); cols <- colnames(data)
  n_genes <- nrow(data); n_arrays <- ncol(data)
  cat(sprintf("  %d genes x %d arrays\n", n_genes, n_arrays))
  
  # --- Execution ---
  cat("\n--- Array clustering ---\n")
  arr_merges <- centroid_linkage_cluster3(t(data), "arrays")
  arr_order  <- cluster3_sort(arr_merges, n_arrays)
  
  cat("\n--- Gene clustering ---\n")
  gene_merges <- centroid_linkage_cluster3(data, "genes")
  gene_order  <- cluster3_sort(gene_merges, n_genes)
  
  # Prepare results list
  res <- list(
    atr = format_tree_df(arr_merges, n_arrays, "ARRY"),
    gtr = format_tree_df(gene_merges, n_genes, "GENE"),
    cdt = data[gene_order, arr_order]
  )
  
  # --- File Writing (Optional) ---
  if (write_files) {
    cat("\n--- Writing output files ---\n")
    write.table(res$atr, paste0(path,prefix, ".atr"), sep="\t", quote=F, row.names=F, col.names=F)
    write.table(res$gtr, paste0(path,prefix, ".gtr"), sep="\t", quote=F, row.names=F, col.names=F)
    
    # CDT custom writer for Java TreeView compatibility
    con <- file(paste0(path,prefix, ".cdt"), "w")
    writeLines(paste(c("GID", "gene_name", "NAME", "GWEIGHT", colnames(res$cdt)), collapse = "\t"), con)
    writeLines(paste(c("AID", "", "", "", sprintf("ARRY%dX", arr_order - 1L)), collapse = "\t"), con)
    writeLines(paste(c("EWEIGHT", "", "", "", rep("1.000000", n_arrays)), collapse = "\t"), con)
    for (ri in seq_len(nrow(res$cdt))) {
      vals <- formatC(res$cdt[ri, ], format = "f", digits = 6)
      writeLines(paste(c(sprintf("GENE%dX", gene_order[ri] - 1L), rownames(res$cdt)[ri], rownames(res$cdt)[ri], "1.000000", vals), collapse = "\t"), con)
    }
    close(con)
    cat(sprintf("Files written: %s.atr, %s.gtr, %s.cdt\n", prefix, prefix, prefix))
  }
  
  cat("\nDone.\n")
  return(res)
}


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

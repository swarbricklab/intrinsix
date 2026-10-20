# Intrinsic clustering + claudin-low identification. Faithful port of
# vendor/clustering_scripts.R (Alec Lobanov, Perou Lab), which replicates Cluster
# 3.0 centroid-linkage hierarchical clustering (Pearson correlation) in base R and
# flags claudin-low samples from the resulting tree. The clustering math is kept
# verbatim; only roxygen and tidier output-path handling were added.

# --- internal: Cluster 3.0 centroid-linkage clustering of rows of `obs` ---
# @noRd
.centroid_linkage_cluster3 <- function(obs, label = "") {
  n <- nrow(obs)
  pfx <- if (nzchar(label)) sprintf("[%s] ", label) else ""

  sim <- stats::cor(t(obs))
  diag(sim) <- -Inf

  max_nodes <- 2L * n - 1L
  full_sim <- matrix(-Inf, max_nodes, max_nodes)
  full_sim[1:n, 1:n] <- sim
  full_cent <- matrix(0, max_nodes, ncol(obs))
  full_cent[1:n, ] <- obs

  sizes <- numeric(max_nodes); sizes[1:n] <- 1
  stored_sim <- numeric(max_nodes); stored_sim[1:n] <- Inf
  active <- logical(max_nodes); active[1:n] <- TRUE

  merges <- vector("list", n - 1L)
  next_id <- n + 1L

  for (step in seq_len(n - 1L)) {
    active_idx <- which(active)
    m <- length(active_idx)
    sub_ut <- full_sim[active_idx, active_idx]
    sub_ut[!upper.tri(sub_ut)] <- -Inf

    flat <- which.max(sub_ut)
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
    sizes[new_id] <- new_sz
    stored_sim[new_id] <- new_stored

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
  merges
}

# --- internal: Cluster 3.0 treesort leaf ordering ---
# @noRd
.cluster3_sort <- function(merges, n_items) {
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
    if (top$is_leaf) {
      order_out[pos] <- top$idx; pos <- pos + 1L
    } else {
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
  order_out
}

# @noRd
.format_tree_df <- function(merges, n_items, type = "GENE") {
  data.frame(
    NODE = paste0("NODE", seq_along(merges), "X"),
    LEFT = sapply(merges, function(m) if (m$left <= n_items) paste0(type, m$left - 1, "X") else paste0("NODE", m$left - n_items, "X")),
    RIGHT = sapply(merges, function(m) if (m$right <= n_items) paste0(type, m$right - 1, "X") else paste0("NODE", m$right - n_items, "X")),
    SIMILARITY = sapply(merges, function(m) m$sim),
    stringsAsFactors = FALSE
  )
}

#' Cluster 3.0 centroid-linkage clustering of an expression matrix
#'
#' Replicates Cluster 3.0 (Pearson correlation, centroid linkage) on both genes
#' (rows) and arrays (columns), producing TreeView-compatible `.atr`/`.gtr`/`.cdt`
#' structures. Used to drive claudin-low identification and the intrinsic heatmap.
#'
#' @param input Numeric matrix (genes x samples), or a path to a TSV of one.
#'   For claudin-low work this should be the median-centred intrinsic matrix.
#' @param prefix Output file prefix (used when `write_files = TRUE`).
#' @param write_files Logical; write `.atr`/`.gtr`/`.cdt` to `path`.
#' @param path Output directory for the cluster files.
#' @return A list with `atr`, `gtr` (data frames of the array/gene merge trees)
#'   and `cdt` (the reordered matrix).
#' @export
centroid_clustering <- function(input, prefix = deparse(substitute(input)),
                                write_files = TRUE, path = ".") {
  if (is.character(input)) {
    data <- as.matrix(read.table(input, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE))
  } else {
    data <- as.matrix(input)
  }
  genes <- rownames(data); cols <- colnames(data)
  n_genes <- nrow(data); n_arrays <- ncol(data)

  arr_merges <- .centroid_linkage_cluster3(t(data), "arrays")
  arr_order <- .cluster3_sort(arr_merges, n_arrays)
  gene_merges <- .centroid_linkage_cluster3(data, "genes")
  gene_order <- .cluster3_sort(gene_merges, n_genes)

  res <- list(
    atr = .format_tree_df(arr_merges, n_arrays, "ARRY"),
    gtr = .format_tree_df(gene_merges, n_genes, "GENE"),
    cdt = data[gene_order, arr_order]
  )

  if (write_files) {
    out <- function(ext) file.path(path, paste0(prefix, ext))
    write.table(res$atr, out(".atr"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(res$gtr, out(".gtr"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    con <- file(out(".cdt"), "w")
    on.exit(close(con))
    writeLines(paste(c("GID", "gene_name", "NAME", "GWEIGHT", colnames(res$cdt)), collapse = "\t"), con)
    writeLines(paste(c("AID", "", "", "", sprintf("ARRY%dX", arr_order - 1L)), collapse = "\t"), con)
    writeLines(paste(c("EWEIGHT", "", "", "", rep("1.000000", n_arrays)), collapse = "\t"), con)
    for (ri in seq_len(nrow(res$cdt))) {
      vals <- formatC(res$cdt[ri, ], format = "f", digits = 6)
      writeLines(paste(c(sprintf("GENE%dX", gene_order[ri] - 1L), rownames(res$cdt)[ri], rownames(res$cdt)[ri], "1.000000", vals), collapse = "\t"), con)
    }
  }
  res
}

#' Identify claudin-low samples from a clustering tree
#'
#' Finds the gene node containing CLDN3/CLDN4/CLDN7 (subject to correlation and
#' size thresholds), then flags sample clusters whose median claudin expression
#' falls below `median_expr_cutoff`. Faithful port of the pipeline logic.
#'
#' @param cluster_results Output of [centroid_clustering()] (uses `atr`/`gtr`).
#' @param expr_matrix The median-centred intrinsic expression matrix (genes x
#'   samples) that was clustered.
#' @param array_corr_cutoff Minimum array-node correlation. Default 0.5.
#' @param gene_corr_cutoff Minimum gene-node correlation. Default 0.6.
#' @param max_gene_node_size Maximum genes in the claudin gene node. Default 20.
#' @param median_expr_cutoff Claudin median expression must be below this.
#'   Default -1.
#' @return Sorted character vector of claudin-low sample IDs.
#' @export
identify_claudin_low <- function(cluster_results, expr_matrix,
                                 array_corr_cutoff = 0.5, gene_corr_cutoff = 0.6,
                                 max_gene_node_size = 20, median_expr_cutoff = -1) {
  CLAUDIN_GENES <- c("CLDN3", "CLDN4", "CLDN7")
  atr <- cluster_results$atr
  gtr <- cluster_results$gtr

  gene_names <- rownames(expr_matrix)
  sample_names <- colnames(expr_matrix)
  gene_id_map <- stats::setNames(gene_names, paste0("GENE", seq_along(gene_names) - 1, "X"))
  arry_to_sample <- stats::setNames(sample_names, paste0("ARRY", seq_along(sample_names) - 1, "X"))

  get_leaves <- function(node_id, tree_dict, leaf_prefix) {
    if (startsWith(node_id, leaf_prefix)) return(node_id)
    children <- tree_dict[[node_id]]
    if (is.null(children)) return(character(0))
    c(get_leaves(children[1], tree_dict, leaf_prefix),
      get_leaves(children[2], tree_dict, leaf_prefix))
  }

  gtr_dict <- stats::setNames(Map(c, gtr$LEFT, gtr$RIGHT), gtr$NODE)
  atr_dict <- stats::setNames(Map(c, atr$LEFT, atr$RIGHT), atr$NODE)

  claudin_gids <- names(gene_id_map)[gene_id_map %in% CLAUDIN_GENES]
  if (length(claudin_gids) < length(CLAUDIN_GENES)) {
    stop("Not all claudin genes found in matrix: ",
         paste(setdiff(CLAUDIN_GENES, gene_names), collapse = ", "))
  }

  qualifying_gene_nodes <- lapply(
    which(gtr$SIMILARITY > gene_corr_cutoff),
    function(i) {
      nd <- gtr$NODE[i]; corr <- gtr$SIMILARITY[i]
      leaves <- get_leaves(nd, gtr_dict, "GENE")
      n_g <- length(leaves)
      if (n_g > max_gene_node_size) return(NULL)
      if (!all(claudin_gids %in% leaves)) return(NULL)
      list(node = nd, n_genes = n_g, correlation = corr, gene_names = gene_id_map[leaves])
    }
  )
  qualifying_gene_nodes <- Filter(Negate(is.null), qualifying_gene_nodes)
  if (length(qualifying_gene_nodes) == 0) {
    stop("No gene node found containing all claudins with correlation > ",
         gene_corr_cutoff, " and size <= ", max_gene_node_size)
  }
  n_genes_vec <- sapply(qualifying_gene_nodes, `[[`, "n_genes")
  best_gene_node <- qualifying_gene_nodes[[which.max(n_genes_vec)]]

  claudin_low_samples <- character(0)
  for (i in which(atr$SIMILARITY > array_corr_cutoff)) {
    nd <- atr$NODE[i]
    arry_leaves <- get_leaves(nd, atr_dict, "ARRY")
    samp_names <- arry_to_sample[arry_leaves]
    samp_names <- samp_names[!is.na(samp_names)]
    if (length(samp_names) == 0) next
    claudin_expr <- expr_matrix[CLAUDIN_GENES, samp_names, drop = FALSE]
    gene_medians <- apply(claudin_expr, 1, stats::median, na.rm = TRUE)
    if (all(gene_medians < median_expr_cutoff)) {
      claudin_low_samples <- union(claudin_low_samples, samp_names)
    }
  }
  sort(claudin_low_samples)
}

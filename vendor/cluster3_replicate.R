# #!/usr/bin/env Rscript
# =============================================================================
# cluster3_replicate.R
#
# Replicates Cluster 3.0 hierarchical clustering output:
#   - Similarity metric : Pearson correlation (centered)
#   - Clustering method : Centroid linkage
#   - Clusters both genes (rows) and arrays (columns)
#   - Output            : .atr, .gtr, .cdt files compatible with Java TreeView
#
# Usage:
#   Rscript cluster3_replicate.R input.tsv [output_prefix]
#
#   If output_prefix is omitted, the input filename stem is used.
#   Example:
#     Rscript cluster3_replicate.R intrinsic_expr_medctr.tsv
#     -> produces intrinsic_expr_medctr.atr, .gtr, .cdt
#
# Requirements: base R only (no external packages)
#
# Algorithm notes:
#   Centroid linkage: new centroid = weighted mean of member vectors.
#   Merge selection : max Pearson correlation between current centroids.
#   Stored similarity: min(centroid_sim, left_stored, right_stored) — monotone
#     decreasing along the tree, matching Cluster 3.0's exact ATR/GTR values.
#   Left/right assignment: upper-triangle column-major scan (R default), which
#     matches the C Clustering Library loop  for(is=0; is<n-1; is++)
#     for(js=is+1; js<n; js++) — left = smaller active index, right = larger.
#   CDT ordering: Cluster 3.0 treesort — weighted-average leaf index orients each
#     subtree; left child first on ties.
# =============================================================================

# ---------------------------------------------------------------------------
# Core clustering
# ---------------------------------------------------------------------------

centroid_linkage_cluster3 <- function(obs, label = "") {
  # obs : matrix (n_samples x n_features)
  # Returns a list of (n_samples - 1) merge records, each:
  #   $left, $right  — 1-based indices into a (2*n-1)-length node array
  #                    (1..n = original leaves, n+1.. = internal nodes)
  #   $sim           — Cluster 3.0 stored similarity value

  n   <- nrow(obs)
  pfx <- if (nzchar(label)) sprintf("[%s] ", label) else ""

  cat(sprintf("%sBuilding %dx%d correlation matrix...\n", pfx, n, n))
  t0  <- proc.time()["elapsed"]

  # cor(t(obs)) gives n x n sample-sample Pearson correlation matrix
  sim <- cor(t(obs))
  diag(sim) <- -Inf
  cat(sprintf("%sDone in %.1fs\n", pfx, proc.time()["elapsed"] - t0))

  # Expand to full (2n-1) x (2n-1) matrix to accommodate internal nodes
  max_nodes <- 2L * n - 1L
  full_sim  <- matrix(-Inf, max_nodes, max_nodes)
  full_sim[1:n, 1:n] <- sim
  rm(sim)

  full_cent <- matrix(0, max_nodes, ncol(obs))
  full_cent[1:n, ] <- obs

  sizes      <- numeric(max_nodes);  sizes[1:n]      <- 1
  stored_sim <- numeric(max_nodes);  stored_sim[1:n] <- Inf   # leaves unconstrained
  active     <- logical(max_nodes);  active[1:n]     <- TRUE

  merges  <- vector("list", n - 1L)
  next_id <- n + 1L

  cat(sprintf("%sClustering %d items...\n", pfx, n))
  t1 <- proc.time()["elapsed"]

  for (step in seq_len(n - 1L)) {
    if (n > 300L && step %% 200L == 0L)
      cat(sprintf("%s  Step %d/%d  (%.1fs)\n", pfx, step, n - 1L,
                  proc.time()["elapsed"] - t1))

    active_idx <- which(active)
    m          <- length(active_idx)
    sub        <- full_sim[active_idx, active_idx]

    # Upper-triangle col-major scan matches C Clustering Library loop order:
    # for(is=0; is<n-1; is++) for(js=is+1; js<n; js++)
    # => left = smaller active position, right = larger active position
    sub_ut        <- sub
    sub_ut[!upper.tri(sub)] <- -Inf
    flat    <- which.max(sub_ut)
    col_pos <- (flat - 1L) %/% m + 1L   # col index within sub (1-based)
    row_pos <- (flat - 1L) %% m + 1L    # row index within sub (1-based)

    ci <- active_idx[row_pos]   # left  child (smaller active index)
    cj <- active_idx[col_pos]   # right child (larger  active index)
    cv <- full_sim[ci, cj]

    new_id  <- next_id;  next_id <- next_id + 1L
    sz_i    <- sizes[ci];  sz_j <- sizes[cj];  new_sz <- sz_i + sz_j
    new_cent   <- (full_cent[ci, ] * sz_i + full_cent[cj, ] * sz_j) / new_sz

    # Cluster 3.0 stored similarity: monotone minimum along tree path
    new_stored <- min(cv, stored_sim[ci], stored_sim[cj])

    merges[[step]] <- list(left = ci, right = cj, sim = new_stored)

    active[ci] <- FALSE;  active[cj] <- FALSE;  active[new_id] <- TRUE
    full_cent[new_id, ] <- new_cent
    sizes[new_id]       <- new_sz
    stored_sim[new_id]  <- new_stored

    # Vectorised centroid similarity update
    ar <- which(active)
    ar <- ar[ar != new_id]
    if (length(ar) > 0L) {
      others <- full_cent[ar, , drop = FALSE]
      a      <- new_cent - mean(new_cent)
      b      <- others  - rowMeans(others)
      da     <- sqrt(sum(a^2));   if (da < 1e-15) da <- 1e-15
      db     <- sqrt(rowSums(b^2));  db[db < 1e-15] <- 1e-15
      r      <- as.vector((b %*% a) / (da * db))
      full_sim[new_id, ar] <- r
      full_sim[ar, new_id] <- r
    }
    full_sim[new_id, new_id] <- -Inf
  }

  cat(sprintf("%sDone in %.1fs\n", pfx, proc.time()["elapsed"] - t1))
  merges
}

# ---------------------------------------------------------------------------
# Leaf ordering (Cluster 3.0 treesort)
# ---------------------------------------------------------------------------

cluster3_sort <- function(merges, n_items) {
  # Compute Cluster 3.0 leaf ordering from merge list.
  # Uses weighted-average leaf index to orient each subtree;
  # left child (as recorded) comes first on ties.
  #
  # merges : list of length n_items-1, each with $left, $right (1-based, leaf=1..n)
  # Returns integer vector of original item indices (1-based) in CDT display order.

  n_nodes    <- length(merges)
  node_order <- numeric(n_nodes)
  node_count <- numeric(n_nodes)

  # Build node_order (weighted-average position) bottom-up
  for (i in seq_len(n_nodes)) {
    m  <- merges[[i]]
    l  <- m$left;   r <- m$right;   n_items_local <- n_items
    if (l <= n_items) { lo <- as.numeric(l); lc <- 1 } else {
      ni <- l - n_items; lo <- node_order[ni]; lc <- node_count[ni] }
    if (r <= n_items) { ro <- as.numeric(r); rc <- 1 } else {
      ni <- r - n_items; ro <- node_order[ni]; rc <- node_count[ni] }
    node_order[i] <- (lo * lc + ro * rc) / (lc + rc)
    node_count[i] <- lc + rc
  }

  # Iterative DFS — stack holds (is_leaf, index)
  order_out <- integer(n_items)
  pos       <- 1L
  stack     <- list(list(is_leaf = FALSE, idx = n_nodes))

  while (length(stack) > 0L) {
    top   <- stack[[length(stack)]]
    stack <- stack[-length(stack)]

    if (top$is_leaf) {
      order_out[pos] <- top$idx
      pos <- pos + 1L
    } else {
      m   <- merges[[top$idx]]
      l   <- m$left;  r <- m$right
      lo  <- if (l <= n_items) as.numeric(l) else node_order[l - n_items]
      ro  <- if (r <= n_items) as.numeric(r) else node_order[r - n_items]
      # Left child first on tie
      if (lo <= ro) { first <- l;  second <- r }
      else          { first <- r;  second <- l }
      # Push second first so first is processed first (LIFO)
      for (child in c(second, first)) {
        if (child <= n_items)
          stack <- c(stack, list(list(is_leaf = TRUE,  idx = child)))
        else
          stack <- c(stack, list(list(is_leaf = FALSE, idx = child - n_items)))
      }
    }
  }
  order_out
}

# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

write_atr <- function(merges, n_arrays, path) {
  lines <- character(length(merges))
  for (i in seq_along(merges)) {
    m  <- merges[[i]]
    l  <- m$left;  r <- m$right
    ls <- if (l <= n_arrays) sprintf("ARRY%dX", l - 1L) else sprintf("NODE%dX", l - n_arrays)
    rs <- if (r <= n_arrays) sprintf("ARRY%dX", r - 1L) else sprintf("NODE%dX", r - n_arrays)
    lines[i] <- sprintf("NODE%dX\t%s\t%s\t%f", i, ls, rs, m$sim)
  }
  writeLines(lines, path)
}

write_gtr <- function(merges, n_genes, path) {
  lines <- character(length(merges))
  for (i in seq_along(merges)) {
    m  <- merges[[i]]
    l  <- m$left;  r <- m$right
    ls <- if (l <= n_genes) sprintf("GENE%dX", l - 1L) else sprintf("NODE%dX", l - n_genes)
    rs <- if (r <= n_genes) sprintf("GENE%dX", r - 1L) else sprintf("NODE%dX", r - n_genes)
    lines[i] <- sprintf("NODE%dX\t%s\t%s\t%f", i, ls, rs, m$sim)
  }
  writeLines(lines, path)
}

write_cdt <- function(data, genes, cols, arr_order, gene_order, path) {
  # arr_order, gene_order: 1-based R indices
  data_reord  <- data[gene_order, arr_order]
  cols_reord  <- cols[arr_order]
  genes_reord <- genes[gene_order]
  arry_ids    <- sprintf("ARRY%dX", arr_order - 1L)
  gene_ids    <- sprintf("GENE%dX", gene_order - 1L)

  con <- file(path, "w")
  on.exit(close(con))

  # Row 1: column names
  writeLines(paste(c("GID", "gene_name", "NAME", "GWEIGHT", cols_reord), collapse = "\t"), con)
  # Row 2: ARRY IDs
  writeLines(paste(c("AID", "", "", "", arry_ids), collapse = "\t"), con)
  # Row 3: EWEIGHT
  writeLines(paste(c("EWEIGHT", "", "", "", rep("1.000000", length(arr_order))), collapse = "\t"), con)
  # Data rows
  for (ri in seq_len(nrow(data_reord))) {
    vals <- formatC(data_reord[ri, ], format = "f", digits = 6)
    writeLines(paste(c(gene_ids[ri], genes_reord[ri], genes_reord[ri], "1.000000", vals),
                     collapse = "\t"), con)
  }
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  cat("Usage: Rscript cluster3_replicate.R input.tsv [output_prefix]\n")
  quit(status = 1L)
}

input_path <- args[1]
prefix     <- if (length(args) >= 2L) args[2] else sub("\\.[^.]+$", "", input_path)

cat(sprintf("Loading %s...\n", input_path))
tsv   <- read.table(input_path, sep = "\t", header = TRUE, row.names = 1,
                    check.names = FALSE)
data  <- as.matrix(tsv)
genes <- rownames(data)
cols  <- colnames(data)
n_genes  <- nrow(data)
n_arrays <- ncol(data)
cat(sprintf("  %d genes x %d arrays\n", n_genes, n_arrays))

# --- Cluster arrays ---
cat("\n--- Array clustering ---\n")
arr_merges <- centroid_linkage_cluster3(t(data), "arrays")  # obs=(arrays x genes)
arr_order  <- cluster3_sort(arr_merges, n_arrays)
write_atr(arr_merges, n_arrays, paste0(prefix, ".atr"))
cat(sprintf("Written: %s.atr\n", prefix))

# --- Cluster genes ---
cat("\n--- Gene clustering ---\n")
gene_merges <- centroid_linkage_cluster3(data, "genes")  # obs=(genes x arrays)
gene_order  <- cluster3_sort(gene_merges, n_genes)
write_gtr(gene_merges, n_genes, paste0(prefix, ".gtr"))
cat(sprintf("Written: %s.gtr\n", prefix))

# --- Write CDT ---
cat("\n--- Writing CDT ---\n")
write_cdt(data, genes, cols, arr_order, gene_order, paste0(prefix, ".cdt"))
cat(sprintf("Written: %s.cdt\n", prefix))

cat("\nDone.\n")


#' Differential expression analysis between two group of cells using `limma`
#'
#' Performs differential expression analysis between two groups of cells using the
#' \code{limma} package. Supports optional observation-level weights and filters genes
#' based on minimum expression percent in either group.
#'
#' @param exp_mat A normalized gene expression matrix with genes as rows and cells as columns. For example, the log-normalized data in a Seurat object.
#' @param cells_reference A vector of cell IDs for the reference (baseline) group.
#' @param cells_target A vector of cell IDs for the target (comparison) group.
#' @param weights Optional numeric vector of observation-level weights for each cell (same length as combined groups). For example, the scaled spatial distance between cells and the boundary.
#' @param adj_p.value Adjusted p-value threshold for significance to obtain differentially expressed genes. Default is 0.05.
#' @param min.pct Minimum proportion of cells expressing the gene in at least one group. Default is 0.
#'
#' @return A data frame of differentially expressed genes between two cell groups
#'
#' @export
RunLimmaDE <- function(exp_mat = NULL,
                       cells_reference = NULL,
                       cells_target = NULL,
                       weights = NULL,
                       adj_p.value = 0.05,
                       min.pct = 0) {

  # --- Input checks ---
  if (is.null(exp_mat) || !inherits(exp_mat, c("matrix", "dgCMatrix"))) {
    stop("'exp_mat' must be a non-null numeric matrix (genes by cells).")
  }

  all_cells <- c(cells_reference, cells_target)

  if (is.null(cells_reference) || is.null(cells_target)) {
    stop("Both 'cells_reference' and 'cells_target' must be provided.")
  }
  if (!all(all_cells %in% colnames(exp_mat))) {
    stop("Some cell IDs not found in column names of 'exp_mat'.")
  }

  if (!is.null(weights)){
    if( is.null(names(weights)) || (length(weights) != length(all_cells)) ) {
      stop("'weights' should be a named vector, and length of 'weights' must match the total number of selected cells.")
    }

    # Make the weights of correct order
    weights <- weights[all_cells]
  }

  # --- Convert exp_mat into sparse matrix if it is a matrix ---
  if( inherits(exp_mat,"matrix") ){
    exp_mat <- Matrix::Matrix(exp_mat, sparse = TRUE)
    }

  # --- Subset expression matrix ---
  expr <- exp_mat[, all_cells, drop = FALSE]

  # --- Compute percentage of expressing cells per group ---
  df_pct <- data.frame(
    pct.reference = Matrix::rowMeans(expr[, cells_reference] > 0),
    pct.target = Matrix::rowMeans(expr[, cells_target] > 0)
  )
  rownames(df_pct) <- rownames(expr)

  # --- Filter genes by min.pct ---
  keep_genes <- rowSums(df_pct >= min.pct) >= 1
  expr <- expr[keep_genes, ]
  df_pct <- df_pct[keep_genes, ]

  if (nrow(expr) == 0) {
    stop("No genes passed the 'min.pct' filter.")
  }

  # --- Design matrix and group labels ---
  group <- factor(c(rep(1, length(cells_reference)), rep(2, length(cells_target))))
  design <- stats::model.matrix(~ group)

  # --- Run limma ---
  fit <- limma::lmFit(expr, design, weights = weights)
  fit <- limma::eBayes(fit)
  tab <- limma::topTable(fit, coef = 2, number = Inf, p.value = adj_p.value)

  if (nrow(tab) == 0) {
    warning("No genes passed the adjusted p-value threshold.")
    return(tab)
  }

  # --- Add gene stats ---
  tab <- cbind(tab, df_pct[rownames(tab), ])
  tab$gene <- rownames(tab)
  tab <- tab[order(-abs(tab$logFC)), ]

  return(tab)
}



#' Differential expression analysis along spatial distance gradients
#'
#'
#' @inheritParams RunLimmaDE
#' @param cells_ids A vector of cell IDs in `exp_mat` for spatial differential expression analysis.
#' @param spatial_distance A numeric vector of the spatial distance for each cell. Scaled distance rather than real distance is is recommended.
#'
#' @return A data.frame of differentially expressed genes along spatial distance. Positive trend means positive correlation of the gene expression with the spatial distance (Up regulated), and negative trend means negative correlation (Down regulated).
#'
#' @export
RunSpatialDE <- function(exp_mat = NULL,
                         cells_ids = NULL,
                         spatial_distance = NULL,
                         adj_p.value = 0.05
                         ) {

  # --- Input checks ---
  if (is.null(exp_mat) || !inherits(exp_mat, c("matrix", "dgCMatrix"))) {
    stop("'exp_mat' must be a non-null numeric matrix (genes by cells).")
  }

  if (is.null(cells_ids)) {
    stop("'cells_ids' must be provided.")
  }

  if (!all(cells_ids %in% colnames(exp_mat))) {
    stop("Some cell IDs not found in column names of 'exp_mat'.")
  }

  all_cells <- cells_ids

  if (!is.null(spatial_distance)){
    if( is.null(names(spatial_distance)) || (length(spatial_distance) != length(all_cells)) ) {
      stop("'spatial_distance' should be a named vector, and length of 'spatial_distance' must match the total number of selected cells.")
    }

  }

  # --- Make spatial_distance numeric ---
  t1 <- as.numeric(spatial_distance[all_cells])

  # --- Subset expression matrix ---
  expr <- exp_mat[, all_cells, drop = FALSE]

  # ---  Design matrix ---
  X <- splines::ns(as.numeric(t1),df = 3)
  A <- cbind(1,t1,X)
  QR <- qr(A)
  r <- QR$rank
  R_rank <- QR$qr[1:r,1:r]
  Z <- t(backsolve(R_rank,t(A),transpose=TRUE))
  Z <- Z[,-1]
  design <- stats::model.matrix(~ Z)

  # --- Run limma ---
  fit <- limma::lmFit(expr, design)
  fit <- limma::eBayes(fit)
  tab <- limma::topTable(fit, coef = 2:4, number = Inf, p.value = adj_p.value)

  if (nrow(tab) == 0) {
    warning("No genes passed the adjusted p-value threshold.")
    return(tab)
  }

  # --- Add DE gene trend ---
  tab$gene <- rownames(tab)
  tab$trend <- ifelse(tab$Z1 > 0, "Positive", "Negative")

  return(tab)
}



#' Differential expression analysis between two group of cells using `limma`
#'
#' Performs differential expression analysis between two groups of cells using the
#' `limma` package. Supports optional observation-level weights and filters genes
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
#'
#' @examples
#'
#' # Load coordinates of cells in mouse brain
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'                               package = "SpNeigh"))
#' head(coords)
#' table(coords$cluster)
#'
#' # Load log normalized expression data of mouse brain,
#' # only cells in cluster 0 and cluster 6 are included
#' logNorm_expr <- readRDS(system.file("extdata", "LogNormExpr.rds",
#'                               package = "SpNeigh"))
#' class(logNorm_expr)
#'
#' # Obtain reference cells and target cells for DE anallysis
#' cells_ref <- subset(coords, cluster==0)[,"cell"]
#' cells_tar <- subset(coords, cluster==6)[,"cell"]
#'
#' # Perform DE analysis of target (C6) vs reference (C0) cells,
#' # Keep genes expressed >= 25% in either target or reference cells
#' tab <- RunLimmaDE(exp_mat = logNorm_expr,
#'                   cells_reference = cells_ref,
#'                   cells_target = cells_tar,
#'                   min.pct = 0.25)
#' head(tab)
#' head(tab[,c(1,5,7:9)])
#'
#' # Set cell names as rownames for coords
#' rownames(coords) <- coords$cell
#'
#' # Combine reference and target cells
#' all_cells <- c(cells_ref, cells_tar)
#'
#' # Plot expression of top DE genes,
#' # Make sure cell names in data and exp_mat match
#' PlotExpression(data = coords[all_cells,],
#'                exp_mat = logNorm_expr[,all_cells],
#'                genes = tab$gene[1],
#'                split_by = "cluster",
#'                angle_x_label = 45)
#'
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
#' Performs differential expression (DE) testing along a continuous spatial distance gradient.
#' Gene expression is modeled as a smooth function of spatial proximity using a natural spline basis.
#' This approach is useful for detecting genes with spatially varying expression relative to a boundary or centroid.
#'
#' @inheritParams RunLimmaDE
#' @inheritParams SplineDesign
#' @param cell_ids A vector of cell IDs in `exp_mat` for spatial differential expression analysis.
#' @param spatial_distance A numeric vector of the spatial distance for each cell. Scaled distance rather than real distance is is recommended.
#'
#' @return A data frame of differentially expressed genes, including p-values, adjusted p-values,
#'         and a directional trend column. The sign of Z1 gives the direction of the trend.
#'         A positive trend indicates increasing expression with distance;
#'         a negative trend indicates decreasing expression with distance.
#'         The values of Z1, Z2 and Z3 are not meaningful.
#'
#' @export
#'
#' @examples
#'
#' # Load coordinates of cells in mouse brain
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'                               package = "SpNeigh"))
#' head(coords)
#'
#' # Load log normalized expression data of mouse brain,
#' # only cells in cluster 0 and cluster 6 are included
#' logNorm_expr <- readRDS(system.file("extdata", "LogNormExpr.rds",
#'                               package = "SpNeigh"))
#' class(logNorm_expr)
#'
#' # Get and plot boundary of cluster 0
#' bon_c0 <- GetBoundary(data = coords, one_cluster = 0)
#' PlotBoundary(data = coords) + AddBoundary(bon_c0)
#'
#' # Compute boundary weights and plot weights
#' cells_c0 <- subset(coords, cluster==0)[,"cell"]
#' weights_bon <- ComputeBoundaryWeights(data = coords,
#'                                       cell_ids = cells_c0,
#'                                       boundary = bon_c0)
#' PlotWeights(data = coords, weights = weights_bon)
#'
#' # Perform DE analysis for cells in cluster 0
#' # along the distance (weights) to the boundaries
#' tab <- RunSpatialDE(exp_mat = logNorm_expr,
#'                     cell_ids = cells_c0,
#'                     spatial_distance = weights_bon)
#' head(tab)
#'
#' # Plot expression of top DE genes,
#' # Make sure cell names in data and exp_mat match
#' PlotExpression(data = coords[cells_c0,],
#'                exp_mat = logNorm_expr[,cells_c0],
#'                genes = tab$gene[1],
#'                angle_x_label = 45)
#'
RunSpatialDE <- function(exp_mat = NULL,
                         cell_ids = NULL,
                         spatial_distance = NULL,
                         adj_p.value = 0.05,
                         df = 3
                         ) {

  # --- Input checks ---
  if (is.null(exp_mat) || !inherits(exp_mat, c("matrix", "dgCMatrix"))) {
    stop("'exp_mat' must be a non-null numeric matrix (genes by cells).")
  }

  if (is.null(cell_ids)) {
    stop("'cell_ids' must be provided.")
  }

  if (!all(cell_ids %in% colnames(exp_mat))) {
    stop("Some cell IDs not found in column names of 'exp_mat'.")
  }

  all_cells <- cell_ids

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
  Z <- SplineDesign(t1, df = df)

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


#' Generate a spline-based orthonormal design matrix from a numeric covariate
#'
#' This function transforms a numeric input vector (e.g., spatial distance)
#' into a natural cubic spline basis and returns an orthonormalized design matrix
#' suitable for smooth trend modeling.
#' The first column of the output matrix is guaranteed to be positively correlated
#' with the input vector and typically captures the primary linear trend.
#'
#' @param x A numeric vector representing a continuous variable.
#' @param df Degrees of freedom for the spline basis. Default is 3.
#'
#' @return A numeric matrix with orthonormal columns derived from the spline-transformed input.
#'         The first column is aligned to show the main linear trend with respect to `x`.
#'
#' @export
#' @examples
#' x <- runif(100)
#' Z <- SplineDesign(x)
#'
#' # The first column of Z is positvely correlated with x
#' cor(Z[,1], x)
#'
SplineDesign <- function(x, df = 3) {
  if (!is.numeric(x)) {
    stop("'x' must be a numeric vector.")
  }

  # Spline basis
  X <- splines::ns(as.numeric(x), df = df)

  # Augment with intercept and original variable
  A <- cbind(1, x, X)

  # Orthonormalize via QR decomposition
  QR <- qr(A)
  r <- QR$rank
  R_rank <- QR$qr[1:r, 1:r]
  Z <- t(backsolve(R_rank, t(A), transpose = TRUE))

  # Remove intercept column
  Z <- Z[, -1]

  # Align first component with positive trend
  if (stats::cor(Z[, 1], x) < 0) {
    Z[, 1] <- -Z[, 1]
  }

  return(Z)
}


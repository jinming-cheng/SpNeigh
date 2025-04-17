#' Compute spatial enrichment index (SEI) for each gene
#'
#' Computes a spatial enrichment index (SEI) for each gene using precomputed spatial weights.
#' This generalizes to any distance-based weighting scheme (e.g. boundary, centroid, or other).
#'
#' @param exp_mat A gene expression matrix (genes × cells), either of class `matrix` or `dgCMatrix`.
#' @param weights A numeric vector of spatial weights, same length as the number of cells.
#'
#' @return A data frame with columns:
#' \itemize{
#'   \item `gene`: gene name
#'   \item `SEI`: spatial enrichment index (weighted mean expression)
#'   \item `mean_expr`: mean expression across all cells
#'   \item `normalized_SEI`: ratio of SEI to mean expression
#' }
#' Sorted by `normalized_SEI` from high to low.
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
#' # Compute boundary weights for cluster 0 cells
#' bon_c0 <- GetBoundary(data = coords, one_cluster = 0)
#' cells_c0 <- subset(coords, cluster==0)[,"cell"]
#' weights_bon <- ComputeBoundaryWeights(data = coords,
#'                                       cell_ids = cells_c0,
#'                                       boundary = bon_c0)
#'
#' # Compute spatial enrichment index for boundaries
#' library(Matrix)
#' df_sei <- ComputeSpatialEnrichmentIndex(exp_mat = logNorm_expr[,cells_c0],
#'                                         weights = weights_bon)
#' head(df_sei)
#'

ComputeSpatialEnrichmentIndex <- function(exp_mat, weights) {
  # --- Input checks ---
  if (!inherits(exp_mat, c("matrix", "dgCMatrix"))) {
    stop("'exp_mat' must be of class 'matrix' or 'dgCMatrix'")
  }
  if (!is.numeric(weights)) {
    stop("'weights' must be a numeric vector.")
  }
  if (length(weights) != ncol(exp_mat)) {
    stop("Length of 'weights' must match number of columns in 'exp_mat'")
  }

  # --- Compute weighted mean expression (SEI) ---
  if (inherits(exp_mat, "dgCMatrix")) {
    weighted_expr <- exp_mat %*% Matrix::Diagonal(x = weights)
    SEI_scores <- Matrix::rowSums(weighted_expr) / sum(weights)
    mean_expr <- Matrix::rowMeans(exp_mat)
  } else {
    weighted_expr <- sweep(exp_mat, 2, weights, `*`)
    SEI_scores <- rowSums(weighted_expr) / sum(weights)
    mean_expr <- rowMeans(exp_mat)
  }

  # --- Compute normalized SEI ---
  normalized_SEI <- SEI_scores / (mean_expr + 1e-6)

  # --- Assemble and sort result ---
  df <- data.frame(
    gene = rownames(exp_mat),
    SEI = SEI_scores,
    mean_expr = mean_expr,
    normalized_SEI = normalized_SEI,
    row.names = NULL
  )

  df <- df[order(-df$normalized_SEI), ]
  rownames(df) <- NULL

  return(df)
}


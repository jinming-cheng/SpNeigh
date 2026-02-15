#' Compute Spatial Enrichment Index (SEI) for Each Gene
#'
#' Calculates the spatial enrichment index (SEI) for each gene based on
#' a user-supplied set of spatial weights. The SEI reflects the extent to
#' which gene expression is enriched in spatially weighted regions of
#' the tissue (e.g., near boundaries or centroids).
#'
#' This method supports both dense (`matrix`) and sparse (`dgCMatrix`)
#' gene expression formats, and can be applied using any distance-based
#' weighting scheme.
#'
#' @param exp_mat A normalized gene expression matrix with genes as rows and
#'                cells as columns. Should be of class `matrix` or `dgCMatrix`.
#' @param weights A numeric vector of spatial weights (e.g.,
#'                from `computeBoundaryWeights` or `computeCentroidWeights`).
#'                Must be the same length as the number of columns (cells)
#'                in `exp_mat`.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{`gene`}{Gene name}
#'   \item{`SEI`}{Spatial enrichment index: weighted mean expression
#'                across cells}
#'   \item{`mean_expr`}{Mean expression across all cells (unweighted)}
#'   \item{`normalized_SEI`}{Ratio of SEI to mean expression;
#'                           used to compare genes independent of
#'                           baseline expression}
#' }
#' The result is sorted in descending order by `normalized_SEI`.
#'
#' @export
#'
#' @examples
#' # Load spatial coordinates and log-normalized expression
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'     package = "SpNeigh"
#' ))
#' logNorm_expr <- readRDS(system.file("extdata", "LogNormExpr.rds",
#'     package = "SpNeigh"
#' ))
#'
#' # Compute spatial weights and SEI
#' bon_c0 <- getBoundary(data = coords, one_cluster = 0)
#' cells_c0 <- subset(coords, cluster == 0)$cell
#' weights <- computeBoundaryWeights(
#'     data = coords,
#'     cell_ids = cells_c0,
#'     boundary = bon_c0
#' )
#' sei_df <- computeSpatialEnrichmentIndex(logNorm_expr[, cells_c0], weights)
#' head(sei_df)
#'
computeSpatialEnrichmentIndex <- function(exp_mat = NULL, weights = NULL) {
    # --- Input checks ---
    if (!inherits(exp_mat, c("matrix", "dgCMatrix"))) {
        stop("'exp_mat' must be of class 'matrix' or 'dgCMatrix'.")
    }

    if (!is.numeric(weights)) {
        stop("'weights' must be a numeric vector.")
    }

    if (length(weights) != ncol(exp_mat)) {
        stop(
            "Length of 'weights' must match the number ",
            "of columns (cells) in 'exp_mat'."
        )
    }

    if (sum(weights) == 0) {
        stop("'weights' must not sum to zero.")
    }

    # --- Ensure gene names exist ---
    gene_names <- rownames(exp_mat)
    if (is.null(gene_names)) {
        gene_names <- paste0("gene_", seq_len(nrow(exp_mat)))
    }

    # --- Normalize weights ---
    weights <- weights / sum(weights)

    # --- Compute SEI ---
    SEI_scores <- as.numeric(exp_mat %*% weights)
    mean_expr <- rowMeans(exp_mat)

    # --- Normalize SEI (avoid division by zero) ---
    normalized_SEI <- SEI_scores / (mean_expr + 1e-6)

    # --- Assemble result ---
    df <- data.frame(
        gene = gene_names,
        SEI = SEI_scores,
        mean_expr = mean_expr,
        normalized_SEI = normalized_SEI,
        stringsAsFactors = FALSE
    )

    df <- df[order(df$normalized_SEI, decreasing = TRUE), ]
    rownames(df) <- NULL

    return(df)
}



#' @title Compute Spatial Enrichment Index (SEI)
#' @description Alias for \code{\link{computeSpatialEnrichmentIndex}}.
#' @inheritParams computeSpatialEnrichmentIndex
#' @return A data frame containing the spatial enrichment index (SEI) results.
#' @export
#' @examples
#' # Load spatial coordinates and log-normalized expression
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'     package = "SpNeigh"
#' ))
#' logNorm_expr <- readRDS(system.file("extdata", "LogNormExpr.rds",
#'     package = "SpNeigh"
#' ))
#'
#' # Compute spatial weights and SEI
#' bon_c0 <- getBoundary(data = coords, one_cluster = 0)
#' cells_c0 <- subset(coords, cluster == 0)$cell
#' weights <- computeBoundaryWeights(
#'     data = coords,
#'     cell_ids = cells_c0,
#'     boundary = bon_c0
#' )
#' sei_df <- computeSEI(logNorm_expr[, cells_c0], weights)
#' head(sei_df)
computeSEI <- computeSpatialEnrichmentIndex

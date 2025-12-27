#' Compute a spatial neighborhood interaction matrix using
#' K-nearest neighbors (KNN)
#'
#' Computes a spatial interaction matrix where each entry quantifies
#' the number of neighboring cells from a given cluster (columns) that are
#' among the `k`-nearest neighbors of cells in another cluster (rows).
#' This provides a summary of spatial proximity and enrichment of
#' neighboring clusters around each focal cluster.
#'
#' The matrix is built by identifying the `k` nearest neighbors for each cell
#' based on spatial coordinates, and then tabulating the cluster identities
#' of those neighbors with respect to the cluster identity of the focal cell.
#'
#' @inheritParams ExtractCoords
#' @param k Integer. Number of nearest neighbors to use for each cell.
#'          Default is 10.
#'
#' @return A numeric matrix where rows represent focal clusters and columns
#'         represent neighboring clusters. Each cell in the matrix indicates
#'         how frequently a neighbor cluster appears among the k-nearest
#'         neighbors of cells from the focal cluster.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' # Load coordinates
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'     package = "SpNeigh"
#' ))
#'
#' # Compute interaction matrix using all cells
#' interaction_matrix <- ComputeSpatialInteractionMatrix(coords)
#' head(interaction_matrix)
#'
#' # Compute interaction matrix for cells inside boundaries
#' boundary_points <- GetBoundary(
#'     data = coords, one_cluster = 2,
#'     eps = 120, minPts = 10
#' )
#' cells_inside <- GetCellsInside(data = coords, boundary = boundary_points)
#' coords_sub <- subset(coords, cell %in% cells_inside$cell)
#' ComputeSpatialInteractionMatrix(coords_sub)
#'
#' # Compute interaction matrix for cells inside ring region 2
#' ring_regions <- GetRingRegion(boundary = boundary_points, dist = 100)
#' cells_ring <- GetCellsInside(data = coords, boundary = ring_regions[2, ])
#' coords_sub <- subset(coords, cell %in% cells_ring$cell)
#' ComputeSpatialInteractionMatrix(coords_sub)
#'
ComputeSpatialInteractionMatrix <- function(data = NULL,
                                            cluster_col = NULL,
                                            k = 10) {
    # --- Extract coordinates from data ---
    sp_coords <- ExtractCoords(data, cluster_col = cluster_col)

    # --- Build KNN graph (cells x neighbors) ---
    knn <- FNN::get.knn(sp_coords[, c("x", "y")], k = k)

    # Build long-format neighbor dataframe
    knn_df <- lapply(seq_len(nrow(knn$nn.index)), function(i) {
        neighbors <- knn$nn.index[i, ]
        data.frame(
            cell = sp_coords$cell[i],
            neighbor = sp_coords$cell[neighbors]
        )
    }) %>% dplyr::bind_rows()

    # Annotate with cluster identities
    neighbor_annot <- knn_df %>%
        dplyr::left_join(sp_coords, by = "cell") %>%
        dplyr::rename(cell_cluster = "cluster") %>%
        dplyr::left_join(sp_coords, by = c("neighbor" = "cell")) %>%
        dplyr::rename(neighbor_cluster = "cluster")

    # --- Build interaction matrix ---
    interaction_counts <- neighbor_annot %>%
        dplyr::count(.data$cell_cluster, .data$neighbor_cluster, name = "count")

    interaction_matrix <- interaction_counts %>%
        tidyr::pivot_wider(
            names_from = "neighbor_cluster",
            values_from = "count",
            values_fill = 0
        ) %>%
        tibble::column_to_rownames("cell_cluster")

    return(as.matrix(interaction_matrix))
}

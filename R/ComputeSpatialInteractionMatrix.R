
#' Compute a spatial neighborhood interaction matrix based on K-nearest neighbors algorithm
#'
#' The spatial neighborhood interaction matrix is computed based on counting the number of clusters of k-nearest neighbor cells by the cluster a focal cell.
#' @inheritParams ExtractCoords
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @param k Number of nearest neighbors to use. Default is 10.
#'
#' @return A matrix where rows are focal clusters and columns are neighbor clusters.
#' @export
#' @examples
#' # Load coordinates
#' load(system.file("extdata", "MouseBrainTinyCoords.rda",
#'                  package = "SpNeigh"))
#' head(coords)
#'
#' # Compute KNN spatial interaction matrix of all cells
#' interaction_matrix <- ComputeSpatialInteractionMatrix(coords)
#' interaction_matrix
#'
#' # Compute KNN spatial interaction matrix of cells insides boundaries
#' boundary_points <- GetBoundary(data = coords, one_cluster = 2,
#'                                eps = 120, minPts = 10)
#' cells_inside <- GetCellsInside(data = coords, boundary =  boundary_points)
#' coords_sub <- subset(coords, cell %in% cells_inside$cell)
#' interaction_matrix <- ComputeSpatialInteractionMatrix(coords_sub)
#' interaction_matrix
#'
#' # Compute KNN spatial interaction matrix of cells insides rings
#' ring_regions <- GetRingRegion(boundary = boundary_points, dist = 100)
#' cells_ring <- GetCellsInside(data = coords, boundary =  ring_regions)
#' coords_sub <- subset(coords, cell %in% cells_ring$cell)
#' interaction_matrix <- ComputeSpatialInteractionMatrix(coords_sub)
#' interaction_matrix

ComputeSpatialInteractionMatrix <- function(data = NULL, k = 10){

  # extract coordinates from data
  sp_coords <- ExtractCoords(data)

  #  Build KNN graph
  knn <- FNN::get.knn(sp_coords[, c("x", "y")], k = k)
  knn_df <- lapply(1:nrow(knn$nn.index), function(i) {
      neighbors <- knn$nn.index[i, ]
      data.frame( cell = sp_coords$cell[i],
                  neighbor = sp_coords$cell[neighbors]
                  )
      }) %>% dplyr::bind_rows()

  # Annotate neighbor clusters
  neighbor_annot <- knn_df %>%
    dplyr::left_join(sp_coords, by = c("cell" = "cell")) %>%
    dplyr::rename(cell_cluster = .data$cluster) %>%
    dplyr::left_join(sp_coords, by = c("neighbor" = "cell")) %>%
    dplyr::rename(neighbor_cluster = .data$cluster)

  # Build interaction matrix
  interaction_counts <- neighbor_annot %>%
    dplyr::count(.data$cell_cluster, .data$neighbor_cluster, name = "count")

  interaction_matrix <- interaction_counts %>%
    tidyr::pivot_wider(names_from = .data$neighbor_cluster,
                         values_from = .data$count, values_fill = 0) %>%
      tibble::column_to_rownames("cell_cluster")

  as.matrix(interaction_matrix)
}

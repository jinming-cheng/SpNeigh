
#' Make a factor of natural order
#' @param x A vector
#' @export
#' @examples
#' FactorNaturalOrder(10:1)
#' FactorNaturalOrder( c("a11", "a12", "a1","a2", "a") )
#'
FactorNaturalOrder <- function(x) {
  x <- as.character(x)
  ord <- stringr::str_order(x, numeric = TRUE)
  factor(x, levels = unique(x[ord]))
}

#' Heatmap of row-scaled interaction matrix
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @param interaction_matrix A numeric matrix where rows represent focal clusters and columns represent neighbor clusters. For example, the output `interaction_matrix` from `ComputeSpatialInteractionMatrix`.
#' @param low_color Color for low values. Default is blue.
#' @param mid_color Color for mid-range values. Default is white.
#' @param high_color Color for high values. Default: red.
#' @param angle_x_label Angle to rotate the x-axis labels.
#' @param title Title for the heatmap.
#' @export
#' @examples
#' # Load coordinates
#' load(system.file("extdata", "MouseBrainTinyCoords.rda",
#'                  package = "SpNeigh"))
#' head(coords)
#'
#'
#' # Compute KNN spatial interaction matrix of cells insides boundaries
#' boundary_points <- GetBoundary(data = coords, one_cluster = 2,
#'                                eps = 120, minPts = 10)
#' ring_regions <- GetRingRegion(boundary = boundary_points, dist = 100)
#' cells_ring <- GetCellsInside(data = coords, boundary =  ring_regions)
#' coords_sub <- subset(coords, cell %in% cells_ring$cell)
#' interaction_matrix <- ComputeSpatialInteractionMatrix(coords_sub)
#'
#' # Plot interaction matrix
#' PlotInteractionMatrix(interaction_matrix)

PlotInteractionMatrix <- function(interaction_matrix,
                                  low_color = "blue",
                                  mid_color = "white",
                                  high_color = "red",
                                  angle_x_label = 45,
                                  title = "Row-Scaled Interaction Matrix (Z-scores)") {
  # Ensure it's a matrix
  if (!is.matrix(interaction_matrix)) {
    stop("interaction_matrix must be a numeric matrix.")
  }

  # Row-scale (z-score per row)
  scaled_mat <- t(scale(t(interaction_matrix)))

  # Convert to long format for ggplot
  df_long <- as.data.frame(scaled_mat) %>%
    tibble::rownames_to_column("Focal") %>%
    tidyr::pivot_longer(-.data$Focal, names_to = "Neighbor", values_to = "Zscore")

  # Make Focal and Neigbor clusters factors and with natural order levels
  df_long$Focal <- FactorNaturalOrder(df_long$Focal)
  df_long$Neighbor <- FactorNaturalOrder(df_long$Neighbor)

  # Plot with ggplot
  ggplot2::ggplot(df_long, ggplot2::aes(x = .data$Neighbor, y = .data$Focal,
                                        fill = .data$Zscore)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(low = low_color, mid = mid_color,
                                  high = high_color, midpoint = 0) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = angle_x_label, hjust = 1),
                   plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                   panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(title = title, x = "Neighbor Cluster",
                  y = "Focal Cluster", fill = "Z-score")
}

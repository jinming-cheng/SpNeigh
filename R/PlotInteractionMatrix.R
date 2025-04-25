

#' Plot a heatmap of a row-scaled spatial interaction matrix
#'
#' Visualizes a spatial interaction matrix using a heatmap, where rows represent focal clusters
#' and columns represent neighbor clusters. Each row is scaled using z-scores to highlight
#' relative enrichment patterns across neighbor types. This is useful for detecting
#' spatial proximity patterns between cell populations.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @param interaction_matrix A numeric matrix with focal clusters as rows and neighbor clusters as columns.
#'        Typically the output from `ComputeSpatialInteractionMatrix()`.
#' @param low_color Color representing low z-score values. Default is `"blue"`.
#' @param mid_color Color representing the midpoint (z-score = 0). Default is `"white"`.
#' @param high_color Color representing high z-score values. Default is `"red"`.
#' @param angle_x_label Angle (in degrees) to rotate x-axis labels. Default is `45`.
#' @param title Title for the heatmap.
#'
#' @return A `ggplot` object representing the row-scaled heatmap of the interaction matrix.
#'
#' @export
#' @examples
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'                               package = "SpNeigh"))
#'
#' boundary_points <- GetBoundary(data = coords, one_cluster = 2,
#'                                eps = 120, minPts = 10)
#' ring_regions <- GetRingRegion(boundary = boundary_points, dist = 100)
#' cells_ring <- GetCellsInside(data = coords, boundary = ring_regions)
#' coords_sub <- subset(coords, cell %in% cells_ring$cell)
#' interaction_matrix <- ComputeSpatialInteractionMatrix(coords_sub)
#'
#' PlotInteractionMatrix(interaction_matrix)
#'

PlotInteractionMatrix <- function(interaction_matrix,
                                  low_color = "blue",
                                  mid_color = "white",
                                  high_color = "red",
                                  angle_x_label = 45,
                                  title = "Row-Scaled Interaction Matrix (Z-scores)") {
  # Ensure it's a matrix
  if (!is.matrix(interaction_matrix)) {
    stop("`interaction_matrix` must be a numeric matrix.")
  }

  # Row-scale (z-score per row)
  scaled_mat <- t(scale(t(interaction_matrix)))

  # Convert to long format for ggplot
  df_long <- as.data.frame(scaled_mat) %>%
    tibble::rownames_to_column("Focal") %>%
    tidyr::pivot_longer(-.data$Focal, names_to = "Neighbor", values_to = "Zscore")

  # Keep raw orders
  df_long$Focal <- factor(df_long$Focal, levels = rownames(interaction_matrix))
  df_long$Neighbor <- factor(df_long$Neighbor, levels = colnames(interaction_matrix))

  # Plot
  ggplot2::ggplot(df_long, ggplot2::aes(x = .data$Neighbor, y = .data$Focal,
                                        fill = .data$Zscore)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(low = low_color, mid = mid_color,
                                  high = high_color, midpoint = 0) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = angle_x_label, hjust = 1),
                   plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                   panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(title = title,
                  x = "Neighbor Cluster",
                  y = "Focal Cluster",
                  fill = "Z-score")
}

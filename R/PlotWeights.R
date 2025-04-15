
#' Plot spatial weights for cells
#'
#' Plots spatial coordinates of cells colored by weight values. Only cells in the `data` with weights will be plotted.
#'
#' @inheritParams ComputeCentroidWeights
#' @importFrom rlang .data
#' @param weights A named numeric vector where names are cell IDs.
#' @param point_size Point size in the plot. Default is 0.2.
#' @param theme_ggplot Theme for ggplot.
#'
#' @export
#' @examples
#' # Load coordinates
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'                               package = "SpNeigh"))
#' head(coords)
#'
#' # Obtain cells in cluster 2 and focus on these cells
#' cells_c2 <- subset(coords, cluster==2)[,"cell"]
#'
#' # Obtain and plot centroid weights
#' weights_cen <- ComputeCentroidWeights(data = coords, cell_ids = cells_c2)
#' PlotWeights(data = coords, weights = weights_cen)
#'
#' # Obtain and plot boundary weights
#' boundary_points <- GetBoundary(data = coords, one_cluster = 2)
#' weights_bon <- ComputeBoundaryWeights(data = coords, cell_ids = cells_c2,
#'                                       boundary = boundary_points)
#' PlotWeights(data = coords, weights = weights_bon)
#'
PlotWeights <- function(data = NULL,
                        weights = NULL,
                        point_size = 0.2,
                        theme_ggplot = my_theme_ggplot()) {
  # Extract coordinates from data
  if( inherits(data,"Seurat") ){
    sp_coords <- Seurat::GetTissueCoordinates(data)
  }else{
    sp_coords <- data
  }

  if (!all(c("cell", "x", "y") %in% colnames(sp_coords))) {
    stop("data must contain columns: 'cell', 'x', and 'y'")
  }

  if (is.null(names(weights))) {
    stop("weights must be a named numeric vector with cell IDs as names")
  }

  # Match weights to sp_coords using cell IDs
  matched_idx <- match(sp_coords$cell, names(weights))
  matched_weights <- weights[matched_idx]

  # Remove unmatched cells (NA weights)
  keep <- !is.na(matched_weights)
  plot_df <- sp_coords[keep, ]
  plot_df$weights <- matched_weights[keep]

  if (nrow(plot_df) == 0) {
    stop("No matching cell IDs between 'data' and 'weights'")
  }

  # Plot
  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$x, y = .data$y, color = weights)) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_color_viridis_c(option = "plasma", name = "Weight") +
    theme_ggplot
}


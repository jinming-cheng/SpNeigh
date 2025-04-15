
#' Plot cells inside boundries or rings
#'
#' @inheritParams PlotBoundary
#' @importFrom rlang .data
#' @param cells_inside A sf object. `cells_inside` is returned by `GetCellsInside`.
#' @param fixed_aspect_ratio Logical. Whether to plot using `sf::geom_sf()` function with 1:1 aspect ratio. Default is TRUE
#' @export
#' @examples
#' # Load coordinates
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'                               package = "SpNeigh"))
#' head(coords)
#'
#' # Get cells inside boundaries
#' boundary_points <- GetBoundary(data = coords, one_cluster = 2,
#'                                eps = 120, minPts = 10)
#' cells_inside <- GetCellsInside(data = coords, boundary =  boundary_points)
#' PlotCellsInside(cells_inside = cells_inside)
#'
#' # Get cells inside rings
#' ring_regions <- GetRingRegion(boundary = boundary_points, dist = 100)
#' cells_ring <- GetCellsInside(data = coords, boundary =  ring_regions)
#' PlotCellsInside(cells_inside = cells_ring)
#'
PlotCellsInside <- function(cells_inside = NULL,
                            point_size = 0.5,
                            colors = my_colors_15,
                            theme_ggplot = my_theme_ggplot(),
                            legend_size = 2,
                            fixed_aspect_ratio = TRUE){

  # Make sure cluster is a factor
  if(is.null(levels(cells_inside$cluster))){
    cells_inside$cluster <- FactorNaturalOrder(cells_inside$cluster)
  }

  # Add names to colors
  named_colors <- colors[1:nlevels(cells_inside$cluster)]
  names(named_colors) <- levels(cells_inside$cluster)

  # Convert the sf object into a data frame
  coords <- sf::st_coordinates(cells_inside)
  colnames(coords) <- c("x","y")
  df <- sf::st_drop_geometry(cells_inside)
  df <- cbind(df, coords)
  df$cluster <- factor(df$cluster, levels = levels(cells_inside$cluster) )

  # Plot cells colored by cluster
  if(fixed_aspect_ratio){
    p <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = cells_inside, size = point_size,
                       ggplot2::aes(color = .data$cluster))
  }else{
    p <- ggplot2::ggplot() +
      ggplot2::geom_point(data = df, size = point_size,
                          ggplot2::aes(x = .data$x, y = .data$y,
                                       color = .data$cluster))
  }

  # Add theme and colors
  p <- p +  theme_ggplot +
    ggplot2::scale_color_manual(values = named_colors) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=legend_size)))

  p
}


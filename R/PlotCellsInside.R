
#' Plot cells inside boundries or rings
#'
#' @inheritParams PlotBoundary
#' @importFrom rlang .data
#' @param cells_inside A sf object. cells_inside returned by `GetCellsInside`.
#' @param ... Other paramters in `sf::geom_sf()` function.
#' @export
#' @examples
#' # Load coordinates
#' load(system.file("extdata", "MouseBrainTinyCoords.rda",
#'                  package = "SpNeigh"))
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
                            ...){
  # add names to colors
  named_colors <- colors[1:nlevels(cells_inside$cluster)]
  names(named_colors) <- levels(cells_inside$cluster)

  ggplot2::ggplot() +
    ggplot2::geom_sf(data = cells_inside,
                     ggplot2::aes(color = .data$cluster), size = point_size, ...) +
    theme_ggplot +
    ggplot2::scale_color_manual(values = named_colors) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=legend_size)))
}


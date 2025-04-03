
#' Plot cells inside boundries or rings
#'
#' @inheritParams PlotBoundary
#' @importFrom rlang .data
#' @param cells_inside A sf object. cells_inside returned by `GetCellsInside`.
#' @param ... Other paramters in `sf::geom_sf()` function.
#' @export
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


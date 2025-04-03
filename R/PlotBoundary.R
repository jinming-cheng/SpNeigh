
#' A Custom ggplot2 Theme
#'
#' This is a ggplot2 theme defined and used in this package.
#'
#' @return A ggplot2 theme object.
#' @export
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(mpg, wt)) +
#'   geom_point() +
#'   my_theme_ggplot()

my_theme_ggplot <- function() {
  ggplot2::theme_classic() +
    ggplot2::theme(axis.text = ggplot2::element_text(size=12),
                   axis.title = ggplot2::element_text(size=15,face="bold"),
                   strip.text = ggplot2::element_text(size=12, face="bold"),
                   plot.title = ggplot2::element_text(size=15,face="bold",hjust=0.5),
                   plot.margin = ggplot2::margin(0.5, 0.5, 0.5, 0.5, "cm"))}



#' Draw boundary for a cluster or cell population
#'
#' @inheritParams GetBoundary
#' @inheritParams ExtractCoords
#' @importFrom rlang .data
#' @param one_cluster Specify one cluster to obtain its boundary. If `boundary` parameter is provided, it can be omitted. However, it must be provided, if `sub_plot` is TRUE.
#' @param boundary Boundary obtained from `GetBoundary` function. If it is NULL, get the boundary of the specified cluster automatically. If both `boundary` and `one_cluster` are NULL, plot coordinates only.
#' @param colors A vector of colors used in plotting for clusters.
#' @param point_size Point size of cells. Default is 0.5.
#' @param color_boundary Color of the boundary for a cluster. Default is black.
#' @param linewidth_boundary Linewidth of the boundary for a cluster. Default is 1.5.
#' @param sub_plot A logical value indicates whether to only plot cells in one cluster. If FALSE, plot all cells in all clusters. Default is FALSE.
#' @param theme_ggplot Theme for ggplot.
#' @param legend_size Legend size for ggplot.
#' @param ... Parameters in `GetBoundary` function.
#' @export
PlotBoundary <- function(data = NULL,
                         one_cluster = NULL,
                         boundary = NULL,
                         colors = my_colors_15,
                         point_size = 0.5,
                         color_boundary = "black",
                         linewidth_boundary = 1.5,
                         sub_plot = FALSE,
                         theme_ggplot = my_theme_ggplot(),
                         legend_size = 2,
                         ...){

  # extract coordinates from data
  sp_coords <- ExtractCoords(data = data)


  # if the boundary of a cluster is not provided, get the boundary automatically
  if( is.null(boundary) & (!is.null(one_cluster)) ){
    boundary <- GetBoundary(data = sp_coords,one_cluster = one_cluster, ...)
  }

  # Whether to plot one cluster or all clusters
  if(sub_plot){
    # Obtain coordinates of the target cluster
    coords_sub <- dplyr::filter(sp_coords,.data$cluster == one_cluster)
    data_for_plot = coords_sub
  }else{
    data_for_plot = sp_coords
  }

  # add names to colors
  named_colors <- colors[1:nlevels(sp_coords$cluster)]
  names(named_colors) <- levels(sp_coords$cluster)

  # Plot boundary
  p <- ggplot2::ggplot(data_for_plot, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_point(ggplot2::aes(color = .data$cluster), alpha = 1, size = point_size) +
    theme_ggplot +
    ggplot2::scale_color_manual(values = named_colors) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=legend_size)))

  # whether to plot boundary
  if( !is.null(boundary) ){
    p <- p +
      ggplot2::geom_polygon(data = boundary,
                            ggplot2::aes(x = .data$x, y = .data$y, group = .data$region_id),
                            fill = NA,
                            color = color_boundary,
                            linewidth = linewidth_boundary)
  }

  p

}

#' Add a boundary or boundaries to the plot
#'
#' @inheritParams PlotBoundary
#' @importFrom rlang .data
#' @param boundary Boundary obtained from `GetBoundary` function.
#' @export
#'
AddBoundary <- function(boundary = NULL, color_boundary="black", linewidth_boundary=1.5) {
  ggplot2::geom_polygon(data = boundary,
                        ggplot2::aes(x = .data$x, y = .data$y, group = .data$region_id),
                        fill = NA,
                        color = color_boundary,
                        linewidth = linewidth_boundary)

}

#' Add a boundary poly to the plot
#'
#' @inheritParams PlotBoundary
#' @importFrom rlang .data
#' @param boundary_poly A sf object returned by `BuildBoundaryPoly()`, `GetOuterBoundary()` or `GetRingRegion()`.
#' @export
#'
AddBoundaryPoly <- function(boundary_poly, color_boundary="black", linewidth_boundary=1.5) {
  ggplot2::geom_sf(data = boundary_poly,
                   inherit.aes = FALSE,
                   fill=NA,
                   color = color_boundary,
                   linewidth = linewidth_boundary)

}



#' Plot regions inside boundries or rings
#'
#' @inheritParams PlotBoundary
#' @inheritParams AddBoundaryPoly
#' @importFrom rlang .data
#' @param alpha The transparency for the fill aesthetics of a geom in ggplot.
#' @param color_boundary Color of the boundaries. Default is black.
#' @param linewidth_boundary Linewidth of the boundaries. Default is 1.
#' @param ... Other paramters in `sf::geom_sf()` function.
#' @export
PlotRegion <- function(boundary_poly = NULL,
                       alpha = 0.5,
                       color_boundary="black",
                       linewidth_boundary = 1,
                       theme_ggplot = my_theme_ggplot(),
                       ...){
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = boundary_poly, ggplot2::aes(fill = as.factor(.data$region_id)),
            alpha = alpha, color = color_boundary, linewidth = linewidth_boundary, ...) +
    ggplot2::labs(fill = "Region ID") +
    theme_ggplot
}


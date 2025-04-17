
#' A custom ggplot2 theme
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
#'
#' @inheritParams GetBoundary
#' @inheritParams ExtractCoords
#' @importFrom rlang .data
#' @param one_cluster Specify one cluster to obtain its boundary. If `boundary` parameter is provided, it can be omitted. However, it must be provided, if `sub_plot` is TRUE.
#' @param boundary A data frame or sf object containing x, y and region_id information. If it is NULL, get the boundary of the specified cluster automatically. If both `boundary` and `one_cluster` are NULL, plot coordinates only.
#' @param colors A character vector specifying the colors for clusters.
#' @param point_size Point size of cells. Default is 0.5.
#' @param color_boundary Color of the boundary for a cluster. Default is black.
#' @param linewidth_boundary Linewidth of the boundary borders. Default is 1.5.
#' @param sub_plot Logical. Whether to plot cells in one cluster. If FALSE (default), plot all cells in all clusters.
#' @param theme_ggplot Theme for ggplot.
#' @param legend_size Legend size for ggplot.
#' @param ... Parameters in `GetBoundary` function.
#' @export
#' @examples
#' # Load coordinates
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'                               package = "SpNeigh"))
#' head(coords)
#'
#' # Plot coordinates without boundary
#' PlotBoundary(coords)
#'
#' # Plot boundaries of cluster 2
#' PlotBoundary(coords, one_cluster = 2)
#'
#' # Plot boundaries of cluster 2, and adjust subregion numbers
#' # for dbscan method using eps and minPts parameters
#' PlotBoundary(coords, one_cluster = 2,
#'             multi_region = TRUE,
#'             subregion_method = "dbscan",
#'             eps = 120, minPts = 10)
#'
#' # Alternatively, the boundaries can be manually specified
#' boundary_points <- GetBoundary(data = coords, one_cluster = 2,
#'                                subregion_method = "dbscan",
#'                                eps = 120, minPts = 10)
#' PlotBoundary(data = coords, boundary = boundary_points)

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

  # Extract coordinates from data
  sp_coords <- ExtractCoords(data = data)

  # If the boundary of a cluster is not provided, get the boundary automatically
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
    p <- p + AddBoundary(boundary = boundary,
                         color_boundary = color_boundary,
                         linewidth_boundary = linewidth_boundary)
  }

  p

}

#' Add a boundary or boundaries to the plot
#'
#' @inheritParams PlotBoundary
#' @importFrom rlang .data
#' @param boundary A data frame or sf object containing x, y and region_id information.
#' @export
#' @examples
#' # Load coordinates
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'                               package = "SpNeigh"))
#' head(coords)
#'
#' # Add boundaries to an existing plot
#' boundary_points <- GetBoundary(data = coords, one_cluster = 2,
#'                                subregion_method = "dbscan",
#'                                eps = 120, minPts = 10)
#' PlotBoundary(coords) + AddBoundary(boundary_points)
#'
AddBoundary <- function(boundary = NULL, color_boundary="black", linewidth_boundary=1.5) {

  # Check boundary format
  if(!inherits(boundary,c("data.frame","sf"))){
    stop("`boundary` should be a data.frame or an sf object.")
  }

  # Convert boundary polys to boundary points
  if(inherits(boundary,"sf")){
    boundary <- BoundaryPolyToPoints(boundary)
  }

  # Check required columns in boundary
  if( !all(c("x","y","region_id") %in% colnames(boundary)) ){
    stop("`boundary` must contain columns: 'x', 'y' and 'region_id'.")
  }

  ggplot2::geom_path(data = boundary,
                     ggplot2::aes(x = .data$x, y = .data$y, group = .data$region_id),
                     color = color_boundary,
                     linewidth = linewidth_boundary)

}

#' Add a boundary poly to the plot
#'
#' @inheritParams PlotBoundary
#' @importFrom rlang .data
#' @param boundary_poly An sf object containing x, y and region_id information.
#' @export
#' @examples
#' # Load coordinates
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'                               package = "SpNeigh"))
#' head(coords)
#'
#' # Build boundary polygons from the boundary points
#' boundary_points <- GetBoundary(data = coords, one_cluster = 2,
#'                                subregion_method = "dbscan",
#'                                eps = 120, minPts = 10)
#' boundary_polys <- BuildBoundaryPoly(boundary_points)
#'
#' # Add boundary polygons to an existing plot
#' PlotBoundary(coords) + AddBoundaryPoly(boundary_polys,
#'                                        color="blue")
#'

AddBoundaryPoly <- function(boundary_poly, color_boundary="black", linewidth_boundary=1.5) {

  # Check boundary_poly
  if(!inherits(boundary_poly,c("sf"))){
    stop("`boundary_poly` should be an sf object.")
  }

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
#' @param ... Other parameters in `sf::geom_sf()` function.
#' @export
#' @examples
#' # Load coordinates
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'                               package = "SpNeigh"))
#' head(coords)
#'
#' # Plot region inside boundaries
#' boundary_points <- GetBoundary(data = coords, one_cluster = 2,
#'                                eps = 120, minPts = 10)
#' boundary_polys = BuildBoundaryPoly(boundary_points)
#' PlotRegion(boundary_poly = boundary_polys)
#'
#' # Plot region inside rings
#' ring_regions <- GetRingRegion(boundary = boundary_points, dist = 100)
#' PlotRegion(boundary_poly = ring_regions)
#'
PlotRegion <- function(boundary_poly = NULL,
                       alpha = 0.5,
                       color_boundary="black",
                       linewidth_boundary = 1,
                       theme_ggplot = my_theme_ggplot(),
                       ...){
  # Check boundary_poly
  if(!inherits(boundary_poly,c("sf"))){
    stop("`boundary_poly` should be an sf object.")
  }

  # Make region_id a factor of natural orders
  if(is.null(levels(boundary_poly$region_id))){
    boundary_poly$region_id <- FactorNaturalOrder(boundary_poly$region_id)
  }

  ggplot2::ggplot() +
    ggplot2::geom_sf(data = boundary_poly,
                     ggplot2::aes(fill = .data$region_id),
                     alpha = alpha,
                     color = color_boundary,
                     linewidth = linewidth_boundary, ...) +
    ggplot2::labs(fill = "Region ID") +
    theme_ggplot
}

#' Plot boundary edge
#'
#' @inheritParams PlotBoundary
#' @inheritParams AddBoundaryPoly
#' @importFrom rlang .data
#' @param linewidth_boundary Linewidth of the boundaries. Default is 1.
#' @param ... Other parameters in `sf::geom_sf()` function.
#' @export
#' @examples
#' # Load coordinates
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'                               package = "SpNeigh"))
#' head(coords)
#'
#' # Plot region inside boundaries
#' boundary_points <- GetBoundary(data = coords, one_cluster = 2,
#'                                eps = 120, minPts = 10)
#' boundary_polys = BuildBoundaryPoly(boundary_points)
#' PlotEdge(boundary_poly = boundary_polys)
#'
#' # Split boundary polygon 1 into two edges
#' boundary_edges <- SplitBoundaryPolyByAnchor(boundary_polys[1,])
#' PlotEdge(boundary_edges)
#'
PlotEdge <- function(boundary_poly = NULL,
                     linewidth_boundary = 1,
                     theme_ggplot = my_theme_ggplot(),
                     ...){
  # Check boundary_poly
  if(!inherits(boundary_poly,c("sf"))){
    stop("`boundary_poly` should be an sf object.")
  }

  # Make region_id a factor of natural orders
  if(is.null(levels(boundary_poly$region_id))){
    boundary_poly$region_id <- FactorNaturalOrder(boundary_poly$region_id)
  }

  ggplot2::ggplot() +
    ggplot2::geom_sf(data = boundary_poly,
                     ggplot2::aes(color = .data$region_id),
                     fill = NA,
                     linewidth = linewidth_boundary, ...) +
    ggplot2::labs(color = "Region ID") +
    theme_ggplot

}

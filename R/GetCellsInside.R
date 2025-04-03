
#' Get cells inside boundaries or rings
#'
#' @importFrom rlang .data
#' @inheritParams GetBoundary
#' @param boundary A sf object or a data frame. Boundary result returned by `GetBoundary` or `GetBoundaryPoly` functions. Or ring regions returned by `GetRingRegion`.
#' @export
#'
GetCellsInside <- function(data = NULL,
                           boundary = NULL){

  # extract coordinates from data
  sp_coords <- ExtractCoords(data = data)

  spatial_points <- sf::st_as_sf(x = sp_coords,coords = c("x", "y"), crs = NA)

  # Build safe polygons
  if(is(boundary,"sf")){
    boundary_polys <- boundary
  }else{
    boundary_polys <- BuildBoundaryPoly(boundary)
  }

  cells_in_poly <- sf::st_join(spatial_points, boundary_polys, join = sf::st_within)

  # Obtain cells inside the boundary
  cells_inside <- dplyr::filter(cells_in_poly,!is.na(.data$region_id))
  cells_inside
}

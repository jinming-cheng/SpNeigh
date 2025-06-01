#' Identify cells located within spatial boundaries or ring regions
#'
#' Returns the subset of cells from the input data that fall spatially
#' within a given boundary or ring.
#' The boundary can be provided either as raw boundary points
#' (from `GetBoundary()`), as polygons (from `BuildBoundaryPoly()`),
#' or as ring regions (from `GetRingRegion()`).
#' The function uses spatial point-in-polygon matching via `sf::st_within`.
#'
#' @importFrom rlang .data
#' @inheritParams GetBoundary
#' @param boundary An `sf` object (polygon or ring) or a data frame
#'                 of boundary points returned by `GetBoundary()`,
#'                 `BuildBoundaryPoly()`, or `GetRingRegion()`.
#'
#' @return A data frame (tibble) of cells located inside the
#'         given spatial region(s), with region assignment in
#'         a `region_id` column.
#' @export
#' @examples
#' # Load coordinates
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'     package = "SpNeigh"
#' ))
#' head(coords)
#'
#' # Get cells inside boundaries
#' boundary_points <- GetBoundary(
#'     data = coords, one_cluster = 2,
#'     eps = 120, minPts = 10
#' )
#' cells_inside <- GetCellsInside(data = coords, boundary = boundary_points)
#' cells_inside
#'
#' # Get cells inside rings
#' ring_regions <- GetRingRegion(boundary = boundary_points, dist = 100)
#' cells_ring <- GetCellsInside(data = coords, boundary = ring_regions)
#' cells_ring
#'
GetCellsInside <- function(data = NULL, boundary = NULL) {
    # Extract coordinate data
    sp_coords <- ExtractCoords(data = data)

    # Convert to sf point geometry
    spatial_points <- sf::st_as_sf(sp_coords, coords = c("x", "y"), crs = NA)

    # Ensure polygon format
    boundary_polys <- if (inherits(boundary, "sf")) {
        boundary
    } else {
        BuildBoundaryPoly(boundary)
    }

    # Spatial join: find cells inside polygons
    cells_in_poly <- sf::st_join(spatial_points, boundary_polys,
        join = sf::st_within
    )

    # Keep only cells that intersect with a region
    cells_inside <- dplyr::filter(cells_in_poly, !is.na(.data$region_id))

    return(cells_inside)
}

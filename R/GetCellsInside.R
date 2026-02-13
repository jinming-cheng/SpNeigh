#' Identify cells located within spatial boundaries or ring regions
#'
#' Returns the subset of cells from the input data that fall spatially
#' within a given boundary or ring.
#' The boundary can be provided either as raw boundary points
#' (from `getBoundary()`), as polygons (from `buildBoundaryPoly()`),
#' or as ring regions (from `getRingRegion()`).
#' The function uses spatial point-in-polygon matching via `sf::st_within`.
#'
#' @importFrom rlang .data
#' @inheritParams getBoundary
#' @param boundary An `sf` object (polygon or ring) or a data frame
#'                 of boundary points returned by `getBoundary()`,
#'                 `buildBoundaryPoly()`, or `getRingRegion()`.
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
#' boundary_points <- getBoundary(
#'     data = coords, one_cluster = 2,
#'     eps = 120, minPts = 10
#' )
#' cells_inside <- getCellsInside(data = coords, boundary = boundary_points)
#' cells_inside
#'
#' # Get cells inside rings
#' ring_regions <- getRingRegion(boundary = boundary_points, dist = 100)
#' cells_ring <- getCellsInside(data = coords, boundary = ring_regions)
#' cells_ring
#'
getCellsInside <- function(data = NULL, cluster_col = NULL, boundary = NULL) {
    # Extract coordinate data
    sp_coords <- extractCoords(data = data, cluster_col = cluster_col)

    # Convert to sf point geometry
    spatial_points <- sf::st_as_sf(sp_coords, coords = c("x", "y"), crs = NA)

    # Ensure polygon format
    boundary_polys <- if (inherits(boundary, "sf")) {
        boundary
    } else {
        buildBoundaryPoly(boundary)
    }

    # Spatial join: find cells inside polygons
    cells_in_poly <- sf::st_join(spatial_points, boundary_polys,
        join = sf::st_within
    )

    # Keep only cells that intersect with a region
    cells_inside <- dplyr::filter(cells_in_poly, !is.na(.data$region_id))

    return(cells_inside)
}

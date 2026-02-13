#' Compute spatial weights based on distance to the centroid
#'
#' Computes spatial weights for a set of cells based on their Euclidean distance
#' to the centroid of the selected group of cells. This is useful for modeling
#' gradient-like spatial expression centered around a cluster or region.
#' Supports multiple decay methods and optional distance scaling.
#'
#' @importFrom methods is
#' @param data A data frame, Seurat object or SpatialExperiment object
#'             containing spatial coordinates.
#'             Must include columns: `cell`, `x`, and `y`.
#' @param cell_ids A character vector of cell IDs to use for weight computation.
#' @param scale Logical. Whether to scale distances to the range \[0, 1\]
#'              before applying decay functions. Default is `TRUE`.
#' @param method Decay function to convert distances to weights.
#'               One of: "inverse", "gaussian", "linear", or "quadratic".
#' @param sigma Standard deviation for the Gaussian decay
#'              (used only if `method = "gaussian"`).
#'              Default is `0.5` (recommended for scaled distances).
#'
#' @return A named numeric vector of weights
#'         (length = number of cells in `cell_ids`).
#'         Higher weights correspond to cells closer to the centroid.
#'
#' @export
#' @examples
#' # Load spatial coordinates
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'     package = "SpNeigh"
#' ))
#'
#' # Select cells from cluster 2
#' cells_c2 <- subset(coords, cluster == 2)$cell
#'
#' # Compute centroid-based weights using default settings
#' weights <- computeCentroidWeights(data = coords, cell_ids = cells_c2)
#' head(weights)
#'
computeCentroidWeights <- function(
    data = NULL,
    cell_ids = NULL,
    scale = TRUE,
    method = c("inverse", "gaussian", "linear", "quadratic"),
    sigma = 0.5) {
    method <- match.arg(method)

    # Extract coordinates from data
    if (is(data, "Seurat")) {
        sp_coords <- extractCoords(data, extract_cluster = FALSE)
    } else {
        sp_coords <- data
    }

    if (!all(c("cell", "x", "y") %in% colnames(sp_coords))) {
        stop("data must contain columns: 'cell', 'x', and 'y'")
    }

    # Subset coordinates
    sub_coords <- sp_coords[sp_coords$cell %in% cell_ids, ]
    if (nrow(sub_coords) == 0) {
        stop("No matching cells found in 'data' for 'cell_ids'")
    }

    # Compute centroid of the selected cells
    centroid_x <- mean(sub_coords$x)
    centroid_y <- mean(sub_coords$y)

    # Compute Euclidean distances to centroid
    dists <- sqrt((sub_coords$x - centroid_x)^2 + (sub_coords$y - centroid_y)^2)

    # Optional scaling to [0, 1]
    if (scale && max(dists) > 0) {
        dists <- (dists - min(dists)) / (max(dists) - min(dists))
    }

    # Compute weights based on selected method
    weights <- switch(method,
        inverse = 1 / (1 + dists),
        gaussian = exp(-dists^2 / (2 * sigma^2)),
        linear = 1 - dists,
        quadratic = (1 - dists)^2
    )

    # Ensure valid weights
    weights[dists > 1 & scale] <- 0
    weights[weights < 0] <- 0
    names(weights) <- sub_coords$cell

    return(weights)
}



#' Compute spatial weights based on distance to nearest boundary
#'
#' Computes spatial weights for a subset of cells based on their Euclidean
#' distance to the nearest boundary segment. This method supports both
#' closed boundary polygons and open boundary edges (e.g., LINESTRINGs),
#' and is useful for modeling proximity-based enrichment near anatomical
#' or user-defined regions.
#'
#' Supports multiple decay methods for converting distance into weights.
#' Optionally scales distances to the \[0, 1\] range before computing weights.
#'
#' @importFrom methods is
#' @inheritParams computeCentroidWeights
#' @param boundary Either an `sf` object containing boundary
#'                 geometries (`POLYGON` or `LINESTRING`),
#'                 or a data frame of boundary points
#'                 returned from `getBoundary()`.
#'
#' @return A named numeric vector of weights, with names corresponding
#'         to `cell_ids`.
#'
#' @export
#' @examples
#' # Load spatial coordinates
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'     package = "SpNeigh"
#' ))
#'
#' # Get and build boundary polygons from cluster 2
#' boundary_points <- getBoundary(data = coords, one_cluster = 2)
#' boundary_polys <- buildBoundaryPoly(boundary_points)
#'
#' # Compute weights to polygon boundary
#' cells_c2 <- subset(coords, cluster == 2)$cell
#' weights <- computeBoundaryWeights(
#'     data = coords, cell_ids = cells_c2,
#'     boundary = boundary_polys[1, ]
#' )
#'
#' # Compute weights to a specific boundary edge
#' boundary_edges <- splitBoundaryPolyByAnchor(boundary_polys[1, ])
#' weights_edge <- computeBoundaryWeights(
#'     data = coords, cell_ids = cells_c2,
#'     boundary = boundary_edges[2, ]
#' )
#'
computeBoundaryWeights <- function(
    data = NULL,
    cell_ids = NULL,
    boundary = NULL,
    scale = TRUE,
    method = c("inverse", "gaussian", "linear", "quadratic"),
    sigma = 0.5) {
    method <- match.arg(method)

    # Extract coordinates
    if (is(data, "Seurat")) {
        sp_coords <- extractCoords(data, extract_cluster = FALSE)
    } else {
        sp_coords <- data
    }

    if (!all(c("cell", "x", "y") %in% colnames(sp_coords))) {
        stop("`data` must contain columns: 'cell', 'x', and 'y'.")
    }

    if (!inherits(boundary, c("sf", "data.frame"))) {
        stop(
            "`boundary` must be an sf object (POLYGON or LINESTRING)",
            " or a data frame of boundary points."
        )
    }

    # Convert boundary points to polygon if needed
    if (inherits(boundary, "sf")) {
        boundary_sf <- boundary
    } else {
        boundary_sf <- buildBoundaryPoly(boundary)
    }

    # Subset and convert coordinates to sf
    sub_coords <- sp_coords[sp_coords$cell %in% cell_ids, ]
    if (nrow(sub_coords) == 0) {
        stop("No matching cells found in 'data' for 'cell_ids'.")
    }

    cell_sf <- sf::st_as_sf(sub_coords, coords = c("x", "y"), crs = NA)

    # Handle boundary geometry type
    if (all(sf::st_geometry_type(boundary_sf) == "LINESTRING")) {
        boundary_lines <- boundary_sf
    } else {
        boundary_lines <- sf::st_boundary(boundary_sf)
    }

    # Compute distance matrix [cells x boundaries]
    dist_mat <- sf::st_distance(cell_sf, boundary_lines)
    min_dists <- apply(dist_mat, 1, min)
    dists <- as.numeric(min_dists)

    # Optional scaling
    if (scale && max(dists) > 0) {
        dists <- (dists - min(dists)) / (max(dists) - min(dists))
    }

    # Compute weights from distances
    weights <- switch(method,
        inverse   = 1 / (1 + dists),
        gaussian  = exp(-dists^2 / (2 * sigma^2)),
        linear    = 1 - dists,
        quadratic = (1 - dists)^2
    )

    # Clamp for scaled distances (to avoid negative weights)
    weights[dists > 1 & scale] <- 0
    weights[weights < 0] <- 0
    names(weights) <- sub_coords$cell

    return(weights)
}

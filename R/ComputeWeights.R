
#' Compute spatial weights based on distance to centroid
#'
#' Computes spatial weights for a set of cells based on their distance to
#' the centroid of their group. Supports multiple decay methods and optional
#' distance scaling.
#'
#' @param data A data frame or Seurat object containng at least three columns: `x`, `y` and `cell`.
#' @param cell_ids A vector of cell IDs to subset from data.
#' @param scale A logical value indicating whether to scale distances to \[0, 1\] before computing weights. Scaling is recommended, and default is TRUE.
#' @param method A string specifying the decay method. One of \code{"inverse"},
#'   \code{"gaussian"}, \code{"linear"}, or \code{"quadratic"}.
#' @param sigma Standard deviation for the Gaussian decay function. Default is 0.5 for scaled distance. If `scale = FALSE`, choose a value that is about half the region’s width.
#'
#' @return A named numeric vector of weights, one per cell in `cell_ids`.
#' @export
#' @examples
#' # Load coordinates
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'                               package = "SpNeigh"))
#' head(coords)
#'
#' # Compute spatial weights to centroid for cells in cluster 2
#' cells_c2 <- subset(coords, cluster==2)[,"cell"]
#' weights <- ComputeCentroidWeights(data = coords, cell_ids = cells_c2)
#' weights[1:3]
#'
ComputeCentroidWeights <- function(data = NULL,
                                   cell_ids = NULL,
                                   scale = TRUE,
                                   method = c("inverse", "gaussian", "linear", "quadratic"),
                                   sigma = 0.5) {
  method <- match.arg(method)


  # Extract coordinates from data
  if( inherits(data,"Seurat") ){
    sp_coords <- Seurat::GetTissueCoordinates(data)
  }else{
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
#' This function assigns each cell to the nearest boundary and computes
#' a spatial weight based on its distance to that boundary. Supports inverse,
#' Gaussian, linear, and quadratic decay methods.
#' @inheritParams ComputeCentroidWeights
#' @param boundary An sf object (POLYGON or LINESTRING) or a data.frame of boundary points.
#'
#' @return A named numeric vector of weights, one per cell in `cell_ids`.
#'
#' @export
#' @examples
#' # Load coordinates
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'                               package = "SpNeigh"))
#' head(coords)
#'
#' # Get boundaries for cells in cluster 2
#' boundary_points <- GetBoundary(data = coords, one_cluster = 2)
#' boundary_polys <- BuildBoundaryPoly(boundary_points)
#'
#' # Plot boundary regions
#' PlotRegion(boundary_polys)
#'
#' # Compute spatial weights to the boundary of region 1 for cells in cluster 2
#' cells_c2 <- subset(coords, cluster==2)[,"cell"]
#' weights <- ComputeBoundaryWeights(data = coords, cell_ids = cells_c2,
#'                                   boundary = boundary_polys[1,])
#' weights[1:3]
#'
#' # Compute spatial weights to the boundary edge2 of region 1 for cells in cluster 2
#' boundary_edges <- SplitBoundaryPolyByAnchor(boundary_polys[1,])
#' PlotEdge(boundary_edges)
#' weights <- ComputeBoundaryWeights(data = coords, cell_ids = cells_c2,
#'                                   boundary = boundary_edges[2,])
#' weights[1:3]
#' PlotWeights(coords, weights)
#'
ComputeBoundaryWeights <- function(data = NULL,
                                   cell_ids = NULL,
                                   boundary = NULL,
                                   scale = TRUE,
                                   method = c("inverse", "gaussian", "linear", "quadratic"),
                                   sigma = 0.5) {
  method <- match.arg(method)

  # Extract coordinates from data
  if( inherits(data,"Seurat") ){
    sp_coords <- Seurat::GetTissueCoordinates(data)
  }else{
    sp_coords <- data
  }

  if (!all(c("cell", "x", "y") %in% colnames(sp_coords))) {
    stop("data must contain columns: 'cell', 'x', and 'y'")
  }

  if (!inherits(boundary, c("sf", "data.frame"))) {
    stop("'boundary' must be an sf object (POLYGON or LINESTRING) or a data.frame containing boundary points.")
  }

  # Build safe polygons
  if(inherits(boundary, "sf")){
    boundary <- boundary
  }else{
    boundary <- BuildBoundaryPoly(boundary)
  }

  # Subset and convert to sf
  sub_coords <- sp_coords[sp_coords$cell %in% cell_ids, ]
  if (nrow(sub_coords) == 0) {
    stop("No matching cells found in 'data' for 'cell_ids'")
  }

  cell_sf <- sf::st_as_sf(sub_coords, coords = c("x", "y"), crs = NA)


  if (all(sf::st_geometry_type(boundary) == "LINESTRING")) {
    # Already LINESTRING
    boundary_lines <- boundary
  } else {
    # Extract polygon boundary as LINESTRING
    boundary_lines <- sf::st_boundary(boundary)
  }

  # Compute distances to all region boundaries → matrix: [cells × regions]
  dist_mat <- sf::st_distance(cell_sf, boundary_lines)
  min_dists <- apply(dist_mat, 1, min)
  dists <- as.numeric(min_dists)

  # Optional scaling
  if (scale && max(dists) > 0) {
    dists <- (dists - min(dists)) / (max(dists) - min(dists))
  }

  # Compute weights based on decay method
  # Compute weights based on selected method
  weights <- switch(method,
                    inverse = 1 / (1 + dists),
                    gaussian = exp(-dists^2 / (2 * sigma^2)),
                    linear = 1 - dists,
                    quadratic = (1 - dists)^2
  )

  # Ensure valid weights
  weights[dists > 1 & scale] <- 0  #clamp overshoot for linear/quadratic
  weights[weights < 0] <- 0 # Ensure non-negativity
  names(weights) <- sub_coords$cell

  return(weights)
}



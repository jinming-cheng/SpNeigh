
#' Removes spatial outlier cells based on local k-nearest neighbor distances
#'
#' Identifies and removes spatial outliers based on local density. For each cell,
#' the average distance to its k nearest neighbors is computed. Cells with a mean
#' k-NN distance greater than a specified cutoff are considered outliers and removed.
#' This helps retain densely connected cells while filtering isolated ones.
#'
#' @param coords A data frame or matrix of cell coordinates.
#' @param k The number of nearest neighbors to use for computing local distances. Default is 5.
#' @param distance_cutoff A numeric threshold on the average distance to neighbors.
#'        Cells with a higher mean distance will be removed. Default is 30.
#' @export
#' @examples
#' # Set seed for reproducibility
#' set.seed(123)
#'
#' # Generate 20 clustered points
#' n1 <- 20
#' x1 <- rnorm(n1, mean = 5, sd = 0.5)
#' y1 <- rnorm(n1, mean = 5, sd = 0.5)
#'
#' # Generate 5 outliers
#' n2 <- 5
#' x2 <- runif(n2, min = 10, max = 15)
#' y2 <- runif(n2, min = 10, max = 15)
#'
#' # Combine clustered points and outliers
#' coords <- data.frame(x = c(x1, x2), y = c(y1, y2) )
#' dim(coords)
#'
#' # Plot clustered points and outlier points
#' plot(coords$x, coords$y, pch = 16,
#'      col = rep(c("red","green"), times = c(20,5)) )
#'
#' # Remove outlier points with a specified knn distance cutoff
#' new_coords <- RemoveOutliers(coords, k = 5, distance_cutoff = 2)
#' dim(new_coords)
#'
#' # Returns TRUE meaning outliers have been removed in new_coords
#' all((new_coords$x == x1) & (new_coords$y == y1))
#'


RemoveOutliers <- function(coords, k = 5, distance_cutoff = 30) {

  if (!is.data.frame(coords) && !is.matrix(coords)) {
    stop("'coords' must be a data frame or matrix.")
  }
  if ( !(ncol(coords) == 2) ) {
    stop("'coords' must have two columns representing x and y coordinates.")
  }

  knn_res <- FNN::get.knn(coords, k = k)
  mean_dist <- rowMeans(knn_res$nn.dist)
  coords$mean_knn_dist <- mean_dist

  # Keep points with small average neighbor distance
  dplyr::filter(coords, .data$mean_knn_dist < distance_cutoff)
}


#' Extract spatial coordinates and cluster information
#'
#' Extracts a spatial coordinate data frame from either a Seurat object or a user-supplied data frame.
#' If a Seurat object is provided, the function retrieves the tissue coordinates along with cell cluster identities.
#' If a data frame is provided, it is returned as-is, after verifying that it contains the required columns.
#'
#' @importFrom methods is
#' @param data Either a Seurat object or a data frame with columns: `x`, `y`, `cell` and `cluster`.
#' @return A data frame with columns `x`, `y`, `cell`, and `cluster` representing spatial cell coordinates and cluster assignments.
#'
#' @export
#' @examples
#' # Load coordinates
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'                               package = "SpNeigh"))
#'
#' # An exmaple of the coordinate data.frame used for analysis
#' head(coords)
#'
#' # If the data is a data.frame, returns it without doing anything
#' head(ExtractCoords(coords))
#'

ExtractCoords <- function(data = NULL){
  # If data is a Seurat object, extract the coordinates and cluster information
  if( is(data,"Seurat") ){
    sp_coords <- Seurat::GetTissueCoordinates(data)
    sp_coords <- cbind.data.frame(sp_coords, cluster=data$seurat_clusters)
  }else{
    required_cols <- c("x", "y", "cell",  "cluster")
    if (!all(required_cols %in% colnames(data))) {
      stop("Input data frame must contain columns: x, y, cell, cluster")
    }

    # A data frame includes coordinates and cluster information for the spatial data
    sp_coords <- data
  }
  sp_coords
}


#' Extract spatial boundary points for a cluster or cell population
#'
#' Identifies and returns smoothed spatial boundary points for a specified cell cluster or population
#' using a concave hull algorithm. The function supports both single-region and multi-region boundaries.
#' In multi-region mode, subregions can be detected automatically via DBSCAN or manually using k-means clustering.
#' Outlier cells are first removed using a k-nearest neighbor distance filter to ensure boundary smoothness.
#'
#' @inheritParams RemoveOutliers
#' @inheritParams ExtractCoords
#' @importFrom rlang .data
#' @param one_cluster The cluster ID (numeric or character) to extract the boundary for.
#' @param multi_region Logical. If `TRUE`, identifies multiple spatial subregions within the cluster. Default is `TRUE`.
#' @param subregion_method Subregion detection method when `multi_region = TRUE`. Choose from `"dbscan"` (automatic) or `"kmeans"` (manual). Default is `"dbscan"`.
#' @param eps Neighborhood radius for DBSCAN subregion detection. Only used if `subregion_method = "dbscan"`. Default is 80.
#' @param minPts Minimum number of points for DBSCAN core point. Only used if `subregion_method = "dbscan"`. Default is 10.
#' @param n_subregions Number of subregions to use if `subregion_method = "kmeans"`. Default is 3.
#'
#' @return A data frame containing boundary points with columns `x`, `y`, and `region_id`.
#'         If `multi_region = TRUE`, multiple boundaries are returned and labeled by `region_id`.
#'         Otherwise, a single boundary is returned with `region_id = 1`.
#'
#' @export
#' @examples
#' # Load coordinates
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'                               package = "SpNeigh"))
#' head(coords)
#'
#' # Get boundary points of cluster 1 using dbscan method by default
#' boundary_points <- GetBoundary(data = coords, one_cluster = 1)
#'
#' # An example of data in the boundary_points data.frame
#' head(boundary_points)
#'
#' # Number of subregions
#' table(boundary_points$region_id)
#'
#' # Plot coordinates without boundaries
#' PlotBoundary(coords)
#'
#' # Plot boundary points of cluster 1
#' plot(boundary_points[,c("x","y")], pch=16, cex=0.3)
#'
#' # Plot boundary polygon of cluster 1
#' PlotBoundary(coords, one_cluster = 1)
#'
#' # Get boundary points of cluster 1 using kmeans method
#' # and manually specify subregion number
#' boundary_points <- GetBoundary(data = coords, one_cluster = 1,
#'                                subregion_method = "kmeans",
#'                                n_subregions = 2)
#' table(boundary_points$region_id)
#' plot(boundary_points[,c("x","y")], pch=16, cex=0.3)
#'
#' # Get boundary points of cluster 1 using concaveman hull method.
#' # region_id is set to 1 but it is not meaningful for this method.
#' boundary_points <- GetBoundary(data = coords, one_cluster = 1,
#'                                multi_region = FALSE)
#' table(boundary_points$region_id)
#' plot(boundary_points[,c("x","y")], pch=16, cex=0.3)

GetBoundary <- function(data = NULL,
                        one_cluster = NULL,
                        k = 5,
                        distance_cutoff = 30,
                        multi_region = TRUE,
                        subregion_method = c("dbscan","kmeans"),
                        eps = 80,
                        minPts = 10,
                        n_subregions = 3){

  subregion_method <- match.arg(subregion_method)

  # Extract coordinates from data
  sp_coords <- ExtractCoords(data = data)

  # Obtain coordinates of the target cluster
  coords_sub <- dplyr::filter(sp_coords,.data$cluster == one_cluster)

  # Remove outlier points
  clean_df <- RemoveOutliers(coords_sub[, c("x", "y")], k = k,
                              distance_cutoff = distance_cutoff)

  # Get smooth boundary using concave hull
  hull <- concaveman::concaveman(as.matrix(clean_df[, c("x", "y")]))
  hull <- as.data.frame(hull)
  hull$region_id <- 1

  # Rename V1 and V2 as x and y, respectively
  hull <- dplyr::rename(hull, x = "V1", y="V2")

  # If a cluster (or cell population) has multiple regions, obtain multiple boundaries
  if(multi_region){
    if (subregion_method == "dbscan") {
      # automatically detect the regions using dbscan method
      db <- dbscan::dbscan(clean_df[, c("x", "y")], eps = eps, minPts = minPts)
      clean_df$region_id <- db$cluster
      clean_df <- dplyr::filter(clean_df,.data$region_id != 0)  # Remove noise
    } else if (subregion_method == "kmeans") {
      # alternatively, manually specify the region number using kmeans method
      km <- stats::kmeans(clean_df[, c("x", "y")], centers = n_subregions)
      clean_df$region_id <- km$cluster
    }

    # A list of boundary points
    boundary_list <- lapply(unique(clean_df$region_id), function(reg_id) {
      region_df <- dplyr::filter(clean_df,.data$region_id == reg_id)
      hull_coords <- concaveman::concaveman(as.matrix(region_df[, c("x", "y")]))
      hull_df <- as.data.frame(hull_coords)
      hull_df$region_id <- reg_id
      return(hull_df)
    })

    all_boundaries <- dplyr::bind_rows(boundary_list)

    # Rename V1 and V2 as x and y, respectively
    all_boundaries <- dplyr::rename(all_boundaries, x = "V1", y="V2")
  }

  if(multi_region){
    return(all_boundaries)
  }else{
    return(hull)
  }

}


#' Convert boundary points into valid polygon geometries
#'
#' Constructs spatial polygon objects from boundary points generated by the `GetBoundary()` function.
#' Each set of boundary points (grouped by `region_id`) is converted into a closed polygon.
#' This function ensures that each polygon is valid and ready for downstream spatial analysis.
#'
#' @param boundary A data frame of boundary points, typically returned by `GetBoundary()`.
#'                 Must contain columns `x`, `y`, and `region_id`.
#'
#' @return An `sf` object containing one or more valid `POLYGON` geometries with a `region_id` column.
#'
#' @export
#'
#' @examples
#' # Load coordinates
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'                               package = "SpNeigh"))
#'
#' # Get boundary points for cluster 2
#' boundary_points <- GetBoundary(data = coords, one_cluster = 2)
#'
#' # Convert to polygon
#' boundary_polys <- BuildBoundaryPoly(boundary_points)
#'
#' # Plot the resulting polygons
#' plot(boundary_polys)
#'
BuildBoundaryPoly <- function(boundary = NULL) {
  if (is.null(boundary) || !all(c("x", "y", "region_id") %in% colnames(boundary))) {
    stop("Input must be a data frame with columns: x, y, and region_id.")
  }

  # Build polygon list per region
  polygon_list <- lapply(split(boundary, boundary$region_id), function(df) {
    coords <- as.matrix(df[, c("x", "y")])

    # Ensure the polygon is closed
    if (!all(coords[1, ] == coords[nrow(coords), ])) {
      coords <- rbind(coords, coords[1, ])
    }

    poly <- sf::st_polygon(list(coords))
    sf::st_make_valid(poly)
  })

  # Combine into sf object
  boundary_polys <- sf::st_sf(
    region_id = names(polygon_list),
    geometry = sf::st_sfc(polygon_list),
    crs = NA  # Neutral CRS to avoid longlat warnings
  )

  return(boundary_polys)
}


#' Convert boundary polygons to boundary point coordinates
#'
#' Extracts the vertex coordinates (`x`, `y`) of boundary geometries from an `sf` object
#' containing `POLYGON` or `LINESTRING` features. This is useful for recovering the original
#' boundary points from smoothed or labeled polygonal regions.
#'
#' @param boundary_poly An `sf` object containing only `POLYGON` or `LINESTRING` geometries.
#'                      Must include a `region_id` column for labeling subregions.
#'
#' @return A data frame with columns `x`, `y`, and `region_id`, containing the vertex coordinates
#'         of the boundary geometries.
#'
#' @export
#'
#' @examples
#' # Load coordinates and generate boundary
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'                               package = "SpNeigh"))
#' boundary_points <- GetBoundary(data = coords, one_cluster = 2)
#' boundary_polys <- BuildBoundaryPoly(boundary_points)
#'
#' # Convert back to boundary points
#' boundary_pts <- BoundaryPolyToPoints(boundary_polys)
#' head(boundary_pts)
#'
BoundaryPolyToPoints <- function(boundary_poly) {
  if (!inherits(boundary_poly, "sf")) {
    stop("Input must be an 'sf' object.")
  }

  # Ensure geometry type is either POLYGON or LINESTRING
  valid_types <- c("POLYGON", "LINESTRING")
  geom_type <- unique(sf::st_geometry_type(boundary_poly))

  if (!all(geom_type %in% valid_types)) {
    stop("All geometries in the input must be either POLYGON or LINESTRING.")
  }

  # Convert to POINT geometries
  boundary_points <- suppressWarnings(
    sf::st_cast(boundary_poly, "POINT", do_split = TRUE)
  )

  # Extract point coordinates
  coords <- sf::st_coordinates(boundary_points)

  # Return data frame with region_id and (x, y)
  df <- boundary_points |>
    sf::st_drop_geometry() |>
    dplyr::mutate(x = coords[, 1], y = coords[, 2]) |>
    dplyr::select(.data$x, .data$y, .data$region_id)

  return(df)
}


#' Generate an outer or inner boundary polygon by buffering an existing boundary
#'
#' Computes an expanded or shrunken boundary by applying a spatial buffer to an existing polygon.
#' The input can be either boundary points (as returned by `GetBoundary()`) or polygon geometries
#' (as returned by `BuildBoundaryPoly()`). This is useful for defining outer spatial neighborhoods
#' or for shrinking boundaries inward to define inner regions.
#'
#' Be careful: when using a negative buffer distance (for inner boundaries),
#' polygons may collapse, become invalid, or disappear entirely if the buffer width exceeds the shape's interior size.
#'
#'
#' @param boundary A data frame of boundary points (with columns `x`, `y`, `region_id`) or an `sf` object.
#' @param dist A numeric value specifying the buffer distance in spatial units.
#'        Use positive values to expand (outer boundary), and negative values to shrink (inner boundary). Default is 100.
#'
#' @return An `sf` object of expanded outer boundary polygons with the same `region_id` values as the original input.
#'
#' @seealso [GetInnerBoundary()] for a simplified wrapper for inward shrinking.
#'
#' @export
#'
#' @examples
#' # Load coordinates
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'                               package = "SpNeigh"))
#'
#' # Get boundary points of cluster 2
#' boundary_points <- GetBoundary(data = coords, one_cluster = 2)
#'
#' # Build polygons from boundary points
#' boundary_polys <- BuildBoundaryPoly(boundary_points)
#'
#' # Generate outer boundaries with 100-unit buffer
#' outer1 <- GetOuterBoundary(boundary_points, dist = 100)
#' outer2 <- GetOuterBoundary(boundary_polys, dist = 100)
#'
#' # Plot original and expanded boundaries
#' plot(boundary_polys)
#' plot(outer2, add = TRUE, border = "red")
#'
GetOuterBoundary <- function(boundary = NULL, dist = 100) {
  # Convert to sf POLYGON if needed
  if (inherits(boundary, "sf")) {
    boundary_polys <- boundary
  } else {
    boundary_polys <- BuildBoundaryPoly(boundary)
  }

  # Expand boundary using a buffer
  outer_boundaries <- dplyr::mutate(
    boundary_polys,
    geometry = sf::st_buffer(.data$geometry, dist = dist)
  )

  return(outer_boundaries)
}


#' Generate an inner boundary polygon by shrinking an existing boundary inward
#'
#' Computes an inward-shifted (contracted) boundary by applying a negative spatial buffer
#' to an existing boundary polygon. This function is a wrapper around `GetOuterBoundary()`
#' and is useful for defining an inner region around a spatial cluster or structure.
#'
#' Be careful: if the shrinkage distance is too large, the resulting geometry may become invalid
#' or disappear entirely, especially for narrow or irregular shapes.
#'
#' @param boundary A data frame of boundary points (with columns `x`, `y`, `region_id`)
#'        or an `sf` object of `POLYGON` geometries.
#' @param dist A positive numeric value specifying how far inward to shrink the boundary.
#'        This is automatically converted to a negative buffer distance. Default is 50.
#'
#' @return An `sf` object containing the inward-shrunk polygons, one per `region_id`.
#'
#' @seealso [GetOuterBoundary()] for the outward (positive) buffer version.
#'
#' @export
#'
#' @examples
#' # Load coordinates
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'                               package = "SpNeigh"))
#'
#' # Get boundary polygons of cluster 2
#' boundary_points <- GetBoundary(data = coords, one_cluster = 2)
#' boundary_polys <- BuildBoundaryPoly(boundary_points)
#' plot(boundary_polys)
#'
#' # Generate inner boundary with 50-unit buffer for boundary region 1
#' inner_boundary <- GetInnerBoundary(boundary_polys)
#' plot(inner_boundary)
#'

GetInnerBoundary <- function(boundary = NULL, dist = 50) {
  if (!is.numeric(dist) || dist <= 0) {
    stop("'dist' must be a positive number.")
  }

  inner <- GetOuterBoundary(boundary, dist = -abs(dist))

  # Optional: warn about empty geometries (fully collapsed)
  if (any(sf::st_is_empty(inner))) {
    warning("Some inner boundaries may have collapsed or disappeared due to excessive shrinkage.")
  }

  return(inner)
}



#' Generate ring regions between a boundary and its outer buffer
#'
#' Computes spatial ring-shaped regions by subtracting the original boundary polygons
#' from their corresponding outer buffered polygons. If the `outer_boundary` is not supplied,
#' it will be automatically computed using `GetOuterBoundary()`. This is useful for
#' analyzing periphery-enriched cell types or gradient-based features near a boundary.
#'
#' @inheritParams GetOuterBoundary
#' @param outer_boundary Optional `sf` object containing buffered (outer) boundary polygons.
#'        If not provided, it will be automatically computed using `GetOuterBoundary()`.
#' @param ...  Additional arguments passed to `GetOuterBoundary()` if `outer_boundary` is not provided.
#' @export
#' @examples
#' # Load coordinates
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'                               package = "SpNeigh"))
#' head(coords)
#'
#' # Get boundary points of cluster 2
#' boundary_points <- GetBoundary(data = coords, one_cluster = 2)
#'
#' # Automatically compute outer boundary and get rings
#' ring_regions <- GetRingRegion(boundary = boundary_points, dist = 100)
#' plot(ring_regions)
#'
#' # Or provide both inner and outer boundaries explicitly
#' outer <- GetOuterBoundary(boundary_points, dist = 100)
#' rings <- GetRingRegion(boundary = boundary_points, outer_boundary = outer)
#' plot(rings)
#'
GetRingRegion <- function(boundary = NULL, outer_boundary = NULL,...){
  # Build safe polygons
  if(inherits(boundary,"sf")){
    boundary_polys <- boundary
  }else{
    boundary_polys <- BuildBoundaryPoly(boundary)
  }

  # If outer boundary is not provided, get the outer boundary
  if(is.null(outer_boundary)){
    outer_boundaries <- GetOuterBoundary(boundary = boundary_polys,...)
  }else{
    outer_boundaries <- outer_boundary
  }

  # Subtract each original region from its own buffer (row-wise difference)
  # Use purrr::map2 to iterate
  ring_geoms <- purrr::map2(outer_boundaries$geometry,
                            boundary_polys$geometry, ~ sf::st_difference(.x, .y))

  # Create ring_regions as a new sf object
  ring_regions <- sf::st_sf(
    region_id = boundary_polys$region_id,
    geometry = sf::st_sfc(ring_geoms),
    crs = NA)

  ring_regions

}


#' Split a polygon boundary into two parts using anchor points
#'
#' This function takes an sf POLYGON object and two anchor points (or infers them)
#' to split the polygon boundary into two LINESTRING edge segments.
#'
#' @param boundary_poly An sf POLYGON object.
#' @param pt1 A numeric vector of length 2 (x, y) or NULL. If NULL, the leftmost boundary point is used.
#' @param pt2 A numeric vector of length 2 (x, y) or NULL. If NULL, the rightmost boundary point is used.
#'
#' @return An sf object with two LINESTRING features and a 'region_id' column indicating edge1 and edge2.
#'
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
#' # Split boundary polygon 1 into two edges using leftmost and rightmost anchor points
#' boundary_edges <- SplitBoundaryPolyByAnchor(boundary_polys[1,])
#' plot(boundary_edges, lwd = 2)
#'
#' # Split boundary polygon 1 into two edges using two anchor points
#' pt1 <- c(4000, 1500)
#' pt2 <- c(2000, 3000)
#' boundary_edges <- SplitBoundaryPolyByAnchor(boundary_polys[1,],
#'                                             pt1 = pt1,
#'                                             pt2 = pt2)
#' plot(boundary_edges, lwd = 2)
#'
#' @export
SplitBoundaryPolyByAnchor <- function(boundary_poly = NULL,
                                      pt1 = NULL,
                                      pt2 = NULL) {
  if (!inherits(boundary_poly, "sf") || !all(sf::st_geometry_type(boundary_poly) == "POLYGON")) {
    stop("Input must be an sf object with POLYGON geometry.")
  }

  coords <- sf::st_coordinates(sf::st_geometry(boundary_poly)[[1]])[, 1:2]

  # Auto-infer pt1 and pt2 if NULL
  if (is.null(pt1)) {
    pt1 <- coords[which.min(coords[, 1]), ]
  }
  if (is.null(pt2)) {
    pt2 <- coords[which.max(coords[, 1]), ]
  }

  if (!is.numeric(pt1) || length(pt1) != 2 || !is.numeric(pt2) || length(pt2) != 2) {
    stop("pt1 and pt2 must each be a numeric vector of length 2 (x, y) or NULL.")
  }

  # Find nearest boundary points to pt1 and pt2
  d1 <- sqrt(rowSums((coords - matrix(pt1, nrow(coords), 2, byrow = TRUE))^2))
  d2 <- sqrt(rowSums((coords - matrix(pt2, nrow(coords), 2, byrow = TRUE))^2))
  idx1 <- which.min(d1)
  idx2 <- which.min(d2)

  n <- nrow(coords)
  path1 <- coords[c(idx1:idx2), ]
  path2 <- coords[c(idx2:n, 1:idx1), ]

  # Construct LINESTRINGs
  line1 <- sf::st_linestring(path1)
  line2 <- sf::st_linestring(path2)

  lines <- sf::st_sfc(line1, line2, crs = sf::st_crs(boundary_poly))
  sf::st_sf(region_id = c("edge1", "edge2"), geometry = lines)
}

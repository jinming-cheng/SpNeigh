
#' Removes outlier cells that are too far from their neighbors in a cluster or cell population
#'
#' @param coords coordinates of the cells
#' @param k The maximum number of nearest neighbors to search, which is a parameter used in `get.knn` function. Default is set to 5.
#' @param distance_cutoff K-nearest neighbor distance threshold. Default is set to 30. k and distance_cutoff are used to remove outlier points.
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



RemoveOutliers <- function(coords, k = 5, distance_cutoff = 30) {
  knn_res <- FNN::get.knn(coords, k = k)
  mean_dist <- rowMeans(knn_res$nn.dist)
  coords$mean_knn_dist <- mean_dist

  # Keep points with small average neighbor distance
  dplyr::filter(coords, .data$mean_knn_dist < distance_cutoff)
}


#' Extract coordinates from a Seurat object or a data frame
#'
#' @param data A data frame contains coordinates of the cells and cluster information. If the input is a Seurat object, the information will be extracted automatically.
#' @export
#' @examples
#' # Load coordinates
#' load(system.file("extdata", "MouseBrainTinyCoords.rda",
#'                  package = "SpNeigh"))
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
    # A data frame includes coordinates and cluster information for the spatial data
    sp_coords <- data
  }
  sp_coords
}

#' Get boundary points for a cluster or cell population
#'
#' @inheritParams RemoveOutliers
#' @inheritParams ExtractCoords
#' @importFrom methods is
#' @importFrom rlang .data
#' @param one_cluster Specify one cluster to draw the boundary
#' @param multi_region A logical value demonstrates whether a cluster has multiple regions on a spatial plot. Default is TRUE
#' @param subregion_method If `multi_region` is TRUE, set a method to detect the subregions. Default is "dbscan". The "dbscan" method can be used to automatically detect subregions. Alternatively, "kmeans" method can be used to manually set the number of subregions.
#' @param eps A parameter used when `subregion_method` is set to "dbscan". Size of the epsilon neighborhood.
#' @param minPts A parameter used when `subregion_method` is set to "dbscan". Number of minimum points required in the eps neighborhood for core points. The `eps` and `minPts` parameters are used to adjust the subregion numbers for "dbscan" method.
#' @param n_subregions A parameter used when `subregion_method` is set to "kmeans". Manually set the number of subregions. Default is 3.
#' @export
#' @examples
#' # Load coordinates
#' load(system.file("extdata", "MouseBrainTinyCoords.rda",
#'                  package = "SpNeigh"))
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
                        subregion_method = c("dbscan","kmeans")[1],
                        eps = 80,
                        minPts = 10,
                        n_subregions = 3){

  # extract coordinates from data
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

  # rename V1 and V2 as x and y, respectively
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

    # rename V1 and V2 as x and y, respectively
    all_boundaries <- dplyr::rename(all_boundaries, x = "V1", y="V2")
  }

  if(multi_region){
    return(all_boundaries)
  }else{
    return(hull)
  }

}


#' Build boundary polygon based on the boundary points
#'
#' @param boundary Boundary result returned by `GetBoundary` function.
#' @export
#' @examples
#' # Load coordinates
#' load(system.file("extdata", "MouseBrainTinyCoords.rda",
#'                  package = "SpNeigh"))
#' head(coords)
#'
#' # Get boundary points of cluster 2
#' boundary_points <- GetBoundary(data = coords, one_cluster = 2)
#'
#' # Build boundary polygon from the boundary points
#' boundary_polys = BuildBoundaryPoly(boundary_points)
#'
#' # Plot boundary polygons
#' plot(boundary_polys)
#'
BuildBoundaryPoly <- function(boundary = NULL){
  # Build safe polygons
  polygon_list <- lapply(split(boundary,
                               boundary$region_id),
                         function(df) {
                           coords <- as.matrix(df[, c("x", "y")])

                           # Ensure it's closed
                           if (!all(coords[1, ] == coords[nrow(coords), ])) {
                             coords <- rbind(coords, coords[1, ])}

                           # Create polygon and make it valid
                           poly <- sf::st_polygon(list(coords))
                           sf::st_make_valid(poly)}
  )

  # Combine into sf object — use NA CRS or 3857 to avoid longlat issues
  boundary_polys <- sf::st_sf(
    region_id = names(polygon_list),
    geometry = sf::st_sfc(polygon_list),
    crs = NA)

  boundary_polys
}


#' Get outer boundary of the original boundary
#'
#' @param boundary A data frame or sf object. Boundary result returned by `GetBoundary` or `GetBoundaryPoly` function.
#' @param dist Distance of the outer boundary from the existing boundary.
#' @export
#' @examples
#' # Load coordinates
#' load(system.file("extdata", "MouseBrainTinyCoords.rda",
#'                  package = "SpNeigh"))
#' head(coords)
#'
#' # Get boundary points of cluster 2
#' boundary_points <- GetBoundary(data = coords, one_cluster = 2)
#'
#' # Build boundary polygon from the boundary points
#' boundary_polys = BuildBoundaryPoly(boundary_points)
#'
#' # Get outer boundary ploygons that are of a distance of 100 units from the original boundaries
#' # both boundary points and boundary polygons can be used as input
#' outer_boundary_polys <- GetOuterBoundary(boundary_points, dist = 100)
#' outer_boundary_polys <- GetOuterBoundary(boundary_polys, dist = 100)
#'
#' # Plot original boundary polygons
#' plot(boundary_polys)
#'
#' # Plot outer boundary polygons
#' plot(outer_boundary_polys)
#'
GetOuterBoundary <- function(boundary = NULL, dist = 100){
  # Build a safe polygon
  if(is(boundary,"sf")){
    boundary_polys <- boundary
  }else{
    boundary_polys <- BuildBoundaryPoly(boundary)
  }

  outer_boundaries <- dplyr::mutate(boundary_polys,
                                    geometry = sf::st_buffer(.data$geometry, dist = dist))

  outer_boundaries

}


#' Get ring region between the original boundary and outer boundary
#'
#' @inheritParams GetOuterBoundary
#' @param outer_boundary A sf object. Outer boundary returned by `GetOuterBoundary` function.
#' @param ... Parameters in `GetOuterBoundary`
#' @export
#' @examples
#' # Load coordinates
#' load(system.file("extdata", "MouseBrainTinyCoords.rda",
#'                  package = "SpNeigh"))
#' head(coords)
#'
#' # Get boundary points of cluster 2
#' boundary_points <- GetBoundary(data = coords, one_cluster = 2)
#'
#' # Get ring regions (outer regions with be automatically obtained)
#' ring_regions <- GetRingRegion(boundary = boundary_points, dist = 100)
#'
#' # Plot ring regions
#' plot(ring_regions)
#'
#' # Alternatively, get ring region using the provided boundary and outer_boundary
#' outer_boundary_polys <- GetOuterBoundary(boundary_points, dist = 100)
#' ring_regions <- GetRingRegion(boundary = boundary_points,
#'                               outer_boundary = outer_boundary_polys)
#' plot(ring_regions)
#'
GetRingRegion <- function(boundary = NULL, outer_boundary = NULL,...){
  # Build safe polygons
  if(is(boundary,"sf")){
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


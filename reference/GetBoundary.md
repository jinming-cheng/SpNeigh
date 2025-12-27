# Extract spatial boundary points for a cluster or cell population

Identifies and returns smoothed spatial boundary points for a specified
cell cluster or population using a concave hull algorithm. The function
supports both single-region and multi-region boundaries. In multi-region
mode, subregions can be detected automatically via DBSCAN or manually
using k-means clustering. Outlier cells are first removed using a
k-nearest neighbor distance filter to ensure boundary smoothness.

## Usage

``` r
GetBoundary(
  data = NULL,
  cluster_col = NULL,
  one_cluster = NULL,
  k = 5,
  distance_cutoff = 30,
  multi_region = TRUE,
  subregion_method = c("dbscan", "kmeans"),
  eps = 80,
  minPts = 10,
  n_subregions = 3
)
```

## Arguments

- data:

  A `Seurat` object, a `SpatialExperiment` object, or a data frame
  containing spatial coordinates.

- cluster_col:

  Character scalar specifying the metadata column name containing
  cluster assignments. If `NULL`, a default is used depending on the
  input object type:

  - `"seurat_clusters"` for `Seurat` objects

  - `"cluster"` for `SpatialExperiment` objects

- one_cluster:

  The cluster ID (numeric or character) to extract the boundary for.

- k:

  The number of nearest neighbors to use for computing local distances.
  Default is 5.

- distance_cutoff:

  A numeric threshold on the average distance to neighbors. Cells with a
  higher mean distance will be removed. Default is 30.

- multi_region:

  Logical. If `TRUE`, identifies multiple spatial subregions within the
  cluster. Default is `TRUE`.

- subregion_method:

  Subregion detection method when `multi_region = TRUE`. Choose from
  "dbscan" (automatic) or "kmeans" (manual). Default is "dbscan".

- eps:

  Neighborhood radius for `DBSCAN` subregion detection. Only used if
  `subregion_method = "dbscan"`. Default is 80.

- minPts:

  Minimum number of points for `DBSCAN` core point. Only used if
  `subregion_method = "dbscan"`. Default is 10.

- n_subregions:

  Number of subregions to use if `subregion_method = "kmeans"`. Default
  is 3.

## Value

A data frame containing boundary points with columns `x`, `y`, and
`region_id`. If `multi_region = TRUE`, multiple boundaries are returned
and labeled by `region_id`. Otherwise, a single boundary is returned
with `region_id = 1`.

## Examples

``` r
# Load coordinates
coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
    package = "SpNeigh"
))
head(coords)
#>          x        y cell cluster
#> 1 1898.815 2540.963    1       4
#> 2 1895.305 2532.627    2       4
#> 3 2368.073 2534.409    3       2
#> 4 1903.726 2560.010    4       4
#> 5 1917.481 2543.132    5       4
#> 6 1926.540 2560.044    6       4

# Get boundary points of cluster 1
boundary_points <- GetBoundary(data = coords, one_cluster = 1)
head(boundary_points)
#>          x        y region_id
#> 1 447.7028 3404.517         1
#> 2 456.0788 3429.117         1
#> 3 475.1725 3440.692         1
#> 4 475.9918 3453.765         1
#> 5 466.2819 3466.096         1
#> 6 454.2985 3474.720         1

# Number of subregions
table(boundary_points$region_id)
#> 
#>   1   2   3   4   5 
#> 657 146  17  40   7 

# Get boundary points of cluster 1 using kmeans method
# and manually specify subregion number
boundary_points <- GetBoundary(
    data = coords,
    one_cluster = 1,
    subregion_method = "kmeans",
    n_subregions = 2
)
table(boundary_points$region_id)
#> 
#>   1   2 
#> 353 538 

# Get boundary points of cluster 1 without multiple regions
boundary_points <- GetBoundary(
    data = coords,
    one_cluster = 1,
    multi_region = FALSE
)
table(boundary_points$region_id)
#> 
#>   1 
#> 853 
```

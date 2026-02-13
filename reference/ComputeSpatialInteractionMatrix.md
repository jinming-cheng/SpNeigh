# Compute a spatial neighborhood interaction matrix using K-nearest neighbors (KNN)

Computes a spatial interaction matrix where each entry quantifies the
number of neighboring cells from a given cluster (columns) that are
among the `k`-nearest neighbors of cells in another cluster (rows). This
provides a summary of spatial proximity and enrichment of neighboring
clusters around each focal cluster.

## Usage

``` r
computeSpatialInteractionMatrix(data = NULL, cluster_col = NULL, k = 10)
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

- k:

  Integer. Number of nearest neighbors to use for each cell. Default is
  10.

## Value

A numeric matrix where rows represent focal clusters and columns
represent neighboring clusters. Each cell in the matrix indicates how
frequently a neighbor cluster appears among the k-nearest neighbors of
cells from the focal cluster.

## Details

The matrix is built by identifying the `k` nearest neighbors for each
cell based on spatial coordinates, and then tabulating the cluster
identities of those neighbors with respect to the cluster identity of
the focal cell.

## Examples

``` r
# Load coordinates
coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
    package = "SpNeigh"
))

# Compute interaction matrix using all cells
interaction_matrix <- computeSpatialInteractionMatrix(coords)
head(interaction_matrix)
#>       0     1     2     3     4     5    6   7   8   9   10  11
#> 0 26326  6590  5670  7847  4608  4212 2814 548 364 542  456 193
#> 1  6540 33843  2707  7455   607    15 5119   0   0   0   11 203
#> 2  5956  2774 30886  4694  5792  1908 1914 219  79 230 1393 155
#> 3  8572  7701  4895 11808  5558   901 3666 881 628 980  146 194
#> 4  4729   624  5754  5605 23932     0 3589   0   0   0    0   7
#> 5  4337    12  1682   849     0 29000   27   6   0   1    6   0

# Compute interaction matrix for cells inside boundaries
boundary_points <- getBoundary(
    data = coords, one_cluster = 2,
    eps = 120, minPts = 10
)
cells_inside <- getCellsInside(data = coords, boundary = boundary_points)
coords_sub <- subset(coords, cell %in% cells_inside$cell)
computeSpatialInteractionMatrix(coords_sub)
#>       0     2   3   4   5   6   10
#> 0  5115  1965  94  65  24  17  410
#> 2  2199 27625 923 558 187 201 1337
#> 3   100   893  50 136  22  16  133
#> 4    70   526 137 581   0  26    0
#> 5    25   129  18   0 942   0    6
#> 6    20   211  16  25   0  18   50
#> 10  416  1274 123   0   6  54 3987

# Compute interaction matrix for cells inside ring region 2
ring_regions <- getRingRegion(boundary = boundary_points, dist = 100)
cells_ring <- getCellsInside(data = coords, boundary = ring_regions[2, ])
coords_sub <- subset(coords, cell %in% cells_ring$cell)
computeSpatialInteractionMatrix(coords_sub)
#>      0   2   3    5 6 7
#> 0 2634 219 230  476 8 3
#> 2  253 183  51  527 6 0
#> 3  304  52 230   66 7 1
#> 5  492 384  54 7479 1 0
#> 6    6   5   7    2 0 0
#> 7    4   0   6    0 0 0
```

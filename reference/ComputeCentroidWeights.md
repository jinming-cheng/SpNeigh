# Compute spatial weights based on distance to the centroid

Computes spatial weights for a set of cells based on their Euclidean
distance to the centroid of the selected group of cells. This is useful
for modeling gradient-like spatial expression centered around a cluster
or region. Supports multiple decay methods and optional distance
scaling.

## Usage

``` r
computeCentroidWeights(
  data = NULL,
  cell_ids = NULL,
  scale = TRUE,
  method = c("inverse", "gaussian", "linear", "quadratic"),
  sigma = 0.5
)
```

## Arguments

- data:

  A data frame, Seurat object or SpatialExperiment object containing
  spatial coordinates. Must include columns: `cell`, `x`, and `y`.

- cell_ids:

  A character vector of cell IDs to use for weight computation.

- scale:

  Logical. Whether to scale distances to the range \[0, 1\] before
  applying decay functions. Default is `TRUE`.

- method:

  Decay function to convert distances to weights. One of: "inverse",
  "gaussian", "linear", or "quadratic".

- sigma:

  Standard deviation for the Gaussian decay (used only if
  `method = "gaussian"`). Default is `0.5` (recommended for scaled
  distances).

## Value

A named numeric vector of weights (length = number of cells in
`cell_ids`). Higher weights correspond to cells closer to the centroid.

## Examples

``` r
# Load spatial coordinates
coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
    package = "SpNeigh"
))

# Select cells from cluster 2
cells_c2 <- subset(coords, cluster == 2)$cell

# Compute centroid-based weights using default settings
weights <- computeCentroidWeights(data = coords, cell_ids = cells_c2)
head(weights)
#>         3         8         9        11        18        20 
#> 0.7829382 0.7033367 0.7067281 0.7090737 0.7185800 0.8012193 
```

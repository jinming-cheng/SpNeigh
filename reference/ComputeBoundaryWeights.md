# Compute spatial weights based on distance to nearest boundary

Computes spatial weights for a subset of cells based on their Euclidean
distance to the nearest boundary segment. This method supports both
closed boundary polygons and open boundary edges (e.g., LINESTRINGs),
and is useful for modeling proximity-based enrichment near anatomical or
user-defined regions.

## Usage

``` r
computeBoundaryWeights(
  data = NULL,
  cell_ids = NULL,
  boundary = NULL,
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

- boundary:

  Either an `sf` object containing boundary geometries (`POLYGON` or
  `LINESTRING`), or a data frame of boundary points returned from
  [`getBoundary()`](https://github.com/jinming-cheng/SpNeigh/reference/GetBoundary.md).

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

A named numeric vector of weights, with names corresponding to
`cell_ids`.

## Details

Supports multiple decay methods for converting distance into weights.
Optionally scales distances to the \[0, 1\] range before computing
weights.

## Examples

``` r
# Load spatial coordinates
coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
    package = "SpNeigh"
))

# Get and build boundary polygons from cluster 2
boundary_points <- getBoundary(data = coords, one_cluster = 2)
boundary_polys <- buildBoundaryPoly(boundary_points)

# Compute weights to polygon boundary
cells_c2 <- subset(coords, cluster == 2)$cell
weights <- computeBoundaryWeights(
    data = coords, cell_ids = cells_c2,
    boundary = boundary_polys[1, ]
)

# Compute weights to a specific boundary edge
boundary_edges <- splitBoundaryPolyByAnchor(boundary_polys[1, ])
weights_edge <- computeBoundaryWeights(
    data = coords, cell_ids = cells_c2,
    boundary = boundary_edges[2, ]
)
```

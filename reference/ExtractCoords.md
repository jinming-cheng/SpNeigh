# Extract spatial coordinates and cluster information

This function extracts a spatial coordinate data frame from either a
Seurat object or a user-supplied data frame. When a Seurat object is
provided, it retrieves tissue coordinates (`x`, `y`) and optionally the
cluster identities (e.g., `seurat_clusters`). When a data frame is
provided, it is returned as-is after checking that it contains the
required columns.

## Usage

``` r
ExtractCoords(data = NULL, extract_cluster = TRUE)
```

## Arguments

- data:

  A Seurat object or a data frame with columns: `x`, `y`, `cell`, and
  `cluster`. If a Seurat object is provided, the `seurat_clusters`
  metadata column will be used as the `cluster`.

- extract_cluster:

  Logical. Whether to extract the `seurat_clusters` column from the
  Seurat object when available. Ignored if `data` is a data frame.
  Default is `TRUE`.

## Value

A data frame with columns `x`, `y`, `cell`, and optionally `cluster`,
representing spatial coordinates and cluster assignments.

## Examples

``` r
# Load coordinates
coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
    package = "SpNeigh"
))

# An exmaple of the coordinate data.frame used for analysis
head(coords)
#>          x        y cell cluster
#> 1 1898.815 2540.963    1       4
#> 2 1895.305 2532.627    2       4
#> 3 2368.073 2534.409    3       2
#> 4 1903.726 2560.010    4       4
#> 5 1917.481 2543.132    5       4
#> 6 1926.540 2560.044    6       4

# If the data is a data.frame, returns it without doing anything
head(ExtractCoords(coords))
#>          x        y cell cluster
#> 1 1898.815 2540.963    1       4
#> 2 1895.305 2532.627    2       4
#> 3 2368.073 2534.409    3       2
#> 4 1903.726 2560.010    4       4
#> 5 1917.481 2543.132    5       4
#> 6 1926.540 2560.044    6       4
```

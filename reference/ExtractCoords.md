# Extract spatial coordinates and cluster information

Extract spatial coordinates and (optionally) cluster assignments from a
`Seurat` object, a `SpatialExperiment` object, or a user-supplied data
frame.

## Usage

``` r
extractCoords(data, cluster_col = NULL, extract_cluster = TRUE, ...)
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

- extract_cluster:

  Logical indicating whether to extract cluster information. If `FALSE`,
  only spatial coordinates and cell IDs are returned. Default is `TRUE`.

- ...:

  Additional arguments (unused).

## Value

A data frame with columns `x`, `y`, and `cell`, and optionally
`cluster`.

## Details

This is an S3 generic with methods implemented for `Seurat`,
`SpatialExperiment`, and `data.frame`.

## Examples

``` r
coords <- readRDS(system.file(
    "extdata", "MouseBrainCoords.rds",
    package = "SpNeigh"
))

head(extractCoords(coords))
#>   cell        x        y cluster
#> 1    1 1898.815 2540.963       4
#> 2    2 1895.305 2532.627       4
#> 3    3 2368.073 2534.409       2
#> 4    4 1903.726 2560.010       4
#> 5    5 1917.481 2543.132       4
#> 6    6 1926.540 2560.044       4
```

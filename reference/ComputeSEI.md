# Compute Spatial Enrichment Index (SEI)

Alias for
[`computeSpatialEnrichmentIndex`](https://github.com/jinming-cheng/SpNeigh/reference/ComputeSpatialEnrichmentIndex.md).

## Usage

``` r
computeSEI(exp_mat = NULL, weights = NULL)
```

## Arguments

- exp_mat:

  A normalized gene expression matrix with genes as rows and cells as
  columns. Should be of class `matrix` or `dgCMatrix`.

- weights:

  A numeric vector of spatial weights (e.g., from
  `computeBoundaryWeights` or `computeCentroidWeights`). Must be the
  same length as the number of columns (cells) in `exp_mat`.

## Value

A data frame containing the spatial enrichment index (SEI) results.

## Examples

``` r
# Load spatial coordinates and log-normalized expression
coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
    package = "SpNeigh"
))
logNorm_expr <- readRDS(system.file("extdata", "LogNormExpr.rds",
    package = "SpNeigh"
))

# Compute spatial weights and SEI
bon_c0 <- getBoundary(data = coords, one_cluster = 0)
cells_c0 <- subset(coords, cluster == 0)$cell
weights <- computeBoundaryWeights(
    data = coords,
    cell_ids = cells_c0,
    boundary = bon_c0
)
sei_df <- computeSEI(logNorm_expr[, cells_c0], weights)
head(sei_df)
#>      gene      SEI mean_expr normalized_SEI
#> 1    Fmod 1.742086  1.562370       1.115027
#> 2    Gjb2 1.325063  1.193263       1.110452
#> 3 Aldh1a2 1.820161  1.642933       1.107872
#> 4 Slc13a4 2.106610  1.906511       1.104955
#> 5  Col1a1 1.842083  1.675953       1.099125
#> 6  Cyp1b1 1.495145  1.365958       1.094575
```

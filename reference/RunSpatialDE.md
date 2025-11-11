# Differential expression along spatial distance gradients using splines

Performs spatially-aware differential expression (DE) analysis by
modeling gene expression as a smooth function of a continuous spatial
covariate (e.g., distance to a boundary or centroid). Natural spline
basis functions are used to capture non-linear trends in expression
relative to spatial distance. This method is suitable for identifying
genes whose expression varies continuously across spatial structures.

## Usage

``` r
RunSpatialDE(
  exp_mat = NULL,
  cell_ids = NULL,
  spatial_distance = NULL,
  adj_p.value = 0.05,
  df = 3
)
```

## Arguments

- exp_mat:

  A normalized gene expression matrix (genes x cells), either a `matrix`
  or `dgCMatrix`. Typically log-normalized counts, e.g., from a Seurat
  object.

- cell_ids:

  A character vector of cell IDs (column names of `exp_mat`) used for DE
  analysis.

- spatial_distance:

  A named numeric vector containing the spatial distance (or weights)
  for each cell. Must be the same length as `cell_ids`. Scaled distances
  are recommended.

- adj_p.value:

  Adjusted p-value threshold for reporting differentially expressed
  genes. Default is `0.05`.

- df:

  Integer. Degrees of freedom for the spline basis. Default is 3.

## Value

A data frame of differentially expressed genes, including:

- AveExpr, F, P.Value, adj.P.Val:

  limma differential expression outputs

- Z1, Z2, Z3:

  Spline coefficients (Z1 typically corresponds to linear trend)

- gene:

  Gene name (from `exp_mat`)

- trend:

  "Positive" or "Negative" trend based on the sign of `Z1`

The first spline coefficient (`Z1`) captures the main expression trend
along the spatial distance.

## Examples

``` r
# Load example data
coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
    package = "SpNeigh"
))
logNorm_expr <- readRDS(system.file("extdata", "LogNormExpr.rds",
    package = "SpNeigh"
))

# Identify cluster-specific cells and compute spatial weights
cells_c0 <- subset(coords, cluster == 0)$cell
bon_c0 <- GetBoundary(data = coords, one_cluster = 0)
weights <- ComputeBoundaryWeights(
    data = coords, cell_ids = cells_c0,
    boundary = bon_c0
)

# Run spatial DE
result <- RunSpatialDE(
    exp_mat = logNorm_expr, cell_ids = cells_c0,
    spatial_distance = weights
)
head(result)
#>                Z1        Z2        Z3  AveExpr         F P.Value adj.P.Val
#> Aldh1a2  90.19564 -72.83701  18.56803 1.642933 1156.7476       0         0
#> Car4    -86.85851  37.72830 -30.91311 3.007433  590.8437       0         0
#> Col1a1   84.54759 -60.30661  21.47013 1.675953  885.1705       0         0
#> Dcn     112.15109 -72.28618  28.47200 2.337544 1295.7527       0         0
#> Fmod     91.46180 -63.54536  20.99449 1.562370 1083.7186       0         0
#> Gjb2     67.07620 -52.43818  14.55533 1.193263  742.6493       0         0
#>            gene    trend
#> Aldh1a2 Aldh1a2 Positive
#> Car4       Car4 Negative
#> Col1a1   Col1a1 Positive
#> Dcn         Dcn Positive
#> Fmod       Fmod Positive
#> Gjb2       Gjb2 Positive
```

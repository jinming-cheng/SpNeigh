# Differential expression analysis between two groups of cells using limma

Performs differential expression analysis between two groups of cells
using the `limma` linear modeling framework. Supports optional
observation-level weights (e.g., spatial weights) and filters genes by
minimum expression threshold across groups.

## Usage

``` r
runLimmaDE(
  exp_mat = NULL,
  cells_reference = NULL,
  cells_target = NULL,
  weights = NULL,
  adj_p.value = 0.05,
  min.pct = 0
)
```

## Arguments

- exp_mat:

  A normalized gene expression matrix (genes x cells), either a `matrix`
  or `dgCMatrix`. Typically log-normalized counts, e.g., from a Seurat
  object.

- cells_reference:

  A character vector of cell IDs to use as the reference (baseline)
  group.

- cells_target:

  A character vector of cell IDs to use as the target (comparison)
  group.

- weights:

  Optional numeric vector of observation-level weights. Must be named
  with cell IDs and match the length of
  `cells_reference + cells_target`.

- adj_p.value:

  Adjusted p-value threshold for reporting differentially expressed
  genes. Default is `0.05`.

- min.pct:

  Minimum proportion of cells expressing the gene in either group
  (values between 0 and 1). Genes not meeting this threshold are
  excluded before testing. Default is `0`.

## Value

A data frame with differentially expressed genes, sorted by absolute log
fold change. Includes columns:

- logFC:

  Log2 fold change of expression (target vs. reference)

- AveExpr:

  Average expression across both groups

- t, P.Value, adj.P.Val, B:

  Statistical results from `limma`

- pct.reference:

  Proportion of reference cells expressing the gene

- pct.target:

  Proportion of target cells expressing the gene

- gene:

  Gene name (rownames from `exp_mat`)

## Examples

``` r
# Load coordinates and log-normalized expression data
coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
    package = "SpNeigh"
))
logNorm_expr <- readRDS(system.file("extdata", "LogNormExpr.rds",
    package = "SpNeigh"
))

# Subset cells from cluster 0 and 2
cells_ref <- subset(coords, cluster == 0)$cell
cells_tar <- subset(coords, cluster == 2)$cell

# Run differential expression with minimum expression threshold
tab <- runLimmaDE(
    exp_mat = logNorm_expr,
    cells_reference = cells_ref,
    cells_target = cells_tar,
    min.pct = 0.25
)

head(tab[, c("gene", "logFC", "adj.P.Val", "pct.reference", "pct.target")])
#>          gene     logFC adj.P.Val pct.reference pct.target
#> Gjc3     Gjc3  4.502163         0    0.25993020  0.9000000
#> Sox10   Sox10  4.147574         0    0.08259930  0.8078571
#> Fn1       Fn1 -3.769727         0    0.79707495  0.1539286
#> Opalin Opalin  3.694963         0    0.06847266  0.6757143
#> Cldn5   Cldn5 -3.472197         0    0.68140269  0.0900000
#> Ly6a     Ly6a -3.290134         0    0.77829483  0.2583929
```

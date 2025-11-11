# Summarize cell counts and proportions inside spatial regions

Computes the number and proportion of cells from each cluster inside
boundaries or ring regions. This function is useful for downstream
visualizations such as bar plots or pie charts showing the spatial
composition of cell types per region.

## Usage

``` r
StatsCellsInside(cells_inside = NULL)
```

## Arguments

- cells_inside:

  An `sf` object of cells returned by
  [`GetCellsInside()`](https://github.com/jinming-cheng/SpNeigh/reference/GetCellsInside.md).
  Must contain `cluster` and `region_id` columns.

## Value

A data frame with one row per cluster per region, containing the
following columns:

- `region_id`: Identifier for each spatial region.

- `cluster`: Cluster label of the cells.

- `count`: Number of cells from the given cluster in the region.

- `proportion`: Proportion of cells from the given cluster relative to
  the total number of cells in the region.

## Examples

``` r
# Load coordinates data
coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
    package = "SpNeigh"
))

# Get boundary and cells inside
boundary_points <- GetBoundary(
    data = coords, one_cluster = 2,
    eps = 120, minPts = 10
)
cells_inside <- GetCellsInside(data = coords, boundary = boundary_points)

# Summarize cluster statistics per region
stats_cells <- StatsCellsInside(cells_inside)
head(stats_cells)
#> # A tibble: 6 × 4
#>   region_id cluster count proportion
#>   <chr>     <fct>   <int>      <dbl>
#> 1 1         0         456   0.102   
#> 2 1         2        3144   0.700   
#> 3 1         3         132   0.0294  
#> 4 1         4         134   0.0299  
#> 5 1         5           3   0.000668
#> 6 1         6          34   0.00757 
```

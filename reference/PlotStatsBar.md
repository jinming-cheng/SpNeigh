# Bar plot of cluster statistics for cells inside boundaries or ring regions

Creates a bar plot to visualize the distribution of cells inside spatial
regions (e.g., boundaries or rings), either as raw counts or proportions
per cluster. The plot is faceted by `region_id` to show statistics
across multiple spatial subregions.

## Usage

``` r
plotStatsBar(
  cell_stats = NULL,
  stat_column = c("proportion", "count"),
  colors = colors15_cheng,
  angle_x_label = 0,
  theme_ggplot = theme_spneigh()
)
```

## Arguments

- cell_stats:

  A data frame containing summarized cell statistics, typically the
  output from
  [`statsCellsInside()`](https://github.com/jinming-cheng/SpNeigh/reference/StatsCellsInside.md).
  Must include columns `region_id`, `cluster`, and the specified
  `stat_column`.

- stat_column:

  Character. Column name in `cell_stats` to use for the y-axis. Options
  are `"count"` (number of cells) or `"proportion"` (relative fraction
  per region).

- colors:

  A vector of cluster colors. Default uses `colors15_cheng`.

- angle_x_label:

  Numeric angle (in degrees) to rotate the x-axis labels. Useful for
  improving label readability in faceted or dense plots. Default is 0
  (no rotation).

- theme_ggplot:

  A ggplot2 theme object. Default is
  [`theme_spneigh()`](https://github.com/jinming-cheng/SpNeigh/reference/theme_spneigh.md).

## Value

A `ggplot2` object showing a faceted bar plot of cell statistics per
region.

## Examples

``` r
coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
    package = "SpNeigh"
))

boundary_points <- getBoundary(
    data = coords, one_cluster = 2,
    eps = 120, minPts = 10
)
cells_inside <- getCellsInside(data = coords, boundary = boundary_points)
stats_cells <- statsCellsInside(cells_inside)

plotStatsBar(stats_cells, stat_column = "proportion")

plotStatsBar(stats_cells, stat_column = "count")

```

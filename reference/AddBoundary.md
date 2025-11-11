# Add boundary outlines to a spatial ggplot

Overlays spatial boundaries (from points or polygons) onto a `ggplot2`
object. The input can be either a data frame of boundary points or an
`sf` object of `POLYGON` geometries. This function is typically used as
an additive layer (`+ AddBoundary(...)`) in conjunction with a base plot
created using
[`PlotBoundary()`](https://github.com/jinming-cheng/SpNeigh/reference/PlotBoundary.md).

## Usage

``` r
AddBoundary(
  boundary = NULL,
  color_boundary = "black",
  linewidth_boundary = 1.5
)
```

## Arguments

- boundary:

  A data frame with columns `x`, `y`, and `region_id` or an `sf` object
  of `POLYGON` or `LINESTRING` geometries.

- color_boundary:

  Color for boundary lines. Default is `"black"`.

- linewidth_boundary:

  Numeric. Line width for boundary outlines. Default is 1.5.

## Value

A
[`ggplot2::geom_path`](https://ggplot2.tidyverse.org/reference/geom_path.html)
layer that can be added to an existing plot.

## Examples

``` r
coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
    package = "SpNeigh"
))

# Automatically get boundary for a cluster
boundary_points <- GetBoundary(
    data = coords, one_cluster = 2,
    subregion_method = "dbscan",
    eps = 120, minPts = 10
)

# Plot with boundary overlay
PlotBoundary(coords) + AddBoundary(boundary_points)

```

# Add boundary polygons or linestrings to a spatial plot

Overlays boundary polygons or linestrings on a spatial `ggplot2` plot.
This function adds an `sf` geometry layer to display complete boundary
shapes for visualizing spatial clusters, rings, or enriched zones.

## Usage

``` r
AddBoundaryPoly(
  boundary_poly,
  color_boundary = "black",
  linewidth_boundary = 1.5
)
```

## Arguments

- boundary_poly:

  An `sf` object containing `POLYGON` or `LINESTRING` geometries and a
  `region_id` column.

- color_boundary:

  Color for boundary lines. Default is `"black"`.

- linewidth_boundary:

  Numeric. Line width for boundary outlines. Default is 1.5.

## Value

A
[`ggplot2::geom_sf`](https://ggplot2.tidyverse.org/reference/ggsf.html)
layer that can be added to an existing plot.

## Examples

``` r
coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
    package = "SpNeigh"
))

boundary_points <- GetBoundary(
    data = coords, one_cluster = 2,
    subregion_method = "dbscan",
    eps = 120, minPts = 10
)
boundary_polys <- BuildBoundaryPoly(boundary_points)

# Add boundary polygons to the plot
PlotBoundary(coords) +
    AddBoundaryPoly(boundary_polys, color_boundary = "blue")

```

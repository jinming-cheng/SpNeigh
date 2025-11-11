# Plot boundary edges or segmented boundary lines

Plots boundary edges as colored `LINESTRING` or `POLYGON` outlines using
`ggplot2`. This function is especially useful for visualizing specific
edges extracted from a polygon using
[`SplitBoundaryPolyByAnchor()`](https://github.com/jinming-cheng/SpNeigh/reference/SplitBoundaryPolyByAnchor.md).

## Usage

``` r
PlotEdge(
  boundary_poly = NULL,
  linewidth_boundary = 1,
  theme_ggplot = my_theme_ggplot(),
  ...
)
```

## Arguments

- boundary_poly:

  An `sf` object containing `POLYGON` or `LINESTRING` geometries and a
  `region_id` column.

- linewidth_boundary:

  Numeric value specifying the line width of the edge. Default is `1`.

- theme_ggplot:

  A ggplot2 theme object. Default is
  [`my_theme_ggplot()`](https://github.com/jinming-cheng/SpNeigh/reference/my_theme_ggplot.md).

- ...:

  Additional arguments passed to
  [`ggplot2::geom_sf()`](https://ggplot2.tidyverse.org/reference/ggsf.html).

## Value

A `ggplot` object displaying the edge outlines colored by `region_id`.

## Examples

``` r
coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
    package = "SpNeigh"
))

# Build boundary polygon and plot its outline
boundary_points <- GetBoundary(
    data = coords, one_cluster = 2,
    eps = 120, minPts = 10
)
boundary_polys <- BuildBoundaryPoly(boundary_points)
PlotEdge(boundary_poly = boundary_polys)


# Split a polygon into edge segments and plot
boundary_edges <- SplitBoundaryPolyByAnchor(boundary_polys[1, ])
PlotEdge(boundary_poly = boundary_edges)

```

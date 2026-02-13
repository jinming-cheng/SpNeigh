# Plot boundary edges or segmented boundary lines

Plots boundary edges as colored `LINESTRING` or `POLYGON` outlines using
`ggplot2`. This function is especially useful for visualizing specific
edges extracted from a polygon using
[`splitBoundaryPolyByAnchor()`](https://github.com/jinming-cheng/SpNeigh/reference/SplitBoundaryPolyByAnchor.md).

## Usage

``` r
plotEdge(
  boundary_poly = NULL,
  linewidth_boundary = 1,
  theme_ggplot = theme_spneigh(),
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
  [`theme_spneigh()`](https://github.com/jinming-cheng/SpNeigh/reference/theme_spneigh.md).

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
boundary_points <- getBoundary(
    data = coords, one_cluster = 2,
    eps = 120, minPts = 10
)
boundary_polys <- buildBoundaryPoly(boundary_points)
plotEdge(boundary_poly = boundary_polys)


# Split a polygon into edge segments and plot
boundary_edges <- splitBoundaryPolyByAnchor(boundary_polys[1, ])
plotEdge(boundary_poly = boundary_edges)

```

# Plot filled spatial regions inside boundaries or rings

Creates a ggplot of spatial regions (e.g., subregions or ring areas)
using filled polygons. Each region is automatically assigned a fill
color based on its `region_id`. This function is commonly used to
visualize the area inside spatial boundaries or surrounding rings, such
as those created by
[`BuildBoundaryPoly()`](https://github.com/jinming-cheng/SpNeigh/reference/BuildBoundaryPoly.md)
or
[`GetRingRegion()`](https://github.com/jinming-cheng/SpNeigh/reference/GetRingRegion.md).

## Usage

``` r
PlotRegion(
  boundary_poly = NULL,
  alpha = 0.5,
  color_boundary = "black",
  linewidth_boundary = 1,
  theme_ggplot = my_theme_ggplot(),
  ...
)
```

## Arguments

- boundary_poly:

  An `sf` object of `POLYGON` geometries containing a `region_id`
  column. Typically created by
  [`BuildBoundaryPoly()`](https://github.com/jinming-cheng/SpNeigh/reference/BuildBoundaryPoly.md)
  or
  [`GetRingRegion()`](https://github.com/jinming-cheng/SpNeigh/reference/GetRingRegion.md).

- alpha:

  Numeric value controlling the transparency of the filled regions.
  Default is `0.5`.

- color_boundary:

  Color of the region outlines. Default is `"black"`.

- linewidth_boundary:

  Numeric line width of the region borders. Default is `1`.

- theme_ggplot:

  A ggplot2 theme object. Default is
  [`my_theme_ggplot()`](https://github.com/jinming-cheng/SpNeigh/reference/my_theme_ggplot.md).

- ...:

  Additional arguments passed to
  [`ggplot2::geom_sf()`](https://ggplot2.tidyverse.org/reference/ggsf.html).

## Value

A `ggplot` object displaying filled spatial regions by `region_id`.

## Examples

``` r
coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
    package = "SpNeigh"
))

# Plot filled boundary regions
boundary_points <- GetBoundary(
    data = coords, one_cluster = 2,
    eps = 120, minPts = 10
)
boundary_polys <- BuildBoundaryPoly(boundary_points)
PlotRegion(boundary_poly = boundary_polys)


# Plot filled ring regions
ring_regions <- GetRingRegion(boundary = boundary_points, dist = 100)
PlotRegion(boundary_poly = ring_regions)

```

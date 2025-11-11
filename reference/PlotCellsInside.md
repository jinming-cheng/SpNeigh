# Plot cells located within spatial boundaries or ring regions

Visualizes cells that fall within defined spatial regions (boundaries or
rings), typically obtained using the
[`GetCellsInside()`](https://github.com/jinming-cheng/SpNeigh/reference/GetCellsInside.md)
function. The cells are colored by cluster, and the function offers two
plotting modes: using `geom_sf()` (with fixed 1:1 aspect ratio) or
`geom_point()` (with more flexible layout).

## Usage

``` r
PlotCellsInside(
  cells_inside = NULL,
  point_size = 0.5,
  colors = my_colors_15,
  theme_ggplot = my_theme_ggplot(),
  legend_size = 2,
  fixed_aspect_ratio = TRUE
)
```

## Arguments

- cells_inside:

  An `sf` object of cells returned by
  [`GetCellsInside()`](https://github.com/jinming-cheng/SpNeigh/reference/GetCellsInside.md).
  Must contain `cluster` and `region_id` columns.

- point_size:

  Numeric. Size of the points representing cells. Default is 0.5.

- colors:

  A vector of cluster colors. Default uses `my_colors_15`.

- theme_ggplot:

  A ggplot2 theme object. Default is
  [`my_theme_ggplot()`](https://github.com/jinming-cheng/SpNeigh/reference/my_theme_ggplot.md).

- legend_size:

  Numeric. Size of legend keys. Default is 2.

- fixed_aspect_ratio:

  Logical. If `TRUE`, uses `geom_sf()` to preserve spatial scale. If
  `FALSE`, uses `geom_point()` with extracted coordinates. Default is
  `TRUE`.

## Value

A `ggplot` object showing the cells colored by cluster within spatial
regions.

## Examples

``` r
coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
    package = "SpNeigh"
))

# Plot cells inside boundaries
boundary_points <- GetBoundary(data = coords, one_cluster = 2)
cells_inside <- GetCellsInside(data = coords, boundary = boundary_points)
PlotCellsInside(cells_inside)


# Plot cells inside rings
ring_regions <- GetRingRegion(boundary = boundary_points, dist = 100)
cells_ring <- GetCellsInside(data = coords, boundary = ring_regions)
PlotCellsInside(cells_ring, fixed_aspect_ratio = FALSE)

```

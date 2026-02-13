# Convert boundary polygons to boundary point coordinates

Extracts the vertex coordinates (`x`, `y`) of boundary geometries from
an `sf` object containing `POLYGON` or `LINESTRING` features. This is
useful for recovering the original boundary points from smoothed or
labeled polygonal regions.

## Usage

``` r
boundaryPolyToPoints(boundary_poly = NULL)
```

## Arguments

- boundary_poly:

  An `sf` object containing only `POLYGON` or `LINESTRING` geometries.
  Must include a `region_id` column for labeling subregions.

## Value

A data frame with columns `x`, `y`, and `region_id`, containing the
vertex coordinates of the boundary geometries.

## Examples

``` r
# Load coordinates and generate boundary
coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
    package = "SpNeigh"
))
boundary_points <- getBoundary(data = coords, one_cluster = 2)
boundary_polys <- buildBoundaryPoly(boundary_points)

# Convert back to boundary points
boundary_pts <- boundaryPolyToPoints(boundary_polys)
#> Warning: repeating attributes for all sub-geometries for which they may not be constant
head(boundary_pts)
#>            x        y region_id
#> 1   1698.220 3196.813         1
#> 1.1 1718.113 3201.256         1
#> 1.2 1718.264 3211.133         1
#> 1.3 1732.328 3227.695         1
#> 1.4 1707.091 3262.311         1
#> 1.5 1698.292 3289.946         1
```

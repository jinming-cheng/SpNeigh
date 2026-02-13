# Generate an outer or inner boundary polygon by buffering an existing boundary

Computes an expanded or shrunken boundary by applying a spatial buffer
to an existing polygon. The input can be either boundary points (as
returned by
[`getBoundary()`](https://github.com/jinming-cheng/SpNeigh/reference/GetBoundary.md))
or polygon geometries (as returned by
[`buildBoundaryPoly()`](https://github.com/jinming-cheng/SpNeigh/reference/BuildBoundaryPoly.md)).
This is useful for defining outer spatial neighborhoods or for shrinking
boundaries inward to define inner regions.

## Usage

``` r
getOuterBoundary(boundary = NULL, dist = 100)
```

## Arguments

- boundary:

  A data frame of boundary points (with columns `x`, `y`, `region_id`)
  or an `sf` object.

- dist:

  A numeric value specifying the buffer distance in spatial units. Use
  positive values to expand (outer boundary), and negative values to
  shrink (inner boundary). Default is 100.

## Value

An `sf` object of expanded outer boundary polygons with the same
`region_id` values as the original input.

## Details

Be careful: when using a negative buffer distance (for inner
boundaries), polygons may collapse, become invalid, or disappear
entirely if the buffer width exceeds the shape's interior size.

## See also

[`getInnerBoundary()`](https://github.com/jinming-cheng/SpNeigh/reference/GetInnerBoundary.md)
for a simplified wrapper for inward shrinking.

## Examples

``` r
# Load coordinates
coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
    package = "SpNeigh"
))

# Get boundary points of cluster 2
boundary_points <- getBoundary(data = coords, one_cluster = 2)

# Build polygons from boundary points
boundary_polys <- buildBoundaryPoly(boundary_points)

# Generate outer boundaries with 100-unit buffer
outer1 <- getOuterBoundary(boundary_points, dist = 100)
outer2 <- getOuterBoundary(boundary_polys, dist = 100)

# Plot original and expanded boundaries
plot(boundary_polys)
plot(outer2, add = TRUE, border = "red")

```

# Generate an inner boundary polygon by shrinking an existing boundary inward

Computes an inward-shifted (contracted) boundary by applying a negative
spatial buffer to an existing boundary polygon. This function is a
wrapper around
[`GetOuterBoundary()`](https://github.com/jinming-cheng/SpNeigh/reference/GetOuterBoundary.md)
and is useful for defining an inner region around a spatial cluster or
structure.

## Usage

``` r
GetInnerBoundary(boundary = NULL, dist = 50)
```

## Arguments

- boundary:

  A data frame of boundary points (with columns `x`, `y`, `region_id`)
  or an `sf` object of `POLYGON` geometries.

- dist:

  A positive numeric value specifying how far inward to shrink the
  boundary. This is automatically converted to a negative buffer
  distance. Default is 50.

## Value

An `sf` object containing the inward-shrunk polygons, one per
`region_id`.

## Details

Be careful: if the shrinkage distance is too large, the resulting
geometry may become invalid or disappear entirely, especially for narrow
or irregular shapes.

## See also

[`GetOuterBoundary()`](https://github.com/jinming-cheng/SpNeigh/reference/GetOuterBoundary.md)
for the outward (positive) buffer version.

## Examples

``` r
# Load coordinates
coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
    package = "SpNeigh"
))

# Get boundary polygons of cluster 2
boundary_points <- GetBoundary(data = coords, one_cluster = 2)
boundary_polys <- BuildBoundaryPoly(boundary_points)
plot(boundary_polys)


# Generate inner boundary with 50-unit buffer for boundary region 1
inner_boundary <- GetInnerBoundary(boundary_polys)
#> Warning: Some inner boundaries may have collapsed or disappeared due to excessive shrinkage.
plot(inner_boundary)

```

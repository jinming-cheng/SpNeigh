# Identify cells located within spatial boundaries or ring regions

Returns the subset of cells from the input data that fall spatially
within a given boundary or ring. The boundary can be provided either as
raw boundary points (from
[`GetBoundary()`](https://github.com/jinming-cheng/SpNeigh/reference/GetBoundary.md)),
as polygons (from
[`BuildBoundaryPoly()`](https://github.com/jinming-cheng/SpNeigh/reference/BuildBoundaryPoly.md)),
or as ring regions (from
[`GetRingRegion()`](https://github.com/jinming-cheng/SpNeigh/reference/GetRingRegion.md)).
The function uses spatial point-in-polygon matching via
[`sf::st_within`](https://r-spatial.github.io/sf/reference/geos_binary_pred.html).

## Usage

``` r
GetCellsInside(data = NULL, cluster_col = NULL, boundary = NULL)
```

## Arguments

- data:

  A `Seurat` object, a `SpatialExperiment` object, or a data frame
  containing spatial coordinates.

- cluster_col:

  Character scalar specifying the metadata column name containing
  cluster assignments. If `NULL`, a default is used depending on the
  input object type:

  - `"seurat_clusters"` for `Seurat` objects

  - `"cluster"` for `SpatialExperiment` objects

- boundary:

  An `sf` object (polygon or ring) or a data frame of boundary points
  returned by
  [`GetBoundary()`](https://github.com/jinming-cheng/SpNeigh/reference/GetBoundary.md),
  [`BuildBoundaryPoly()`](https://github.com/jinming-cheng/SpNeigh/reference/BuildBoundaryPoly.md),
  or
  [`GetRingRegion()`](https://github.com/jinming-cheng/SpNeigh/reference/GetRingRegion.md).

## Value

A data frame (tibble) of cells located inside the given spatial
region(s), with region assignment in a `region_id` column.

## Examples

``` r
# Load coordinates
coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
    package = "SpNeigh"
))
head(coords)
#>          x        y cell cluster
#> 1 1898.815 2540.963    1       4
#> 2 1895.305 2532.627    2       4
#> 3 2368.073 2534.409    3       2
#> 4 1903.726 2560.010    4       4
#> 5 1917.481 2543.132    5       4
#> 6 1926.540 2560.044    6       4

# Get cells inside boundaries
boundary_points <- GetBoundary(
    data = coords, one_cluster = 2,
    eps = 120, minPts = 10
)
cells_inside <- GetCellsInside(data = coords, boundary = boundary_points)
cells_inside
#> Simple feature collection with 5073 features and 3 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 1590.613 ymin: 1640.136 xmax: 5443.607 ymax: 3531.379
#> CRS:           NA
#> First 10 features:
#>    cell cluster region_id                  geometry
#> 3     3       2         1 POINT (2368.073 2534.409)
#> 20   20       2         1 POINT (2456.959 2560.494)
#> 31   31       2         1 POINT (2152.998 2554.538)
#> 43   43       2         1 POINT (2312.183 2555.867)
#> 44   44       2         1  POINT (2323.21 2551.805)
#> 50   50       2         1 POINT (2356.391 2544.839)
#> 52   52       2         1 POINT (2387.844 2537.827)
#> 53   53       2         1 POINT (2401.216 2526.458)
#> 54   54       2         1 POINT (2392.535 2555.696)
#> 55   55       2         1 POINT (2407.538 2545.443)

# Get cells inside rings
ring_regions <- GetRingRegion(boundary = boundary_points, dist = 100)
cells_ring <- GetCellsInside(data = coords, boundary = ring_regions)
cells_ring
#> Simple feature collection with 4362 features and 3 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 1482.997 ymin: 1530.545 xmax: 5441.679 ymax: 3532.463
#> CRS:           NA
#> First 10 features:
#>    cell cluster region_id                  geometry
#> 21   21       4         1 POINT (2080.017 2559.984)
#> 22   22       4         1 POINT (2088.025 2526.528)
#> 23   23       2         1 POINT (2085.119 2546.843)
#> 24   24       2         1 POINT (2101.906 2549.501)
#> 25   25       4         1 POINT (2102.112 2540.185)
#> 26   26       4         1 POINT (2089.812 2539.796)
#> 27   27       4         1  POINT (2095.253 2533.08)
#> 28   28       4         1  POINT (2100.02 2527.252)
#> 29   29       4         1 POINT (2197.507 2524.511)
#> 30   30       3         1 POINT (2209.931 2542.314)
```

# Identify cells located within spatial boundaries or ring regions

Returns the subset of cells from the input data that fall spatially
within a given boundary or ring. The boundary can be provided either as
raw boundary points (from
[`getBoundary()`](https://github.com/jinming-cheng/SpNeigh/reference/GetBoundary.md)),
as polygons (from
[`buildBoundaryPoly()`](https://github.com/jinming-cheng/SpNeigh/reference/BuildBoundaryPoly.md)),
or as ring regions (from
[`getRingRegion()`](https://github.com/jinming-cheng/SpNeigh/reference/GetRingRegion.md)).
The function uses spatial point-in-polygon matching via
[`sf::st_within`](https://r-spatial.github.io/sf/reference/geos_binary_pred.html).

## Usage

``` r
getCellsInside(data = NULL, cluster_col = NULL, boundary = NULL)
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

  - `"cluster"` for `data.frame` objects

- boundary:

  An `sf` object (polygon or ring) or a data frame of boundary points
  returned by
  [`getBoundary()`](https://github.com/jinming-cheng/SpNeigh/reference/GetBoundary.md),
  [`buildBoundaryPoly()`](https://github.com/jinming-cheng/SpNeigh/reference/BuildBoundaryPoly.md),
  or
  [`getRingRegion()`](https://github.com/jinming-cheng/SpNeigh/reference/GetRingRegion.md).

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
boundary_points <- getBoundary(
    data = coords, one_cluster = 2,
    eps = 120, minPts = 10
)

# Select regions of interests if needed (Optional)
boundary_points <- subset(boundary_points, region_id == 2)

cells_inside <- getCellsInside(data = coords, boundary = boundary_points)
cells_inside
#> Simple feature collection with 584 features and 3 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 4665.024 ymin: 2692.795 xmax: 5163.274 ymax: 3317.818
#> CRS:           NA
#> First 10 features:
#>      cell cluster region_id                  geometry
#> 1860 1860       2         2  POINT (4937.139 3062.65)
#> 1863 1863       2         2  POINT (4926.83 3260.598)
#> 1865 1865       2         2 POINT (4959.216 3264.934)
#> 1875 1875       5         2 POINT (4961.168 3277.609)
#> 1880 1880       2         2 POINT (4939.041 3272.613)
#> 1891 1891       2         2 POINT (4955.413 3256.343)
#> 1894 1894       2         2 POINT (4932.086 3195.083)
#> 1897 1897       2         2 POINT (4934.451 3038.603)
#> 1898 1898       2         2  POINT (4942.686 3216.76)
#> 1902 1902       2         2 POINT (4953.803 3244.153)

# Get cells inside rings
ring_regions <- getRingRegion(boundary = boundary_points, dist = 100)
cells_ring <- getCellsInside(data = coords, boundary = ring_regions)
cells_ring
#> Simple feature collection with 1369 features and 3 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 4572.773 ymin: 2600.433 xmax: 5264.996 ymax: 3421.448
#> CRS:           NA
#> First 10 features:
#>      cell cluster region_id                  geometry
#> 1861 1861       5         2 POINT (4933.268 3097.647)
#> 1862 1862       2         2 POINT (4934.709 3075.913)
#> 1864 1864       5         2 POINT (4954.823 3070.896)
#> 1866 1866       5         2 POINT (4951.338 3054.113)
#> 1867 1867       5         2 POINT (4952.474 3130.228)
#> 1868 1868       5         2 POINT (4958.769 3086.888)
#> 1869 1869       5         2  POINT (4951.975 3185.69)
#> 1870 1870       5         2 POINT (4950.835 3113.254)
#> 1871 1871       5         2 POINT (4949.809 3047.716)
#> 1872 1872       2         2  POINT (4966.83 3247.072)
```

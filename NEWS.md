
# SpNeigh 0.99.43

- Updated CITATION file
- Added Bioconductor installation instructions to the vignette and README


# SpNeigh 0.99.42

- Updated `extractCoords.data.frame()` implementation


# SpNeigh 0.99.41

- Renamed functions to follow lowerCamelCase convention
- Renamed `my_theme_ggplot()` to `theme_spneigh()`
- Renamed `my_color_15` to `colors15_cheng`
- Implemented `extractCoords()` as an S3 generic
- Updated faceting implementation


# SpNeigh 0.99.2

- Added CITATION file


# SpNeigh 0.99.1

- Made functions compatible with `SpatialExperiment` objects
- Added script in `inst/scripts` to document data in `inst/extdata`


# SpNeigh 0.99.0

*Initial Bioconductor submission*

## New Features

- Introduced the `RunSpatialDE()` function for spatial differential expression using spatial weights.
- Added `ComputeSpatialEnrichmentIndex()` for enrichment quantification.
- Implemented `GetBoundary()`, `GetRingRegion()`, and `RemoveOutliers()` for spatial region detection.
- Added `PlotSpatialExpression()` and other visualization tools.


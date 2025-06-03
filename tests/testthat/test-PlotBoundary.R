test_that("Test PlotBoundary related functions", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    coords$cluster <- as.character(coords$cluster)

    expect_silent(PlotBoundary(coords, one_cluster = 2, split_by = "cluster"))

    expect_silent(PlotBoundary(coords, one_cluster = 2, sub_plot = TRUE))

    expect_error(PlotBoundary(coords, one_cluster = NULL, sub_plot = TRUE))

    boundary_points <- GetBoundary(
        data = coords, one_cluster = 2
    )

    boundary_polys <- BuildBoundaryPoly(boundary_points)

    expect_silent(AddBoundary(boundary_points))

    expect_warning(AddBoundary(boundary_polys))

    expect_error(AddBoundary(boundary = NULL))

    expect_error(AddBoundary(boundary = boundary_points[, c("x", "y")]))

    expect_silent(AddBoundaryPoly(boundary_poly = boundary_polys))

    expect_error(AddBoundaryPoly(boundary_poly = boundary_points))

    ring_regions <- GetRingRegion(boundary = boundary_points)

    ring_regions$region_id <- as.character(ring_regions$region_id)

    expect_silent(PlotRegion(boundary_poly = ring_regions))

    expect_error(PlotRegion(boundary_poly = NULL))

    expect_silent(PlotEdge(boundary_poly = boundary_polys))

    boundary_edges <- SplitBoundaryPolyByAnchor(boundary_polys[1, ])

    boundary_edges$region_id <- as.character(boundary_edges$region_id)

    expect_silent(PlotEdge(boundary_edges))

    expect_error(PlotEdge(boundary_poly = NULL))
})

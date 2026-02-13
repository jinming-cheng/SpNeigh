test_that("Test plotBoundary related functions", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    coords$cluster <- as.character(coords$cluster)

    expect_silent(plotBoundary(coords, one_cluster = 2, split_by = "cluster"))

    expect_silent(plotBoundary(coords, one_cluster = 2, sub_plot = TRUE))

    expect_error(plotBoundary(coords, one_cluster = NULL, sub_plot = TRUE))

    boundary_points <- getBoundary(
        data = coords, one_cluster = 2
    )

    boundary_polys <- buildBoundaryPoly(boundary_points)

    expect_silent(addBoundary(boundary_points))

    expect_warning(addBoundary(boundary_polys))

    expect_error(addBoundary(boundary = NULL))

    expect_error(addBoundary(boundary = boundary_points[, c("x", "y")]))

    expect_silent(addBoundaryPoly(boundary_poly = boundary_polys))

    expect_error(addBoundaryPoly(boundary_poly = boundary_points))

    ring_regions <- getRingRegion(boundary = boundary_points)

    ring_regions$region_id <- as.character(ring_regions$region_id)

    expect_silent(plotRegion(boundary_poly = ring_regions))

    expect_error(plotRegion(boundary_poly = NULL))

    expect_silent(plotEdge(boundary_poly = boundary_polys))

    boundary_edges <- splitBoundaryPolyByAnchor(boundary_polys[1, ])

    boundary_edges$region_id <- as.character(boundary_edges$region_id)

    expect_silent(plotEdge(boundary_edges))

    expect_error(plotEdge(boundary_poly = NULL))
})

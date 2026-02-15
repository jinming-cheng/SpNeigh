test_that("Test plotBoundary related functions", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    coords$cluster <- as.character(coords$cluster)

    p <- plotBoundary(coords, one_cluster = 2, split_by = "cluster")
    expect_s3_class(p, "ggplot")

    p <- plotBoundary(coords, one_cluster = 2, sub_plot = TRUE)
    expect_s3_class(p, "ggplot")

    expect_error(plotBoundary(coords, one_cluster = NULL, sub_plot = TRUE))

    boundary_points <- getBoundary(
        data = coords, one_cluster = 2
    )

    boundary_polys <- buildBoundaryPoly(boundary_points)

    p_layer <- addBoundary(boundary_points)

    p <- ggplot2::ggplot() + p_layer
    expect_s3_class(p, "ggplot")

    expect_warning(addBoundary(boundary_polys))

    expect_error(addBoundary(boundary = NULL))

    expect_error(addBoundary(boundary = boundary_points[, c("x", "y")]))

    p_layer <- addBoundaryPoly(boundary_poly = boundary_polys)

    p <- ggplot2::ggplot() + p_layer
    expect_s3_class(p, "ggplot")

    expect_error(addBoundaryPoly(boundary_poly = boundary_points))

    ring_regions <- getRingRegion(boundary = boundary_points)

    ring_regions$region_id <- as.character(ring_regions$region_id)

    p <- plotRegion(boundary_poly = ring_regions)
    expect_s3_class(p, "ggplot")

    expect_error(plotRegion(boundary_poly = NULL))

    p <- plotEdge(boundary_poly = boundary_polys)
    expect_s3_class(p, "ggplot")

    boundary_edges <- splitBoundaryPolyByAnchor(boundary_polys[1, ])

    boundary_edges$region_id <- as.character(boundary_edges$region_id)

    p <- plotEdge(boundary_edges)
    expect_s3_class(p, "ggplot")

    expect_error(plotEdge(boundary_poly = NULL))
})

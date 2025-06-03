test_that("Test ComputeCentroidWeights", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    cells_c2 <- subset(coords, cluster == 2)$cell

    weights <- ComputeCentroidWeights(data = coords, cell_ids = cells_c2)

    expect_true(is.vector(weights))

    expect_true(is.vector(names(weights)))

    expect_error(ComputeCentroidWeights(
        data = coords[, c("x", "y")],
        cell_ids = cells_c2
    ))

    expect_error(ComputeCentroidWeights(
        data = coords,
        cell_ids = NULL
    ))
})



test_that("Test ComputeBoundaryWeights", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    boundary_points <- GetBoundary(data = coords, one_cluster = 2)

    boundary_polys <- BuildBoundaryPoly(boundary_points)

    cells_c2 <- subset(coords, cluster == 2)$cell

    weights <- ComputeBoundaryWeights(
        data = coords, cell_ids = cells_c2,
        boundary = boundary_points
    )

    expect_true(is.vector(weights))

    expect_true(is.vector(names(weights)))

    expect_error(ComputeBoundaryWeights(
        data = coords[, c("x", "y")], cell_ids = cells_c2,
        boundary = boundary_polys[1, ]
    ))

    expect_error(ComputeBoundaryWeights(
        data = coords, cell_ids = cells_c2,
        boundary = "not_a_boundary"
    ))

    expect_error(ComputeBoundaryWeights(
        data = coords, cell_ids = NULL,
        boundary = boundary_polys[1, ]
    ))


    boundary_edges <- SplitBoundaryPolyByAnchor(boundary_polys[1, ])

    weights_edge <- ComputeBoundaryWeights(
        data = coords, cell_ids = cells_c2,
        boundary = boundary_edges[2, ]
    )

    expect_true(is.vector(weights))
})

test_that("Test computeCentroidWeights", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    cells_c2 <- subset(coords, cluster == 2)$cell

    weights <- computeCentroidWeights(data = coords, cell_ids = cells_c2)

    expect_true(is.vector(weights))

    expect_true(is.vector(names(weights)))

    expect_error(computeCentroidWeights(
        data = coords[, c("x", "y")],
        cell_ids = cells_c2
    ))

    expect_error(computeCentroidWeights(
        data = coords,
        cell_ids = NULL
    ))
})



test_that("Test computeBoundaryWeights", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    boundary_points <- getBoundary(data = coords, one_cluster = 2)

    boundary_polys <- buildBoundaryPoly(boundary_points)

    cells_c2 <- subset(coords, cluster == 2)$cell

    weights <- computeBoundaryWeights(
        data = coords, cell_ids = cells_c2,
        boundary = boundary_points
    )

    expect_true(is.vector(weights))

    expect_true(is.vector(names(weights)))

    expect_error(computeBoundaryWeights(
        data = coords[, c("x", "y")], cell_ids = cells_c2,
        boundary = boundary_polys[1, ]
    ))

    expect_error(computeBoundaryWeights(
        data = coords, cell_ids = cells_c2,
        boundary = "not_a_boundary"
    ))

    expect_error(computeBoundaryWeights(
        data = coords, cell_ids = NULL,
        boundary = boundary_polys[1, ]
    ))


    boundary_edges <- splitBoundaryPolyByAnchor(boundary_polys[1, ])

    weights_edge <- computeBoundaryWeights(
        data = coords, cell_ids = cells_c2,
        boundary = boundary_edges[2, ]
    )

    expect_true(is.vector(weights))
})

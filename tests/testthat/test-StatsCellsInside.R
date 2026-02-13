test_that("Test statsCellsInside", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    boundary_points <- getBoundary(
        data = coords, one_cluster = 2, eps = 120, minPts = 10
    )

    boundary_polys <- buildBoundaryPoly(boundary_points)

    cells_inside <- getCellsInside(data = coords, boundary = boundary_polys[2, ])

    stats_cells <- statsCellsInside(cells_inside)

    expect_true(is.data.frame(stats_cells))
})

test_that("Test StatsCellsInside", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    boundary_points <- GetBoundary(
        data = coords, one_cluster = 2, eps = 120, minPts = 10
    )

    boundary_polys <- BuildBoundaryPoly(boundary_points)

    cells_inside <- GetCellsInside(data = coords, boundary = boundary_polys[2, ])

    stats_cells <- StatsCellsInside(cells_inside)

    expect_true(is.data.frame(stats_cells))
})

test_that("Test getCellsInside", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    boundary_points <- getBoundary(
        data = coords, one_cluster = 2, eps = 120, minPts = 10
    )

    boundary_points_sub <- subset(boundary_points, region_id == "2")

    cells_inside <- getCellsInside(data = coords, boundary = boundary_points_sub)

    expect_true(inherits(cells_inside, "sf"))
    expect_true(nrow(cells_inside) > 0)

    boundary_polys <- buildBoundaryPoly(boundary_points)

    cells_inside <- getCellsInside(data = coords, boundary = boundary_polys[2, ])

    expect_true(inherits(cells_inside, "sf"))
    expect_true(nrow(cells_inside) > 0)
})

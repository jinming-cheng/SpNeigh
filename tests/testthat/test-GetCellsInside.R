test_that("Test GetCellsInside", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    boundary_points <- GetBoundary(
        data = coords, one_cluster = 2, eps = 120, minPts = 10
    )

    boundary_points_sub <- subset(boundary_points, region_id == "2")

    cells_inside <- GetCellsInside(data = coords, boundary = boundary_points_sub)

    expect_true(inherits(cells_inside, "sf"))

    boundary_polys <- BuildBoundaryPoly(boundary_points)

    cells_inside <- GetCellsInside(data = coords, boundary = boundary_polys[2, ])

    expect_true(inherits(cells_inside, "sf"))
})

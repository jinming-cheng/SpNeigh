test_that("Test PlotCellsInside", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    boundary_points <- GetBoundary(data = coords, one_cluster = 2)

    boundary_polys <- BuildBoundaryPoly(boundary_points)

    cells_inside <- GetCellsInside(data = coords, boundary = boundary_polys[2, ])

    expect_silent(PlotCellsInside(cells_inside))

    cells_inside$cluster <- as.character(cells_inside$cluster)

    expect_silent(PlotCellsInside(cells_inside, fixed_aspect_ratio = FALSE))
})

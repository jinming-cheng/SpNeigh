test_that("Test plotCellsInside", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    boundary_points <- getBoundary(data = coords, one_cluster = 2)

    boundary_polys <- buildBoundaryPoly(boundary_points)

    cells_inside <- getCellsInside(data = coords, boundary = boundary_polys[2, ])

    p <- expect_s3_class(
        plotCellsInside(cells_inside),
        "ggplot"
    )
    expect_gt(length(p$layers), 0)

    cells_inside$cluster <- as.character(cells_inside$cluster)

    p <- expect_s3_class(
        plotCellsInside(cells_inside, fixed_aspect_ratio = FALSE),
        "ggplot"
    )
    expect_gt(length(p$layers), 0)
})

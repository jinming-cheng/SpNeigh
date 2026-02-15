test_that("Test PlotStats functions", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    boundary_points <- getBoundary(
        data = coords, one_cluster = 2,
        eps = 120, minPts = 10
    )

    boundary_polys <- buildBoundaryPoly(boundary_points)

    cells_inside <- getCellsInside(data = coords, boundary = boundary_polys[2, ])

    stats_cells <- statsCellsInside(cells_inside)

    p <- plotStatsBar(stats_cells, stat_column = "proportion")

    expect_s3_class(p, "ggplot")

    stats_cells$cluster <- as.character(stats_cells$cluster)

    p <- plotStatsBar(stats_cells, stat_column = "count")

    expect_s3_class(p, "ggplot")

    p <- plotStatsPie(stats_cells, plot_donut = TRUE)

    expect_s3_class(p, "ggplot")
})

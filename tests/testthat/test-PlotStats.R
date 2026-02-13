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

    expect_silent(plotStatsBar(stats_cells, stat_column = "proportion"))


    stats_cells$cluster <- as.character(stats_cells$cluster)

    expect_silent(plotStatsBar(stats_cells, stat_column = "count"))

    expect_silent(plotStatsPie(stats_cells, plot_donut = TRUE))
})

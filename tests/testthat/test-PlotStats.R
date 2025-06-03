test_that("Test PlotStats functions", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    boundary_points <- GetBoundary(
        data = coords, one_cluster = 2,
        eps = 120, minPts = 10
    )

    boundary_polys <- BuildBoundaryPoly(boundary_points)

    cells_inside <- GetCellsInside(data = coords, boundary = boundary_polys[2, ])

    stats_cells <- StatsCellsInside(cells_inside)

    expect_silent(PlotStatsBar(stats_cells, stat_column = "proportion"))


    stats_cells$cluster <- as.character(stats_cells$cluster)

    expect_silent(PlotStatsBar(stats_cells, stat_column = "count"))

    expect_silent(PlotStatsPie(stats_cells, plot_donut = TRUE))
})

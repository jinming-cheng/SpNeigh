test_that("Test PlotWeights", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    cells_c2 <- subset(coords, cluster == 2)[, "cell"]

    weights_cen <- ComputeCentroidWeights(data = coords, cell_ids = cells_c2)

    expect_silent(PlotWeights(data = coords, weights = weights_cen))

    expect_error(PlotWeights(data = coords[, c("x", "y")], weights = weights_cen))

    expect_error(PlotWeights(data = coords, weights = as.character(weights_cen)))

    new_weights <- weights_cen

    names(new_weights) <- paste0("new", names(new_weights))

    expect_error(PlotWeights(data = coords, weights = new_weights))
})

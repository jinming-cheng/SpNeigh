test_that("Test plotWeights", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    cells_c2 <- subset(coords, cluster == 2)[, "cell"]

    weights_cen <- computeCentroidWeights(data = coords, cell_ids = cells_c2)

    p <- plotWeights(data = coords, weights = weights_cen)

    expect_s3_class(p, "ggplot")

    expect_error(plotWeights(data = coords[, c("x", "y")], weights = weights_cen))

    expect_error(plotWeights(data = coords, weights = as.character(weights_cen)))

    new_weights <- weights_cen

    names(new_weights) <- paste0("new", names(new_weights))

    expect_error(plotWeights(data = coords, weights = new_weights))
})

test_that("Test ComputeSpatialEnrichmentIndex", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    logNorm_expr <- readRDS(system.file("extdata", "LogNormExpr.rds",
        package = "SpNeigh"
    ))

    bon_c0 <- GetBoundary(data = coords, one_cluster = 0)

    cells_c0 <- subset(coords, cluster == 0)$cell

    weights <- ComputeBoundaryWeights(
        data = coords,
        cell_ids = cells_c0,
        boundary = bon_c0
    )

    exp_mat <- logNorm_expr[, cells_c0]

    sei_df <- ComputeSpatialEnrichmentIndex(exp_mat, weights)

    expect_true(is.data.frame(sei_df))

    sei_df <- ComputeSpatialEnrichmentIndex(as.matrix(exp_mat), weights)

    expect_true(is.data.frame(sei_df))

    expect_error(ComputeSpatialEnrichmentIndex("non_matrix", weights))

    expect_error(ComputeSpatialEnrichmentIndex(exp_mat, weights[1:10]))

    expect_error(ComputeSpatialEnrichmentIndex(exp_mat, weights = "non_numeric"))
})

test_that("Test computeSpatialEnrichmentIndex", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    logNorm_expr <- readRDS(system.file("extdata", "LogNormExpr.rds",
        package = "SpNeigh"
    ))

    bon_c0 <- getBoundary(data = coords, one_cluster = 0)

    cells_c0 <- subset(coords, cluster == 0)$cell

    weights <- computeBoundaryWeights(
        data = coords,
        cell_ids = cells_c0,
        boundary = bon_c0
    )

    exp_mat <- logNorm_expr[, cells_c0]

    sei_df <- computeSpatialEnrichmentIndex(exp_mat, weights)

    expect_true(is.data.frame(sei_df))

    sei_df <- computeSEI(exp_mat, weights)

    expect_true(is.data.frame(sei_df))

    sei_df <- computeSpatialEnrichmentIndex(as.matrix(exp_mat), weights)

    expect_true(is.data.frame(sei_df))

    expect_error(computeSpatialEnrichmentIndex("non_matrix", weights))

    expect_error(computeSpatialEnrichmentIndex(exp_mat, weights[1:10]))

    expect_error(computeSpatialEnrichmentIndex(exp_mat, weights = "non_numeric"))
})

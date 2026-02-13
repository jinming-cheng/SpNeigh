test_that("Test runLimmaDE", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    logNorm_expr <- readRDS(system.file("extdata", "LogNormExpr.rds",
        package = "SpNeigh"
    ))

    cells_ref <- subset(coords, cluster == 0)$cell

    cells_tar <- subset(coords, cluster == 2)$cell

    tab <- runLimmaDE(
        exp_mat = as.matrix(logNorm_expr),
        cells_reference = cells_ref,
        cells_target = cells_tar,
        min.pct = 0.25
    )

    expect_true(is.data.frame(tab))

    expect_error(runLimmaDE(
        exp_mat = "non_matrix",
        cells_reference = cells_ref,
        cells_target = cells_tar,
        min.pct = 0.25
    ))

    expect_error(runLimmaDE(
        exp_mat = logNorm_expr,
        cells_reference = NULL,
        cells_target = cells_tar,
        min.pct = 0.25
    ))

    expect_error(runLimmaDE(
        exp_mat = logNorm_expr,
        cells_reference = c(cells_ref, "unknown_cells"),
        cells_target = cells_tar,
        min.pct = 0.25
    ))

    expect_error(runLimmaDE(
        exp_mat = logNorm_expr,
        cells_reference = cells_ref,
        cells_target = cells_tar,
        weights = 1,
        min.pct = 0.25
    ))

    expect_warning(runLimmaDE(
        exp_mat = logNorm_expr,
        cells_reference = cells_ref,
        cells_target = cells_tar,
        min.pct = 0.25,
        adj_p.value = -1
    ))


    expect_error(runLimmaDE(
        exp_mat = logNorm_expr,
        cells_reference = cells_ref,
        cells_target = cells_tar,
        min.pct = 1
    ))
})


test_that("Test runSpatialDE", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    logNorm_expr <- readRDS(system.file("extdata", "LogNormExpr.rds",
        package = "SpNeigh"
    ))

    cells_c0 <- subset(coords, cluster == 0)$cell

    bon_c0 <- getBoundary(data = coords, one_cluster = 0)

    weights <- computeBoundaryWeights(
        data = coords, cell_ids = cells_c0,
        boundary = bon_c0
    )

    result <- runSpatialDE(
        exp_mat = logNorm_expr, cell_ids = cells_c0,
        spatial_distance = weights
    )

    expect_true(is.data.frame(result))

    expect_error(runSpatialDE(
        exp_mat = "non_matrix", cell_ids = cells_c0,
        spatial_distance = weights
    ))

    expect_error(runSpatialDE(
        exp_mat = logNorm_expr, cell_ids = NULL,
        spatial_distance = weights
    ))

    expect_error(runSpatialDE(
        exp_mat = logNorm_expr, cell_ids = c(cells_c0, "unknown_cells"),
        spatial_distance = weights
    ))

    expect_error(runSpatialDE(
        exp_mat = logNorm_expr, cell_ids = cells_c0,
        spatial_distance = weights[1:10]
    ))

    expect_warning(runSpatialDE(
        exp_mat = logNorm_expr, cell_ids = cells_c0,
        spatial_distance = weights, adj_p.value = -1,
    ))
})


test_that("Test splineDesign", {
    x <- seq(1, 0, length.out = 100)

    Z <- splineDesign(x, df = 3)

    k <- ncol(Z)

    expect_true(is.matrix(Z))

    expect_equal(k, 3)

    expect_error(splineDesign(x = "a"))
})

test_that("Test RunLimmaDE", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    logNorm_expr <- readRDS(system.file("extdata", "LogNormExpr.rds",
        package = "SpNeigh"
    ))

    cells_ref <- subset(coords, cluster == 0)$cell

    cells_tar <- subset(coords, cluster == 2)$cell

    tab <- RunLimmaDE(
        exp_mat = as.matrix(logNorm_expr),
        cells_reference = cells_ref,
        cells_target = cells_tar,
        min.pct = 0.25
    )

    expect_true(is.data.frame(tab))

    expect_error(RunLimmaDE(
        exp_mat = "non_matrix",
        cells_reference = cells_ref,
        cells_target = cells_tar,
        min.pct = 0.25
    ))

    expect_error(RunLimmaDE(
        exp_mat = logNorm_expr,
        cells_reference = NULL,
        cells_target = cells_tar,
        min.pct = 0.25
    ))

    expect_error(RunLimmaDE(
        exp_mat = logNorm_expr,
        cells_reference = c(cells_ref, "unknown_cells"),
        cells_target = cells_tar,
        min.pct = 0.25
    ))

    expect_error(RunLimmaDE(
        exp_mat = logNorm_expr,
        cells_reference = cells_ref,
        cells_target = cells_tar,
        weights = 1,
        min.pct = 0.25
    ))

    expect_warning(RunLimmaDE(
        exp_mat = logNorm_expr,
        cells_reference = cells_ref,
        cells_target = cells_tar,
        min.pct = 0.25,
        adj_p.value = -1
    ))


    expect_error(RunLimmaDE(
        exp_mat = logNorm_expr,
        cells_reference = cells_ref,
        cells_target = cells_tar,
        min.pct = 1
    ))
})


test_that("Test RunSpatialDE", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    logNorm_expr <- readRDS(system.file("extdata", "LogNormExpr.rds",
        package = "SpNeigh"
    ))

    cells_c0 <- subset(coords, cluster == 0)$cell

    bon_c0 <- GetBoundary(data = coords, one_cluster = 0)

    weights <- ComputeBoundaryWeights(
        data = coords, cell_ids = cells_c0,
        boundary = bon_c0
    )

    result <- RunSpatialDE(
        exp_mat = logNorm_expr, cell_ids = cells_c0,
        spatial_distance = weights
    )

    expect_true(is.data.frame(result))

    expect_error(RunSpatialDE(
        exp_mat = "non_matrix", cell_ids = cells_c0,
        spatial_distance = weights
    ))

    expect_error(RunSpatialDE(
        exp_mat = logNorm_expr, cell_ids = NULL,
        spatial_distance = weights
    ))

    expect_error(RunSpatialDE(
        exp_mat = logNorm_expr, cell_ids = c(cells_c0, "unknown_cells"),
        spatial_distance = weights
    ))

    expect_error(RunSpatialDE(
        exp_mat = logNorm_expr, cell_ids = cells_c0,
        spatial_distance = weights[1:10]
    ))

    expect_warning(RunSpatialDE(
        exp_mat = logNorm_expr, cell_ids = cells_c0,
        spatial_distance = weights, adj_p.value = -1,
    ))
})


test_that("Test SplineDesign", {
    x <- seq(1, 0, length.out = 100)

    Z <- SplineDesign(x, df = 3)

    k <- ncol(Z)

    expect_true(is.matrix(Z))

    expect_equal(k, 3)

    expect_error(SplineDesign(x = "a"))
})

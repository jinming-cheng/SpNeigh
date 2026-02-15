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

test_that("runLimmaDE detects correct DE direction", {
    set.seed(123)

    # 6 cells total
    cells_reference <- paste0("r", 1:3)
    cells_target <- paste0("t", 1:3)
    all_cells <- c(cells_reference, cells_target)

    # gene1: higher in target
    gene1 <- c(1, 1, 1, 5, 5, 5)

    # gene2: higher in reference
    gene2 <- c(5, 5, 5, 1, 1, 1)

    # gene3: equal
    gene3 <- rep(2, 6)

    exp_mat <- rbind(
        gene1 = gene1,
        gene2 = gene2,
        gene3 = gene3
    )

    colnames(exp_mat) <- all_cells

    res <- runLimmaDE(
        exp_mat = exp_mat,
        cells_reference = cells_reference,
        cells_target = cells_target,
        adj_p.value = 1, # keep all genes
        min.pct = 0
    )

    expect_s3_class(res, "data.frame")
    expect_true(all(c("gene", "logFC") %in% colnames(res)))

    res_named <- res
    rownames(res_named) <- res_named$gene

    # gene1 should have positive logFC (target > reference)
    expect_true(res_named["gene1", "logFC"] > 0)

    # gene2 should have negative logFC
    expect_true(res_named["gene2", "logFC"] < 0)
})


test_that("runLimmaDE min.pct filter works", {
    cells_reference <- c("r1", "r2")
    cells_target <- c("t1", "t2")
    all_cells <- c(cells_reference, cells_target)

    # gene1 expressed only in 1 cell
    gene1 <- c(0, 0, 1, 0)

    # gene2 expressed in all cells
    gene2 <- c(1, 1, 1, 1)

    exp_mat <- rbind(gene1 = gene1, gene2 = gene2)
    colnames(exp_mat) <- all_cells

    res <- suppressWarnings(runLimmaDE(
        exp_mat = exp_mat,
        cells_reference = cells_reference,
        cells_target = cells_target,
        adj_p.value = 1,
        min.pct = 0.51
    ))

    expect_s3_class(res, "data.frame")
    expect_true(all(c("logFC", "AveExpr") %in% colnames(res)))

    expect_equal(nrow(res), 1)
    expect_false("gene1" %in% res$gene)
})


test_that("runLimmaDE warns when no genes pass threshold", {
    cells_reference <- c("r1", "r2")
    cells_target <- c("t1", "t2")
    all_cells <- c(cells_reference, cells_target)

    exp_mat <- matrix(rnorm(4 * 4), nrow = 4)
    rownames(exp_mat) <- paste0("gene", 1:4)
    colnames(exp_mat) <- all_cells

    expect_warning(
        res <- runLimmaDE(
            exp_mat = exp_mat,
            cells_reference = cells_reference,
            cells_target = cells_target,
            adj_p.value = 1e-10
        )
    )

    expect_equal(nrow(res), 0)
})


test_that("runLimmaDE accepts weights", {
    cells_reference <- c("r1", "r2")
    cells_target <- c("t1", "t2")
    all_cells <- c(cells_reference, cells_target)

    exp_mat <- matrix(rnorm(4 * 4), nrow = 4)
    rownames(exp_mat) <- paste0("gene", 1:4)
    colnames(exp_mat) <- all_cells

    weights <- setNames(c(1, 1, 2, 2), all_cells)

    expect_silent(
        runLimmaDE(
            exp_mat = exp_mat,
            cells_reference = cells_reference,
            cells_target = cells_target,
            weights = weights,
            adj_p.value = 1
        )
    )
})


test_that("runLimmaDE input validation works", {
    exp_mat <- matrix(rnorm(4 * 4), nrow = 4)
    colnames(exp_mat) <- c("c1", "c2", "c3", "c4")

    expect_error(
        runLimmaDE(
            exp_mat = NULL,
            cells_reference = "c1",
            cells_target = "c2"
        )
    )

    expect_error(
        runLimmaDE(
            exp_mat = exp_mat,
            cells_reference = NULL,
            cells_target = "c2"
        )
    )

    expect_error(
        runLimmaDE(
            exp_mat = exp_mat,
            cells_reference = "x",
            cells_target = "y"
        )
    )
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

test_that("runSpatialDE detects spatial trends correctly", {
    set.seed(123)

    # 10 cells
    cell_ids <- paste0("c", 1:10)
    spatial_distance <- setNames(seq(0, 1, length.out = 10), cell_ids)

    # gene1: positive trend
    gene1 <- spatial_distance * 5

    # gene2: negative trend
    gene2 <- rev(spatial_distance) * 5

    # gene3: no trend
    gene3 <- rep(1, 10)

    exp_mat <- rbind(
        gene1 = gene1,
        gene2 = gene2,
        gene3 = gene3
    )

    colnames(exp_mat) <- cell_ids

    res <- runSpatialDE(
        exp_mat = exp_mat,
        cell_ids = cell_ids,
        spatial_distance = spatial_distance,
        adj_p.value = 1, # keep all genes
        df = 3
    )

    expect_s3_class(res, "data.frame")

    # Check required columns exist
    expect_true(all(c("gene", "trend") %in% colnames(res)))

    # Check gene names returned
    expect_true(all(c("gene1", "gene2") %in% res$gene))

    # Check trend direction
    res_named <- res
    rownames(res_named) <- res_named$gene

    expect_equal(res_named["gene1", "trend"], "Positive")
    expect_equal(res_named["gene2", "trend"], "Negative")
})

test_that("runSpatialDE filters genes by adjusted p-value", {
    cell_ids <- paste0("c", 1:6)
    spatial_distance <- setNames(seq(0, 1, length.out = 6), cell_ids)

    exp_mat <- matrix(rnorm(3 * 6), nrow = 3)
    rownames(exp_mat) <- paste0("gene", 1:3)
    colnames(exp_mat) <- cell_ids

    expect_warning(
        res <- runSpatialDE(
            exp_mat = exp_mat,
            cell_ids = cell_ids,
            spatial_distance = spatial_distance,
            adj_p.value = 1e-10
        )
    )

    expect_equal(nrow(res), 0)
})

test_that("runSpatialDE input validation works", {
    exp_mat <- matrix(rnorm(6), nrow = 2)
    colnames(exp_mat) <- c("c1", "c2", "c3")

    spatial_distance <- c(c1 = 0, c2 = 1, c3 = 2)

    # Missing cell_ids
    expect_error(
        runSpatialDE(exp_mat = exp_mat, spatial_distance = spatial_distance)
    )

    # Wrong spatial_distance length
    expect_error(
        runSpatialDE(
            exp_mat = exp_mat,
            cell_ids = colnames(exp_mat),
            spatial_distance = c(1, 2)
        )
    )

    # Cell IDs not found
    expect_error(
        runSpatialDE(
            exp_mat = exp_mat,
            cell_ids = c("a", "b", "c"),
            spatial_distance = spatial_distance
        )
    )
})



test_that("Test splineDesign", {
    x <- seq(1, 0, length.out = 100)

    Z <- splineDesign(x, df = 3)

    k <- ncol(Z)

    expect_true(is.matrix(Z))

    expect_equal(k, 3)

    expect_error(splineDesign(x = "a"))
})

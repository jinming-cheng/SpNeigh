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


test_that("computeSpatialEnrichmentIndex computes correct SEI and normalized SEI", {
    exp_mat <- matrix(
        c(
            1, 2, 3,
            4, 5, 6
        ),
        nrow = 2,
        byrow = TRUE
    )
    # genes x cells
    rownames(exp_mat) <- c("gene1", "gene2")
    colnames(exp_mat) <- c("cell1", "cell2", "cell3")

    weights <- c(0.2, 0.3, 0.5)

    res <- computeSpatialEnrichmentIndex(exp_mat, weights)

    # Expect columns and gene names preserved
    expect_true(all(c("gene", "SEI", "mean_expr", "normalized_SEI") %in%
        colnames(res)))
    expect_setequal(res$gene, c("gene1", "gene2"))

    # Manual expected values
    w <- weights / sum(weights)

    expected_sei <- c(
        gene1 = 1 * w[1] + 2 * w[2] + 3 * w[3],
        gene2 = 4 * w[1] + 5 * w[2] + 6 * w[3]
    )

    expected_mean <- c(
        gene1 = mean(c(1, 2, 3)),
        gene2 = mean(c(4, 5, 6))
    )

    expected_norm <- expected_sei / (expected_mean + 1e-6)

    # Match output by gene name (robust to sorting)
    res_named <- res
    rownames(res_named) <- res_named$gene

    expect_equal(
        unname(res_named["gene1", "SEI"]),
        unname(expected_sei["gene1"]),
        tolerance = 1e-12
    )
    expect_equal(
        unname(res_named["gene2", "SEI"]),
        unname(expected_sei["gene2"]),
        tolerance = 1e-12
    )

    expect_equal(
        unname(res_named["gene1", "mean_expr"]),
        unname(expected_mean["gene1"]),
        tolerance = 1e-12
    )
    expect_equal(
        unname(res_named["gene2", "mean_expr"]),
        unname(expected_mean["gene2"]),
        tolerance = 1e-12
    )

    expect_equal(
        unname(res_named["gene1", "normalized_SEI"]),
        unname(expected_norm["gene1"]),
        tolerance = 1e-12
    )
    expect_equal(
        unname(res_named["gene2", "normalized_SEI"]),
        unname(expected_norm["gene2"]),
        tolerance = 1e-12
    )
})

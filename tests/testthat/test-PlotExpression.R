test_that("Test plotExpression", {
    df <- data.frame(
        x = c(rnorm(100, 1), rnorm(100, 5)),
        y = c(rnorm(100, 1), rnorm(100, 5)),
        cell = 1:200,
        cluster = rep(1:2, each = 100)
    )

    exp_mat <- data.frame(
        gene1 = c(runif(100, 4, 6), runif(100, 0, 2)),
        gene2 = c(runif(100, 0, 1), runif(100, 5, 8))
    )

    exp_mat <- t(exp_mat)

    colnames(exp_mat) <- df$cell

    expect_silent(plotExpression(
        data = df, exp_mat = exp_mat, shuffle = TRUE,
        genes = c("gene1", "gene2")
    ))

    expect_silent(plotExpression(
        data = df, exp_mat = exp_mat,
        sub_plot = TRUE, sub_cells = df$cell,
        genes = c("gene1")
    ))

    expect_error(plotExpression(
        data = df, exp_mat = "non_matrix",
        genes = c("gene1", "gene2")
    ))

    expect_error(plotExpression(
        data = df, exp_mat = exp_mat,
        genes = c("gene1", "gene2"),
        split_by = "not_in_columns"
    ))

    expect_error(plotExpression(
        data = df, exp_mat = exp_mat,
        genes = c("gene_not_exsit")
    ))

    expect_error(plotExpression(
        data = df[-c(1:10), ], exp_mat = exp_mat,
        genes = c("gene1", "gene2")
    ))

    p <- plotExpression(
        data = df,
        exp_mat = exp_mat,
        split_by = "cluster",
        sub_plot = TRUE,
        one_cluster = 1,
        genes = c("gene1", "gene2"),
        return_list = TRUE
    )

    expect_true(is.list(p))
})



test_that("Test plotSpatialExpression", {
    exp_mat <- matrix(runif(1000), nrow = 10)

    rownames(exp_mat) <- paste0("Gene", 1:10)

    spatial_distance <- runif(100)

    p <- plotSpatialExpression(
        exp_mat = exp_mat,
        spatial_distance = spatial_distance,
        genes = rownames(exp_mat)[1:5]
    )

    expect_true(inherits(p, "ggplot"))

    p1 <- plotSpatialExpression(
        exp_mat = exp_mat,
        spatial_distance = spatial_distance,
        genes = rownames(exp_mat)[1:5],
        scale_method = "minmax"
    )

    expect_true(inherits(p1, "ggplot"))

    expect_error(plotSpatialExpression(
        exp_mat = NULL,
        spatial_distance = spatial_distance,
        genes = rownames(exp_mat)[1:5]
    ))

    expect_error(plotSpatialExpression(
        exp_mat = exp_mat,
        row_gap = 2,
        spatial_distance = spatial_distance,
        genes = rownames(exp_mat)[1:5]
    ))

    expect_error(plotSpatialExpression(
        exp_mat = exp_mat,
        column_gap = 2,
        spatial_distance = spatial_distance,
        genes = rownames(exp_mat)[1:5]
    ))


    expect_error(plotSpatialExpression(
        exp_mat = exp_mat,
        spatial_distance = spatial_distance[1:3],
        genes = rownames(exp_mat)[1:5]
    ))

    expect_error(plotSpatialExpression(
        exp_mat = exp_mat,
        spatial_distance = spatial_distance,
        genes = c(rownames(exp_mat)[1:5], "unknown_gene")
    ))
})

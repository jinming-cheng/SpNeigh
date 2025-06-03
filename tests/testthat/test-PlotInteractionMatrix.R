test_that("Test PlotInteractionMatrix", {
    mat <- matrix(1:25, nrow = 5, ncol = 5)

    colnames(mat) <- paste0("C", 1:5)

    rownames(mat) <- paste0("C", 1:5)

    expect_silent(PlotInteractionMatrix(mat))

    expect_error(PlotInteractionMatrix(1:10))
})

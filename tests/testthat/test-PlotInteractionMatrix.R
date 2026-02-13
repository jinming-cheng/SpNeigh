test_that("Test plotInteractionMatrix", {
    mat <- matrix(1:25, nrow = 5, ncol = 5)

    colnames(mat) <- paste0("C", 1:5)

    rownames(mat) <- paste0("C", 1:5)

    expect_silent(plotInteractionMatrix(mat))

    expect_error(plotInteractionMatrix(1:10))
})

test_that("Test computeSpatialInteractionMatrix", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    interaction_matrix <- computeSpatialInteractionMatrix(coords)

    expect_true(is.matrix(interaction_matrix))
})


test_that("computeSpatialInteractionMatrix computes correct counts", {
    # Deterministic coordinates
    df <- data.frame(
        x = c(0, 1, 5, 6),
        y = c(0, 0, 0, 0),
        cell = c("A1", "A2", "B1", "B2"),
        cluster = c("A", "A", "B", "B")
    )

    res <- computeSpatialInteractionMatrix(df, k = 1)

    # Expected matrix
    expected <- matrix(
        c(
            2, 0,
            0, 2
        ),
        nrow = 2,
        byrow = TRUE
    )
    rownames(expected) <- c("A", "B")
    colnames(expected) <- c("A", "B")

    # Ensure same ordering
    res <- res[rownames(expected), colnames(expected)]

    expect_equal(res, expected)
})

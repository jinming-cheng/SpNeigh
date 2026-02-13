test_that("Test computeSpatialInteractionMatrix", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    interaction_matrix <- computeSpatialInteractionMatrix(coords)

    expect_true(is.matrix(interaction_matrix))
})

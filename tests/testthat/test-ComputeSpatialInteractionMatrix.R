test_that("Test ComputeSpatialInteractionMatrix", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    interaction_matrix <- ComputeSpatialInteractionMatrix(coords)

    expect_true(is.matrix(interaction_matrix))
})

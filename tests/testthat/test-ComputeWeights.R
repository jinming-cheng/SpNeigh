test_that("Test computeCentroidWeights", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    cells_c2 <- subset(coords, cluster == 2)$cell

    weights <- computeCentroidWeights(data = coords, cell_ids = cells_c2)

    expect_true(is.vector(weights))

    expect_true(is.vector(names(weights)))

    expect_error(computeCentroidWeights(
        data = coords[, c("x", "y")],
        cell_ids = cells_c2
    ))

    expect_error(computeCentroidWeights(
        data = coords,
        cell_ids = NULL
    ))
})

test_that("computeCentroidWeights computes correct inverse weights", {
    df <- data.frame(
        cell = c("c1", "c2", "c3", "c4"),
        x = c(0, 1, 0, 1),
        y = c(0, 0, 1, 1)
    )

    weights <- computeCentroidWeights(
        data = df,
        cell_ids = df$cell,
        scale = FALSE,
        method = "inverse"
    )

    # Expected distance
    d <- sqrt(0.5)

    expected <- rep(1 / (1 + d), 4)
    names(expected) <- df$cell

    expect_equal(weights, expected, tolerance = 1e-12)
})


test_that("computeCentroidWeights handles scaling correctly", {
    df <- data.frame(
        cell = c("c1", "c2", "c3", "c4"),
        x = c(0, 1, 0, 1),
        y = c(0, 0, 1, 1)
    )

    weights_scaled <- computeCentroidWeights(
        data = df,
        cell_ids = df$cell,
        scale = TRUE,
        method = "inverse"
    )

    expect_true(all(weights_scaled == 1))
})


test_that("computeCentroidWeights linear method works", {
    df <- data.frame(
        cell = c("c1", "c2", "c3", "c4"),
        x = c(0, 1, 0, 1),
        y = c(0, 0, 1, 1)
    )

    weights <- computeCentroidWeights(
        data = df,
        cell_ids = df$cell,
        scale = TRUE,
        method = "linear"
    )

    expected <- rep(1, 4)
    names(expected) <- df$cell

    expect_equal(weights, expected)
})

test_that("computeCentroidWeights gaussian method works", {
    df <- data.frame(
        cell = c("c1", "c2", "c3", "c4"),
        x = c(0, 1, 0, 1),
        y = c(0, 0, 1, 1)
    )

    sigma <- 0.5
    d <- sqrt(0.5)
    expected_val <- exp(-d^2 / (2 * sigma^2))

    weights <- computeCentroidWeights(
        data = df,
        cell_ids = df$cell,
        scale = FALSE,
        method = "gaussian",
        sigma = sigma
    )

    expected <- rep(expected_val, 4)
    names(expected) <- df$cell

    expect_equal(weights, expected, tolerance = 1e-12)
})


test_that("Test computeBoundaryWeights", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    boundary_points <- getBoundary(data = coords, one_cluster = 2)

    boundary_polys <- buildBoundaryPoly(boundary_points)

    cells_c2 <- subset(coords, cluster == 2)$cell

    weights <- computeBoundaryWeights(
        data = coords, cell_ids = cells_c2,
        boundary = boundary_points
    )

    expect_true(is.vector(weights))

    expect_true(is.vector(names(weights)))

    expect_error(computeBoundaryWeights(
        data = coords[, c("x", "y")], cell_ids = cells_c2,
        boundary = boundary_polys[1, ]
    ))

    expect_error(computeBoundaryWeights(
        data = coords, cell_ids = cells_c2,
        boundary = "not_a_boundary"
    ))

    expect_error(computeBoundaryWeights(
        data = coords, cell_ids = NULL,
        boundary = boundary_polys[1, ]
    ))


    boundary_edges <- splitBoundaryPolyByAnchor(boundary_polys[1, ])

    weights_edge <- computeBoundaryWeights(
        data = coords, cell_ids = cells_c2,
        boundary = boundary_edges[2, ]
    )

    expect_true(is.vector(weights))
})


test_that("computeBoundaryWeights works on simple geometry", {
    # ---- Simple square dataset ----
    df <- data.frame(
        cell = c("c1", "c2", "c3", "c4"),
        x = c(0, 1, 0, 1),
        y = c(0, 0, 1, 1)
    )

    # ---- Boundary: vertical line x = 0 ----
    boundary_line <- sf::st_sfc(
        sf::st_linestring(
            matrix(
                c(
                    0, 0,
                    0, 1
                ),
                ncol = 2, byrow = TRUE
            )
        )
    )
    boundary_sf <- sf::st_sf(geometry = boundary_line)

    # ---- Unscaled inverse weights ----
    w_unscaled <- computeBoundaryWeights(
        data = df,
        cell_ids = df$cell,
        boundary = boundary_sf,
        scale = FALSE,
        method = "inverse"
    )

    expected_dists <- c(0, 1, 0, 1)
    expected_weights <- 1 / (1 + expected_dists)

    expect_equal(
        unname(w_unscaled),
        expected_weights,
        tolerance = 1e-12
    )

    # ---- Scaled weights ----
    w_scaled <- computeBoundaryWeights(
        data = df,
        cell_ids = df$cell,
        boundary = boundary_sf,
        scale = TRUE,
        method = "inverse"
    )

    # After scaling:
    # dists become 0 and 1
    expected_scaled_weights <- 1 / (1 + c(0, 1, 0, 1))

    expect_equal(
        unname(w_scaled),
        expected_scaled_weights,
        tolerance = 1e-12
    )

    # ---- Structural checks ----
    expect_equal(length(w_scaled), 4)
    expect_equal(names(w_scaled), df$cell)
    expect_false(any(is.nan(w_scaled)))
    expect_true(all(w_scaled >= 0))
})

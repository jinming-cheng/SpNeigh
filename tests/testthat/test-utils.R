test_that("Test factorNaturalOrder", {
    x <- factorNaturalOrder(10:1)

    expect_true(is.factor(x))
})



test_that("Test theme_spneigh", {
    p <- ggplot2::ggplot() +
        theme_spneigh()

    expect_silent(p)
})

test_that("Test safeColorPalette", {
    n1 <- length(safeColorPalette(5))

    expect_equal(n1, 5)

    expect_message(length(safeColorPalette(20)))
})

test_that("Test FactorNaturalOrder", {
    x <- FactorNaturalOrder(10:1)

    expect_true(is.factor(x))
})



test_that("Test my_theme_ggplot", {
    p <- ggplot2::ggplot() +
        my_theme_ggplot()

    expect_silent(p)
})

test_that("Test SafeColorPalette", {
    n1 <- length(SafeColorPalette(5))

    expect_equal(n1, 5)

    expect_message(length(SafeColorPalette(20)))
})

test_that("Test factorNaturalOrder for numeric input", {
    x <- 10:1
    f <- factorNaturalOrder(x)

    expect_s3_class(f, "factor")

    # Levels should be sorted in increasing order
    expect_equal(levels(f), as.character(1:10))

    # Factor values should preserve original order
    expect_equal(as.character(f), as.character(x))
})

test_that("Test factorNaturalOrder for charater input", {
    x <- c("a3", "a1", "a2", "a10")
    f <- factorNaturalOrder(x)

    expect_equal(levels(f), c("a1", "a2", "a3", "a10"))
})



test_that("Test theme_spneigh", {
    p <- ggplot2::ggplot() +
        theme_spneigh()

    expect_true(inherits(theme_spneigh(), "theme"))

    expect_s3_class(p, "ggplot")
})

test_that("Test safeColorPalette", {
    n1 <- length(safeColorPalette(5))

    expect_equal(n1, 5)

    expect_message(length(safeColorPalette(20)))

    n2 <- length(safeColorPalette(20, verbose = FALSE))
    expect_equal(n2, 20)
})

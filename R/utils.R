#' Create a factor with natural (human-friendly) ordering
#'
#' Converts a character or numeric vector into a factor where the levels are
#' ordered naturally (e.g., `a1`, `a2`, ..., `a10` instead of
#' lexicographically as `a1`, `a10`, `a2`, ...).
#' This is useful for plotting or labeling grouped data where numeric
#' substrings should follow numeric order.
#'
#' @param x A character or numeric vector to convert to a factor with
#'          natural order.
#'
#' @return A factor with levels sorted in natural (human-readable) order.
#'
#' @export
#'
#' @examples
#' # Numeric vector
#' factorNaturalOrder(10:1)
#'
#' # Character vector with embedded numbers
#' factorNaturalOrder(c("a11", "a12", "a1", "a2", "a"))
#'
factorNaturalOrder <- function(x) {
    x <- as.character(x)
    ord <- stringr::str_order(x, numeric = TRUE)
    factor(x, levels = unique(x[ord]))
}



#' Custom ggplot2 theme
#'
#' Provides a consistent and clean visual style for plots generated
#' within this package. This theme builds on `theme_classic()`
#' and adjusts text sizes and margins for better readability in figures.
#'
#' @return A `ggplot2` theme object that can be added to ggplot visualizations.
#'
#' @export
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(mpg, wt)) +
#'     geom_point() +
#'     theme_spneigh()
theme_spneigh <- function() {
    ggplot2::theme_classic() +
        ggplot2::theme(
            axis.text = ggplot2::element_text(size = 12),
            axis.title = ggplot2::element_text(size = 15, face = "bold"),
            strip.text = ggplot2::element_text(size = 12, face = "bold"),
            plot.title = ggplot2::element_text(
                size = 15, face = "bold",
                hjust = 0.5
            ),
            plot.margin = ggplot2::margin(0.5, 0.5, 0.5, 0.5, "cm")
        )
}


#' Generate a safe color palette for discrete clusters
#'
#' Returns a vector of visually distinct colors for plotting discrete clusters.
#' Uses `colors15_cheng` as the base palette. If the number of clusters exceeds
#' the base palette, additional colors are generated using `scales::hue_pal()`.
#'
#' @param n_clusters Number of unique clusters.
#' @param base_colors A character vector of base colors.
#'                    Default is `colors15_cheng`.
#' @param verbose Logical. If `TRUE`, displays a message when
#'                extended colors are used. Default is `TRUE`.
#'
#' @return A character vector of colors of length `n_clusters`.
#' @export
#'
#' @examples
#' safeColorPalette(5)
#' safeColorPalette(20)
#'
safeColorPalette <- function(
    n_clusters = NULL,
    base_colors = colors15_cheng,
    verbose = TRUE) {
    n_base <- length(base_colors)

    if (n_clusters <= n_base) {
        return(base_colors[seq_len(n_clusters)])
    }

    # Generate extended palette using scales::hue_pal()
    extra_needed <- n_clusters - n_base
    extended_colors <- c(
        base_colors,
        scales::hue_pal()(extra_needed)
    )

    if (verbose) {
        msg_fmt <- paste0(
            "Note: %d base colors provided. ",
            "Generated %d additional colors using `scales::hue_pal()` ",
            "(total = %d) to match %d clusters."
        )

        message(sprintf(msg_fmt, n_base, extra_needed, n_clusters, n_clusters))
    }

    return(extended_colors)
}

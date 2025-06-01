#' Plot spatial weights for cells on a spatial plot
#'
#' Visualizes spatial weights by plotting cell coordinates colored by
#' a numeric weight value.
#' Only cells in `data` that match the names in `weights` will be included.
#' This is useful for visualizing spatial trends such as proximity to
#' boundaries or centroids.
#'
#' @inheritParams ComputeCentroidWeights
#' @importFrom rlang .data
#' @importFrom methods is
#' @param weights A named numeric vector of spatial weights,
#'                with cell IDs as names.
#' @param point_size Numeric. Point size of the cells in the plot.
#'                   Default is `0.2`.
#' @param theme_ggplot A ggplot2 theme object. Default is `my_theme_ggplot()`.
#'
#' @return A `ggplot2` object displaying a scatter plot of cells
#'         colored by weights.
#' @export
#' @examples
#' # Load coordinate data
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'     package = "SpNeigh"
#' ))
#'
#' # Cells in cluster 2
#' cells_c2 <- subset(coords, cluster == 2)[, "cell"]
#'
#' # Compute centroid weights and plot
#' weights_cen <- ComputeCentroidWeights(data = coords, cell_ids = cells_c2)
#' PlotWeights(data = coords, weights = weights_cen)
#'
#' # Compute boundary weights and plot
#' boundary_points <- GetBoundary(data = coords, one_cluster = 2)
#' weights_bon <- ComputeBoundaryWeights(
#'     data = coords, cell_ids = cells_c2,
#'     boundary = boundary_points
#' )
#' PlotWeights(data = coords, weights = weights_bon)
#'
PlotWeights <- function(
    data = NULL,
    weights = NULL,
    point_size = 0.2,
    theme_ggplot = my_theme_ggplot()) {
    # Extract spatial coordinates
    if (is(data, "Seurat")) {
        sp_coords <- ExtractCoords(data, extract_cluster = FALSE)
    } else {
        sp_coords <- data
    }

    if (!all(c("cell", "x", "y") %in% colnames(sp_coords))) {
        stop("`data` must contain columns: 'cell', 'x', and 'y'")
    }

    if (is.null(names(weights))) {
        stop("`weights` must be a named numeric vector with cell IDs as names")
    }

    # Match weights to spatial coordinates
    matched_idx <- match(sp_coords$cell, names(weights))
    matched_weights <- weights[matched_idx]

    # Filter out unmatched cells
    keep <- !is.na(matched_weights)
    plot_df <- sp_coords[keep, ]
    plot_df$weights <- matched_weights[keep]

    if (nrow(plot_df) == 0) {
        stop("No matching cell IDs between `data` and `weights`")
    }

    # Generate plot
    ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = .data$x, y = .data$y, color = weights)
    ) +
        ggplot2::geom_point(size = point_size) +
        ggplot2::scale_color_viridis_c(option = "plasma", name = "Weight") +
        theme_ggplot
}

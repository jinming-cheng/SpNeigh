#' Plot cells located within spatial boundaries or ring regions
#'
#' Visualizes cells that fall within defined spatial regions
#' (boundaries or rings), typically obtained using
#' the `getCellsInside()` function. The cells are colored by cluster,
#' and the function offers two plotting modes: using `geom_sf()`
#' (with fixed 1:1 aspect ratio) or `geom_point()` (with more flexible layout).
#'
#' @inheritParams plotBoundary
#' @importFrom rlang .data
#' @param cells_inside An `sf` object of cells returned by `getCellsInside()`.
#'                     Must contain `cluster` and `region_id` columns.
#' @param fixed_aspect_ratio Logical. If `TRUE`, uses `geom_sf()` to preserve
#'                           spatial scale. If `FALSE`, uses `geom_point()`
#'                           with extracted coordinates. Default is `TRUE`.
#'
#' @return A `ggplot` object showing the cells colored by cluster within
#'         spatial regions.
#'
#' @export
#' @examples
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'     package = "SpNeigh"
#' ))
#'
#' # Plot cells inside boundaries
#' boundary_points <- getBoundary(data = coords, one_cluster = 2)
#' cells_inside <- getCellsInside(data = coords, boundary = boundary_points)
#' plotCellsInside(cells_inside)
#'
#' # Plot cells inside rings
#' ring_regions <- getRingRegion(boundary = boundary_points, dist = 100)
#' cells_ring <- getCellsInside(data = coords, boundary = ring_regions)
#' plotCellsInside(cells_ring, fixed_aspect_ratio = FALSE)
#'
plotCellsInside <- function(
    cells_inside = NULL,
    point_size = 0.5,
    colors = colors15_cheng,
    theme_ggplot = theme_spneigh(),
    legend_size = 2,
    fixed_aspect_ratio = TRUE) {
    # Ensure cluster is a factor with ordered levels
    if (is.null(levels(cells_inside$cluster))) {
        cells_inside$cluster <- factorNaturalOrder(cells_inside$cluster)
    }

    # Assign colors to clusters
    n_clusters <- nlevels(cells_inside$cluster)
    named_colors <- safeColorPalette(
        n_clusters = n_clusters,
        base_colors = colors
    )
    names(named_colors) <- levels(cells_inside$cluster)

    # Extract coords if using geom_point
    coords <- sf::st_coordinates(cells_inside)
    colnames(coords) <- c("x", "y")
    df <- sf::st_drop_geometry(cells_inside)
    df <- cbind(df, coords)
    df$cluster <- factor(df$cluster, levels = levels(cells_inside$cluster))

    # Choose plotting method
    if (fixed_aspect_ratio) {
        p <- ggplot2::ggplot() +
            ggplot2::geom_sf(
                data = cells_inside,
                size = point_size,
                ggplot2::aes(color = .data$cluster)
            )
    } else {
        p <- ggplot2::ggplot() +
            ggplot2::geom_point(
                data = df,
                size = point_size,
                ggplot2::aes(x = .data$x, y = .data$y, color = .data$cluster)
            )
    }

    # Add theme and colors
    p <- p + theme_ggplot +
        ggplot2::scale_color_manual(values = named_colors) +
        ggplot2::guides(
            colour = ggplot2::guide_legend(
                override.aes = list(size = legend_size)
            )
        )

    return(p)
}

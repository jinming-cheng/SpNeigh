#' Plot spatial cell coordinates with cluster boundaries
#'
#' Plots spatial cell locations and overlays cluster or population boundaries
#' if available.
#' If `boundary` is not provided and `one_cluster` is specified, the boundary
#' will be automatically generated using the `getBoundary()` function.
#' This function supports plotting either all clusters
#' or a specific cluster using `sub_plot = TRUE`. Boundaries can be overlaid
#' as polygons to visualize spatial subregions.
#'
#' @inheritParams getBoundary
#' @inheritParams extractCoords
#' @importFrom rlang .data
#' @param one_cluster The cluster ID to plot and optionally compute its
#'                    boundary. Required if `sub_plot = TRUE` and `boundary`
#'                    is not provided.
#' @param boundary A data frame with columns `x`, `y`, and `region_id`
#'                 or an `sf` object of `POLYGON` or `LINESTRING` geometries.
#' @param colors A vector of cluster colors. Default uses `colors15_cheng`.
#' @param point_size Numeric. Size of the points representing cells.
#'                   Default is 0.5.
#' @param color_boundary Color for boundary lines. Default is `"black"`.
#' @param linewidth_boundary Numeric. Line width for boundary outlines.
#'                           Default is 1.5.
#' @param sub_plot Logical. If `TRUE`, only cells from the specified
#'                 `one_cluster` are plotted.
#'                 If `FALSE` (default), all clusters are plotted.
#' @param split_by Optional column name in `data` to facet the plot by
#'                 (e.g., sample, condition).
#' @param ncol Number of columns in the faceted plot when `split_by` is used.
#'             Passed to `ggplot2::facet_wrap()`. Default is `NULL`,
#'             which lets ggplot2 determine layout automatically.
#' @param angle_x_label Numeric angle (in degrees) to rotate the x-axis labels.
#'                      Useful for improving label readability in faceted
#'                      or dense plots. Default is 0 (no rotation).
#' @param theme_ggplot A ggplot2 theme object. Default is `theme_spneigh()`.
#' @param legend_size Numeric. Size of legend keys. Default is 2.
#' @param ... Additional arguments passed to \code{\link{getBoundary}}
#'            when auto-generating boundaries.
#'
#' @return A ggplot object displaying the spatial layout of cells,
#'         with optional boundary overlays.
#' @export
#' @examples
#' # Load coordinates
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'     package = "SpNeigh"
#' ))
#' head(coords)
#'
#' # Plot all cells without boundaries
#' plotBoundary(coords)
#'
#' # Plot one cluster and its boundary
#' plotBoundary(coords, one_cluster = 2)
#'
#' # Manually compute boundary and plot
#' boundary_points <- getBoundary(
#'     data = coords, one_cluster = 2,
#'     eps = 120, minPts = 10
#' )
#' plotBoundary(data = coords, boundary = boundary_points)
#'
#' # plotBoundary for a SpatialExperiment object
#' logNorm_expr <- readRDS(system.file("extdata", "LogNormExpr.rds",
#'     package = "SpNeigh"
#' ))
#' coords_sub <- subset(coords, cluster %in% c("0", "2"))
#' coords_sub <- as.matrix(coords_sub[, c("x", "y")])
#' metadata_sub <- subset(
#'     coords[, c("cell", "cluster")],
#'     cluster %in% c("0", "2")
#' )
#'
#' spe <- SpatialExperiment::SpatialExperiment(
#'     assay = list("logcounts" = logNorm_expr),
#'     colData = metadata_sub,
#'     spatialCoords = coords_sub
#' )
#'
#' plotBoundary(data = spe, one_cluster = 2)
#'
#' # plotBoundary for a Seurat object
#' seu_sp <- Seurat::CreateSeuratObject(
#'     assay = "Spatial",
#'     counts = logNorm_expr,
#'     meta.data = metadata_sub
#' )
#' SeuratObject::LayerData(seu_sp,
#'     assay = "Spatial",
#'     layer = "data"
#' ) <- logNorm_expr
#'
#' cents <- SeuratObject::CreateCentroids(coords_sub[, c("x", "y")])
#' fov <- SeuratObject::CreateFOV(
#'     coords = list("centroids" = cents),
#'     type = c("centroids"),
#'     assay = "Spatial"
#' )
#' seu_sp[["fov"]] <- fov
#'
#' seu_sp$seurat_clusters <- seu_sp$cluster
#'
#' plotBoundary(data = seu_sp, one_cluster = 2, eps = 120)
#'
plotBoundary <- function(
    data = NULL,
    cluster_col = NULL,
    one_cluster = NULL,
    boundary = NULL,
    colors = colors15_cheng,
    point_size = 0.5,
    color_boundary = "black",
    linewidth_boundary = 1.5,
    sub_plot = FALSE,
    split_by = NULL,
    ncol = NULL,
    angle_x_label = 0,
    theme_ggplot = theme_spneigh(),
    legend_size = 2,
    ...) {
    # Extract coordinates from data
    sp_coords <- extractCoords(data = data, cluster_col = cluster_col)

    # Ensure cluster is a factor with ordered levels if not provided
    if (is.null(levels(sp_coords$cluster))) {
        sp_coords$cluster <- factorNaturalOrder(sp_coords$cluster)
    }

    # Auto-generate boundary if needed
    if (is.null(boundary) && !is.null(one_cluster)) {
        boundary <- getBoundary(
            data = sp_coords,
            one_cluster = one_cluster, ...
        )
    }

    # Choose plotting data: full or single cluster
    if (sub_plot) {
        if (is.null(one_cluster)) {
            stop("`one_cluster` must be specified when `sub_plot = TRUE`.")
        }
        data_for_plot <- dplyr::filter(sp_coords, .data$cluster == one_cluster)
    } else {
        data_for_plot <- sp_coords
    }

    # Assign colors to clusters
    n_clusters <- nlevels(sp_coords$cluster)
    named_colors <- safeColorPalette(
        n_clusters = n_clusters,
        base_colors = colors
    )
    names(named_colors) <- levels(sp_coords$cluster)

    # Base plot
    p <- ggplot2::ggplot(
        data_for_plot, ggplot2::aes(x = .data$x, y = .data$y)
    ) +
        ggplot2::geom_point(
            ggplot2::aes(color = .data$cluster),
            size = point_size
        ) +
        theme_ggplot +
        ggplot2::scale_color_manual(values = named_colors) +
        ggplot2::guides(colour = ggplot2::guide_legend(
            override.aes = list(size = legend_size)
        )) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(
            angle = angle_x_label, hjust = 1
        ))

    # Add boundary if provided or generated
    if (!is.null(boundary)) {
        p <- p + addBoundary(
            boundary = boundary,
            color_boundary = color_boundary,
            linewidth_boundary = linewidth_boundary
        )
    }

    # Make faceted plots if split_by is provided
    if (!is.null(split_by) && split_by %in% colnames(data_for_plot)) {
        p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[[split_by]]),
            ncol = ncol
        )
    }

    return(p)
}

#' Add boundary outlines to a spatial ggplot
#'
#' Overlays spatial boundaries (from points or polygons) onto
#' a `ggplot2` object.
#' The input can be either a data frame of boundary points or an `sf` object
#' of `POLYGON` geometries.
#' This function is typically used as an additive layer (`+ addBoundary(...)`)
#' in conjunction with a base plot created using `plotBoundary()`.
#'
#' @inheritParams plotBoundary
#' @importFrom rlang .data
#'
#' @return A `ggplot2::geom_path` layer that can be added to an existing plot.
#'
#' @export
#' @examples
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'     package = "SpNeigh"
#' ))
#'
#' # Automatically get boundary for a cluster
#' boundary_points <- getBoundary(
#'     data = coords, one_cluster = 2,
#'     subregion_method = "dbscan",
#'     eps = 120, minPts = 10
#' )
#'
#' # Plot with boundary overlay
#' plotBoundary(coords) + addBoundary(boundary_points)
#'
addBoundary <- function(
    boundary = NULL,
    color_boundary = "black",
    linewidth_boundary = 1.5) {
    # Check boundary type
    if (!inherits(boundary, c("data.frame", "sf"))) {
        stop("`boundary` should be a data.frame or an sf object.")
    }

    # Convert polygons to boundary points if needed
    if (inherits(boundary, "sf")) {
        boundary <- boundaryPolyToPoints(boundary)
    }

    # Ensure required columns
    if (!all(c("x", "y", "region_id") %in% colnames(boundary))) {
        stop("`boundary` must contain columns: 'x', 'y', and 'region_id'.")
    }

    ggplot2::geom_path(
        data = boundary,
        ggplot2::aes(x = .data$x, y = .data$y, group = .data$region_id),
        color = color_boundary,
        linewidth = linewidth_boundary
    )
}


#' Add boundary polygons or linestrings to a spatial plot
#'
#' Overlays boundary polygons or linestrings
#' on a spatial `ggplot2` plot. This function adds an `sf` geometry layer
#' to display complete boundary shapes for visualizing spatial clusters,
#' rings, or enriched zones.
#'
#' @inheritParams plotBoundary
#' @importFrom rlang .data
#' @param boundary_poly An `sf` object containing `POLYGON` or `LINESTRING`
#'                      geometries and a `region_id` column.
#' @return A `ggplot2::geom_sf` layer that can be added to an existing plot.
#'
#' @export
#' @examples
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'     package = "SpNeigh"
#' ))
#'
#' boundary_points <- getBoundary(
#'     data = coords, one_cluster = 2,
#'     subregion_method = "dbscan",
#'     eps = 120, minPts = 10
#' )
#' boundary_polys <- buildBoundaryPoly(boundary_points)
#'
#' # Add boundary polygons to the plot
#' plotBoundary(coords) +
#'     addBoundaryPoly(boundary_polys, color_boundary = "blue")
#'
addBoundaryPoly <- function(
    boundary_poly,
    color_boundary = "black",
    linewidth_boundary = 1.5) {
    # Validate input
    if (!inherits(boundary_poly, "sf")) {
        stop("`boundary_poly` must be an sf object.")
    }

    ggplot2::geom_sf(
        data = boundary_poly,
        inherit.aes = FALSE,
        fill = NA,
        color = color_boundary,
        linewidth = linewidth_boundary
    )
}


#' Plot filled spatial regions inside boundaries or rings
#'
#' Creates a ggplot of spatial regions (e.g., subregions or ring areas)
#' using filled polygons.
#' Each region is automatically assigned a fill color based on its `region_id`.
#' This function is commonly used to visualize the area inside
#' spatial boundaries or surrounding rings, such as those created by
#' `buildBoundaryPoly()` or `getRingRegion()`.
#'
#' @importFrom rlang .data
#' @inheritParams plotBoundary
#' @param boundary_poly An `sf` object of `POLYGON` geometries
#'                      containing a `region_id` column. Typically created by
#'                      `buildBoundaryPoly()` or `getRingRegion()`.
#' @param alpha Numeric value controlling the transparency of the filled
#'              regions. Default is `0.5`.
#' @param color_boundary Color of the region outlines. Default is `"black"`.
#' @param linewidth_boundary Numeric line width of the region borders.
#'                           Default is `1`.
#' @param ... Additional arguments passed to `ggplot2::geom_sf()`.
#'
#' @return A `ggplot` object displaying filled spatial regions by `region_id`.
#' @export
#' @examples
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'     package = "SpNeigh"
#' ))
#'
#' # Plot filled boundary regions
#' boundary_points <- getBoundary(
#'     data = coords, one_cluster = 2,
#'     eps = 120, minPts = 10
#' )
#' boundary_polys <- buildBoundaryPoly(boundary_points)
#' plotRegion(boundary_poly = boundary_polys)
#'
#' # Plot filled ring regions
#' ring_regions <- getRingRegion(boundary = boundary_points, dist = 100)
#' plotRegion(boundary_poly = ring_regions)
#'
plotRegion <- function(
    boundary_poly = NULL,
    alpha = 0.5,
    color_boundary = "black",
    linewidth_boundary = 1,
    theme_ggplot = theme_spneigh(),
    ...) {
    # Check input
    if (!inherits(boundary_poly, c("sf"))) {
        stop("`boundary_poly` should be an sf object.")
    }

    # Ensure region_id is a factor for consistent fill coloring
    if (is.null(levels(boundary_poly$region_id))) {
        boundary_poly$region_id <- factorNaturalOrder(boundary_poly$region_id)
    }

    ggplot2::ggplot() +
        ggplot2::geom_sf(
            data = boundary_poly,
            ggplot2::aes(fill = .data$region_id),
            alpha = alpha,
            color = color_boundary,
            linewidth = linewidth_boundary, ...
        ) +
        ggplot2::labs(fill = "Region ID") +
        theme_ggplot
}


#' Plot boundary edges or segmented boundary lines
#'
#' Plots boundary edges as colored `LINESTRING` or `POLYGON` outlines
#' using `ggplot2`.
#' This function is especially useful for visualizing specific edges
#' extracted from a polygon using `splitBoundaryPolyByAnchor()`.
#'
#' @importFrom rlang .data
#' @inheritParams plotBoundary
#' @inheritParams addBoundaryPoly
#' @param linewidth_boundary Numeric value specifying the line width of
#'                           the edge. Default is `1`.
#' @param ... Additional arguments passed to `ggplot2::geom_sf()`.
#'
#' @return A `ggplot` object displaying the edge outlines colored
#'         by `region_id`.
#' @export
#' @examples
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'     package = "SpNeigh"
#' ))
#'
#' # Build boundary polygon and plot its outline
#' boundary_points <- getBoundary(
#'     data = coords, one_cluster = 2,
#'     eps = 120, minPts = 10
#' )
#' boundary_polys <- buildBoundaryPoly(boundary_points)
#' plotEdge(boundary_poly = boundary_polys)
#'
#' # Split a polygon into edge segments and plot
#' boundary_edges <- splitBoundaryPolyByAnchor(boundary_polys[1, ])
#' plotEdge(boundary_poly = boundary_edges)
#'
plotEdge <- function(
    boundary_poly = NULL,
    linewidth_boundary = 1,
    theme_ggplot = theme_spneigh(),
    ...) {
    # Check input
    if (!inherits(boundary_poly, c("sf"))) {
        stop("`boundary_poly` should be an sf object.")
    }

    # Ensure region_id is a factor for grouped coloring
    if (is.null(levels(boundary_poly$region_id))) {
        boundary_poly$region_id <- factorNaturalOrder(boundary_poly$region_id)
    }

    ggplot2::ggplot() +
        ggplot2::geom_sf(
            data = boundary_poly,
            ggplot2::aes(color = .data$region_id),
            fill = NA,
            linewidth = linewidth_boundary, ...
        ) +
        ggplot2::labs(color = "Region ID") +
        theme_ggplot
}

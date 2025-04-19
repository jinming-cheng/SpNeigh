
#' Bar plot of cluster statistics for cells inside boundaries or ring regions
#'
#' Creates a bar plot to visualize the distribution of cells inside spatial regions (e.g., boundaries or rings),
#' either as raw counts or proportions per cluster. The plot is faceted by `region_id` to show statistics across
#' multiple spatial subregions.
#'
#' @inheritParams PlotBoundary
#' @importFrom rlang .data
#' @param cell_stats A data frame containing summarized cell statistics, typically the output from `StatsCellsInside()`.
#'                   Must include columns `region_id`, `cluster`, and the specified `stat_column`.
#' @param stat_column Character. Column name in `cell_stats` to use for the y-axis.
#'        Options are `"count"` (number of cells) or `"proportion"` (relative fraction per region).
#'
#' @return A `ggplot2` object showing a faceted bar plot of cell statistics per region.
#'
#' @export
#' @examples
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'                               package = "SpNeigh"))
#'
#' boundary_points <- GetBoundary(data = coords, one_cluster = 2,
#'                                eps = 120, minPts = 10)
#' cells_inside <- GetCellsInside(data = coords, boundary = boundary_points)
#' stats_cells <- StatsCellsInside(cells_inside)
#'
#' PlotStatsBar(stats_cells, stat_column = "proportion")
#' PlotStatsBar(stats_cells, stat_column = "count")
#'
PlotStatsBar <- function(cell_stats = NULL,
                         stat_column = c("proportion", "count"),
                         colors = my_colors_15,
                         angle_x_label = 0,
                         theme_ggplot = my_theme_ggplot()) {

  stat_column <- match.arg(stat_column)

  # Ensure cluster is a factor with defined levels
  if (is.null(levels(cell_stats$cluster))) {
    cell_stats$cluster <- FactorNaturalOrder(cell_stats$cluster)
  }

  # Assign cluster colors
  named_colors <- colors[1:nlevels(cell_stats$cluster)]
  names(named_colors) <- levels(cell_stats$cluster)

  # Create plot
  p <- ggplot2::ggplot(cell_stats,
                       ggplot2::aes(x = .data$cluster,
                                    y = .data[[stat_column]],
                                    fill = .data$cluster)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_wrap(~.data$region_id) +
    ggplot2::scale_fill_manual(values = named_colors) +
    theme_ggplot +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = angle_x_label, hjust = 1))

  # If proportions, format y-axis as percentage
  if (stat_column == "proportion") {
    p <- p + ggplot2::scale_y_continuous(labels = scales::percent_format())
  }

  return(p)
}



#' Pie chart or donut chart of cluster proportions inside spatial regions
#'
#' Generates pie or donut charts to visualize the proportion of cells from different clusters
#' within each spatial region (e.g., boundary or ring). The plot is faceted by `region_id`
#' to show the composition of each spatial subregion. Optionally, percentage labels can be added
#' with filtering based on a minimum proportion threshold.
#'
#' @inheritParams PlotBoundary
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @param cell_stats A data frame containing cluster statistics per region, typically
#'        the output from `StatsCellsInside()`. Must include columns `region_id`, `cluster`, and `proportion`.
#' @param plot_donut Logical. If `TRUE`, a donut chart is generated; otherwise, a pie chart. Default is `FALSE`.
#' @param add_labels Logical. If `TRUE`, percentage labels are displayed inside each slice. Default is `TRUE`.
#' @param label_cutoff Numeric. Proportional threshold below which labels are hidden (e.g., `0.01` = 1%). Default is `0.01`.
#' @param label_color Character. Color of the percentage labels. Default is `"white"`.
#' @param label_size Numeric. Text size for percentage labels. Default is `4`.
#' @param label_nudge_x Numeric. Horizontal adjustment for label positioning. Default is `0.1`.
#'
#' @return A `ggplot2` object representing a faceted pie or donut chart per spatial region.
#'
#' @export
#' @examples
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'                               package = "SpNeigh"))
#'
#' boundary_points <- GetBoundary(data = coords, one_cluster = 2,
#'                                eps = 120, minPts = 10)
#' cells_inside <- GetCellsInside(data = coords, boundary = boundary_points)
#' stats_cells <- StatsCellsInside(cells_inside)
#'
#' PlotStatsPie(stats_cells, add_labels = FALSE)
#' PlotStatsPie(stats_cells, label_cutoff = 0)
#' PlotStatsPie(stats_cells, label_cutoff = 0.01)
#' PlotStatsPie(stats_cells, plot_donut = TRUE)
#'
PlotStatsPie <- function(cell_stats = NULL,
                         plot_donut = FALSE,
                         add_labels = TRUE,
                         label_cutoff = 0.01,
                         label_color = "white",
                         label_size = 4,
                         label_nudge_x = 0.1,
                         colors = my_colors_15) {


  # Ensure cluster is a factor with defined levels
  if (is.null(levels(cell_stats$cluster))) {
    cell_stats$cluster <- FactorNaturalOrder(cell_stats$cluster)
  }

  # Assign cluster colors
  named_colors <- colors[1:nlevels(cell_stats$cluster)]
  names(named_colors) <- levels(cell_stats$cluster)

  # Calculate label positions and cumulative slices
  cell_stats <- cell_stats %>%
    dplyr::group_by(.data$region_id) %>%
    dplyr::arrange(.data$region_id, dplyr::desc(.data$cluster)) %>%
    dplyr::mutate(
      cumulative = cumsum(.data$proportion),
      midpoint = .data$cumulative - .data$proportion / 2,
      label_text = ifelse(.data$proportion < label_cutoff, "",
                          paste0(round(.data$proportion * 100), "%"))
    )

  # x-axis dummy for donut chart
  cell_stats$dummy <- if (plot_donut) 2 else 1

  # Base plot
  p <- ggplot2::ggplot(cell_stats,
                       ggplot2::aes(x = .data$dummy,
                                    y = .data$proportion,
                                    fill = .data$cluster)) +
    ggplot2::geom_col(width = 1, color = "white") +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::facet_wrap(~.data$region_id) +
    ggplot2::scale_fill_manual(values = named_colors) +
    ggplot2::theme_void()

  # Optional percentage labels
  if (add_labels) {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(y = .data$midpoint, label = .data$label_text),
        color = label_color,
        size = label_size,
        position = ggplot2::position_nudge(x = label_nudge_x)
      )
  }

  # Add hole if donut
  if (plot_donut) {
    p <- p + ggplot2::xlim(0.5, 2.5)
  }

  return(p)
}

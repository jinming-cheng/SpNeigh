
#' Bar plot of statistics of cells inside boundaries or rings
#'
#' @inheritParams PlotBoundary
#' @importFrom rlang .data
#' @param cell_stats Statistics of cells returned by `StatsCellsInside()`.
#' @param stat_column A column name in `cell_stats` and used for y axis.
#' @param angle_x_label Angle to rotate the x-axis labels.
#' @export
#' @examples
#' # Load coordinates
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'                               package = "SpNeigh"))
#' head(coords)
#'
#' # Obtain statistics of cells inside boundaries
#' boundary_points <- GetBoundary(data = coords, one_cluster = 2,
#'                                eps = 120, minPts = 10)
#' cells_inside <- GetCellsInside(data = coords, boundary =  boundary_points)
#' stats_cells = StatsCellsInside(cells_inside)
#'
#' # Plot proportion of cells in different clusters using bar plot
#' PlotStatsBar(stats_cells, stat_column = "proportion")
#'
#' # Plot count of cells in different clusters using bar plot
#' PlotStatsBar(stats_cells, stat_column = "count")
#'

PlotStatsBar <- function(cell_stats = NULL,
                         stat_column = c("proportion","count"),
                         colors = my_colors_15,
                         angle_x_label = 0,
                         theme_ggplot = my_theme_ggplot()){

  stat_column <- match.arg(stat_column)

  # Add names to colors
  named_colors <- colors[1:nlevels(cell_stats$cluster)]
  names(named_colors) <- levels(cell_stats$cluster)

  p <- ggplot2::ggplot(data = cell_stats,
                       ggplot2::aes(x = .data$cluster, y = .data[[stat_column]],
                                    fill = .data$cluster)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_wrap(~ .data$region_id) +
    ggplot2::scale_fill_manual(values = named_colors) +
    theme_ggplot +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = angle_x_label, hjust = 1))

  if(stat_column=="proportion"){
    p <- p + ggplot2::scale_y_continuous(labels = scales::percent_format()) }

  p
}

#' Pie chart or donut chart of statistics of cells inside boundaries or rings
#'
#' @inheritParams PlotBoundary
#' @importFrom rlang .data
#' @param cell_stats Statistics of cells returned by `StatsCellsInside()`.
#' @param plot_donut A logical value indicates whether to make a donut chart. Default is FALSE.
#' @param add_labels A logical value indicates whether to add percentage labels inside the plot. Default is TRUE.
#' @param label_cutoff Labels of proportion smaller than this cutoff will not be plotted. Default is 0.01 (or 1%). Set it to 0 to plot labels of all percentages.
#' @param label_color Color of the percentage labels.
#' @param label_size Text size of the percentage labels.
#' @param label_nudge_x Nudge x coordinate to adjust positions for the percentage labels.
#' @export
#' @examples
#' # Load coordinates
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'                               package = "SpNeigh"))
#' head(coords)
#'
#' # Obtain statistics of cells inside boundaries
#' boundary_points <- GetBoundary(data = coords, one_cluster = 2,
#'                                eps = 120, minPts = 10)
#' cells_inside <- GetCellsInside(data = coords, boundary =  boundary_points)
#' stats_cells = StatsCellsInside(cells_inside)
#'
#' # Plot proportion using pie chart without adding percentage labels
#' PlotStatsPie(stats_cells, add_labels = FALSE)
#'
#' # Plot proportion using pie chart showing labels of all percentages
#' PlotStatsPie(stats_cells, label_cutoff = 0)
#'
#' # Plot proportion using pie chart showing labels of percentages greater than 1%
#' PlotStatsPie(stats_cells, label_cutoff = 0.01)
#'
#' # Plot proportion using donut chart
#' PlotStatsPie(stats_cells, plot_donut = TRUE)
#'

PlotStatsPie <- function(cell_stats = NULL,
                         plot_donut = FALSE,
                         add_labels = TRUE,
                         label_cutoff = 0.01,
                         label_color = "white",
                         label_size = 4,
                         label_nudge_x = 0.1,
                         colors = my_colors_15){

  # Add names to colors
  named_colors <- colors[1:nlevels(cell_stats$cluster)]
  names(named_colors) <- levels(cell_stats$cluster)

  cell_stats <- cell_stats %>%
    dplyr::group_by(.data$region_id) %>%
    dplyr::arrange(.data$region_id, dplyr::desc(.data$cluster)) %>%
    dplyr::mutate(cumulative = cumsum(.data$proportion),
                  midpoint = .data$cumulative - .data$proportion / 2,
                  label_text = ifelse(.data$proportion < label_cutoff, "",
                                      paste0(round(.data$proportion * 100), "%"))
    )

  if(plot_donut){
    cell_stats$dummy <-2
  }else{
    cell_stats$dummy <- 1
  }

  p <- ggplot2::ggplot(data = cell_stats,
                       ggplot2::aes(x = .data$dummy, y = .data$proportion,
                                    fill = .data$cluster)) +
    ggplot2::geom_col(width = 1, color = "white") +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::facet_wrap(~ .data$region_id) +
    ggplot2::scale_fill_manual(values = named_colors) +
    ggplot2::theme_void()

  if(add_labels){
    p <- p +
      ggplot2::geom_text(ggplot2::aes(y = .data$midpoint, label = .data$label_text),
                         color = label_color, size = label_size,
                         position = ggplot2::position_nudge(x = label_nudge_x) )
  }

  # add the hole to make the donut chart
  if(plot_donut){
    p <- p + ggplot2::xlim(0.5, 2.5)
  }

  p
}

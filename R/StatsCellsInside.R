
#' Summarize cell counts and proportions inside spatial regions
#'
#' Computes the number and proportion of cells from each cluster inside boundaries or ring regions.
#' This function is useful for downstream visualizations such as bar plots or pie charts
#' showing the spatial composition of cell types per region.
#'
#' @inheritParams PlotCellsInside
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return A data frame with one row per cluster per region, containing the following columns:
#' \itemize{
#'   \item \code{region_id}: Identifier for each spatial region.
#'   \item \code{cluster}: Cluster label of the cells.
#'   \item \code{count}: Number of cells from the given cluster in the region.
#'   \item \code{proportion}: Proportion of cells from the given cluster relative to the total number of cells in the region.
#' }
#'
#' @export
#' @examples
#' # Load coordinates data
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'                               package = "SpNeigh"))
#'
#' # Get boundary and cells inside
#' boundary_points <- GetBoundary(data = coords, one_cluster = 2,
#'                                eps = 120, minPts = 10)
#' cells_inside <- GetCellsInside(data = coords, boundary = boundary_points)
#'
#' # Summarize cluster statistics per region
#' stats_cells <- StatsCellsInside(cells_inside)
#' head(stats_cells)
#'
StatsCellsInside <- function(cells_inside = NULL) {
  # Count and proportions of cells in different clusters for each region
  cell_stats <- cells_inside %>%
    sf::st_drop_geometry() %>%
    dplyr::group_by(.data$region_id, .data$cluster) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(.data$region_id) %>%
    dplyr::mutate(proportion = .data$count / sum(.data$count)) %>%
    dplyr::ungroup()

  return(cell_stats)
}


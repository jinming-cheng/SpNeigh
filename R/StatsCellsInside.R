
#' Statistics of cells inside boundaries or rings
#'
#' @inheritParams PlotCellsInside
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
StatsCellsInside <- function(cells_inside = NULL){
  # Count and proportions of cells in different clusters for each region
  cell_stats <- cells_inside %>%
    sf::st_drop_geometry() %>%
    dplyr::group_by(.data$region_id, .data$cluster) %>%
    dplyr::summarise(count = dplyr::n(), .groups = 'drop') %>%
    dplyr::group_by(.data$region_id) %>%
    dplyr::mutate(proportion = .data$count / sum(.data$count))  %>%
    dplyr::ungroup()

  cell_stats
}


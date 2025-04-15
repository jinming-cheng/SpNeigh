
#' Statistics of cells inside boundaries or rings
#'
#' @inheritParams PlotCellsInside
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#' @examples
#' # Load coordinates
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'                               package = "SpNeigh"))
#' head(coords)
#'
#' # Get cells inside boundaries
#' boundary_points <- GetBoundary(data = coords, one_cluster = 2,
#'                                eps = 120, minPts = 10)
#' cells_inside <- GetCellsInside(data = coords, boundary =  boundary_points)
#'
#' # Statistics of cells inside boundaries
#' stats_cells = StatsCellsInside(cells_inside)
#' head(stats_cells)
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


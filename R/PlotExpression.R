
#' Plot gene expression across spatial coordinates
#'
#' Visualizes expression of selected genes on a spatial plot using cell coordinates and expression matrix.
#' Can display expression for all cells or a subset (e.g., by cluster or manually specified cells).
#' Also supports splitting the plot using metadata columns (e.g., `cluster`) and returning either
#' a combined plot or a list of individual ggplot objects.
#'
#' @inheritParams ExtractCoords
#' @inheritParams PlotBoundary
#' @importFrom methods is
#' @importFrom rlang .data
#' @param exp_mat A gene expression matrix (genes x cells) of class `matrix` or `dgCMatrix`.
#'        If `data` is a Seurat object, the matrix will be extracted automatically using `Seurat::GetAssayData()`.
#' @param genes Character vector specifying gene names to be plotted. Must match row names in `exp_mat`.
#' @param sub_plot Logical. If `TRUE`, only a subset of cells is plotted (based on `one_cluster` or `sub_cells`).
#'        Default is `FALSE`.
#' @param one_cluster Optional. Cluster ID to subset cells if `sub_plot = TRUE`.
#' @param sub_cells Optional. Vector of cell IDs to include in the plot if `sub_plot = TRUE`.
#'        If both `one_cluster` and `sub_cells` are provided, the intersection is used.
#' @param split_by Optional. Column name in the metadata (from `data`) to facet the plots (e.g., by `cluster`).
#' @param return_list Logical. If `TRUE`, returns a named list of individual ggplot objects per gene.
#'        If `FALSE` (default), plots are wrapped into a single patchwork layout.
#' @param point_size Numeric. Size of plotted points. Default is `0.2`.
#' @param shuffle Logical. If `TRUE`, shuffles cell order before plotting. Otherwise, cells with higher expression are plotted on top. Default is `FALSE`.
#'
#' @return A `patchwork` object or a list of ggplot objects if `return_list = TRUE`.
#' @export
#' @examples
#' df <- data.frame(x = c(rnorm(100, 1), rnorm(100, 5)),
#'                  y = c(rnorm(100, 1), rnorm(100, 5)),
#'                  cell = 1:200,
#'                  cluster = rep(1:2, each = 100))
#' exp_mat <- data.frame(
#'   gene1 = c(runif(100, 4, 6), runif(100, 0, 2)),
#'   gene2 = c(runif(100, 0, 1), runif(100, 5, 8))
#'   )
#'
#' exp_mat <- t(exp_mat)
#' colnames(exp_mat) <- df$cell
#'
#' set.seed(123)
#' PlotExpression(data = df, exp_mat = exp_mat,
#'                genes = c("gene1", "gene2"), point_size = 2)
#'
#' PlotExpression(data = df, exp_mat = exp_mat,
#'                genes = "gene1", sub_plot = TRUE,
#'                one_cluster = 1, point_size = 2)
#'
PlotExpression <- function(data = NULL,
                           exp_mat = NULL,
                           genes = NULL,
                           sub_plot = FALSE,
                           one_cluster = NULL,
                           sub_cells = NULL,
                           split_by = NULL,
                           ncol = NULL,
                           return_list = FALSE,
                           point_size = 0.2,
                           angle_x_label = 0,
                           shuffle = FALSE,
                           theme_ggplot = my_theme_ggplot()) {

  # Automatically extract expression matrix from Seurat object
  if (is(data, "Seurat")) {
    exp_mat <- Seurat::GetAssayData(data)
  } else {
    if (!inherits(exp_mat, c("dgCMatrix", "matrix"))) {
      stop("`exp_mat` must be a numeric matrix or dgCMatrix (genes x cells).")
    }
  }

  # Validate genes
  if (!all(genes %in% rownames(exp_mat))) {
    stop("One or more specified genes are not found in `exp_mat`.")
  }

  # Extract coordinates and metadata
  sp_coords <- ExtractCoords(data = data)

  # Check that cell names match and are in the same order
  if (!identical(as.character(sp_coords$cell), colnames(exp_mat))) {
    stop("Cell names in `data` and `exp_mat` must match and be in the same order.")
  }

  # Merge expression into coordinates
  sp_coords <- cbind(sp_coords, t(as.matrix(exp_mat[genes, , drop = FALSE])))

  # Subset if needed
  data_for_plot <- sp_coords
  if (sub_plot) {
    if (!is.null(one_cluster)) {
      data_for_plot <- dplyr::filter(data_for_plot, .data$cluster == one_cluster)
    }
    if (!is.null(sub_cells)) {
      data_for_plot <- dplyr::filter(data_for_plot, .data$cell %in% sub_cells)
    }
  }

  # Plot loop over genes
  p <- list()
  for (i in genes) {
    plot_data <- if (shuffle) {
      data_for_plot[sample(nrow(data_for_plot)), ]
    } else {
      dplyr::arrange(data_for_plot, .data[[i]])
    }

    p[[i]] <- ggplot2::ggplot(plot_data,
                              ggplot2::aes(x = .data$x, y = .data$y,
                                           color = !!rlang::sym(i))) +
      ggplot2::geom_point(size = point_size) +
      theme_ggplot +
      ggplot2::scale_color_viridis_c(option = "plasma") +
      ggplot2::theme(axis.text.x =
                       ggplot2::element_text(angle = angle_x_label, hjust = 1))

    if (!is.null(split_by)) {
      if (!(split_by %in% colnames(plot_data))) {
        stop("`split_by` must be a column in the input data.")
      }
      p[[i]] <- p[[i]] + ggplot2::facet_wrap(stats::as.formula(paste("~", split_by)), ncol = ncol)
    }
  }

  if (return_list) {
    return(p)
  } else {
    return(patchwork::wrap_plots(p))
  }
}


#' Plot average gene expression along spatial distance
#'
#' Visualizes how the average expression of specified genes varies along a spatial distance gradient.
#' The spatial distance is binned, and the average expression within each bin is plotted as a heatmap.
#'
#' @inheritParams RunSpatialDE
#' @inheritParams PlotBoundary
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @param spatial_distance A named numeric vector containing the spatial distance (or weights) for each cell.
#' @param genes Character vector specifying gene names to be plotted. Must match row names in `exp_mat`.
#' @param n_bins Integer. Number of bins to divide the spatial distance into. Default is `50`.
#' @param scale_method A string indicating how to scale the average expression values across bins for each gene.
#' Options are:
#' \describe{
#'   \item{\code{"none"}}{No scaling (default). The average expression is plotted as-is.}
#'   \item{\code{"zscore"}}{Standardize expression (mean 0, SD 1) per gene using \code{scale()}.}
#'   \item{\code{"minmax"}}{Normalize expression to \[0, 1\] range per gene using \code{scales::rescale()}.}
#' }
#' @param n_labels Integer. Number of axis labels to show along the distance axis. Default is `6`.
#' @param row_gap Numeric between 0 (inclusive) and 1 (exclusive). Gap between rows (genes) in the plot. Default is `0.1`.
#' @param column_gap Numeric between 0 (inclusive) and 1 (exclusive). Gap between columns (distance bins) in the plot. Default is `0`.
#' @param label_x Character. Label for the x-axis. Default is "Spatial distance".
#' @param label_y Character. Label for the y-axis. Default is "Gene".
#'
#' @return A ggplot2 object displaying a heatmap of binned average gene expression across spatial distances.
#'
#' @export
#'
#' @examples
#' # Example spatial expression heatmap
#' set.seed(1)
#' exp_mat <- matrix(runif(1000), nrow = 10)
#' rownames(exp_mat) <- paste0("Gene", 1:10)
#' spatial_distance <- runif(100)
#' PlotSpatialExpression(exp_mat = exp_mat,
#'                       spatial_distance = spatial_distance,
#'                       genes = rownames(exp_mat)[1:5])
#'

PlotSpatialExpression <- function(exp_mat = NULL,
                                  spatial_distance = NULL,
                                  genes = NULL,
                                  n_bins = 50,
                                  scale_method = c("none", "zscore", "minmax"),
                                  n_labels = 6,
                                  row_gap = 0.1,
                                  column_gap = 0,
                                  label_x = "Spatial distance",
                                  label_y = "Gene",
                                  theme_ggplot = my_theme_ggplot()) {

  # --- Checks ---
  if (is.null(exp_mat) || !inherits(exp_mat, c("matrix", "dgCMatrix"))) {
    stop("`exp_mat` must be a numeric matrix or dgCMatrix (genes x cells).")
  }
  if (!(row_gap >= 0 & row_gap < 1)) {
    stop("`row_gap` must be a numeric value in the interval [0,1).")
  }
  if (!(column_gap >= 0 & column_gap < 1)) {
    stop("`column_gap` must be a numeric value in the interval [0,1).")
  }
  if (ncol(exp_mat) != length(spatial_distance)) {
    stop("Number of columns in `exp_mat` must match length of `spatial_distance`.")
  }
  if (!all(genes %in% rownames(exp_mat))) {
    stop("One or more specified genes are not found in `exp_mat`.")
  }

  scale_method <- match.arg(scale_method)

  # --- Prepare Data ---
  df <- data.frame(t(as.matrix(exp_mat[genes, , drop = FALSE])),
                   distance = spatial_distance)

  df_long <- df %>%
    tidyr::pivot_longer(cols = -.data$distance,
                        names_to = "gene", values_to = "expression") %>%
    dplyr::mutate(bin = cut(.data$distance,
                            breaks = n_bins, include.lowest = TRUE)) %>%
    dplyr::group_by(.data$gene, .data$bin) %>%
    dplyr::summarise(avg_expression = mean(.data$expression),
                     avg_distance = mean(.data$distance), .groups = "drop")
  df_long$gene <- factor(df_long$gene, levels = genes)

  # Optionally scale average expression
  if (scale_method != "none") {
    df_long <- df_long %>% dplyr::group_by(.data$gene) %>%
      dplyr::mutate(avg_expression = dplyr::case_when(
        scale_method == "zscore" ~ as.numeric(scale(.data$avg_expression)),
        scale_method == "minmax" ~ scales::rescale(.data$avg_expression, to = c(0, 1)),
        TRUE ~ .data$avg_expression
      )) %>% dplyr::ungroup()
  }

  # --- Label setup for x-axis ---
  df_labels <- df_long %>% dplyr::select(.data$bin, .data$avg_distance) %>%
    dplyr::distinct() %>% dplyr::arrange(.data$avg_distance)

  bin_indices <- round(seq(1, nrow(df_labels), length.out = n_labels))
  df_labels$label <- ""
  df_labels$label[bin_indices] <- round(df_labels$avg_distance[bin_indices], 2)
  df_labels <- dplyr::filter(df_labels, .data$label != "")

  # --- Plot ---
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = .data$bin, y = .data$gene,
                                             fill = .data$avg_expression)) +
    ggplot2::geom_tile(width = 1 - column_gap, height = 1 - row_gap) +
    ggplot2::scale_x_discrete(breaks = df_labels$bin, labels = df_labels$label) +
    ggplot2::scale_fill_viridis_c(option = "plasma") +
    ggplot2::labs(x = label_x, y = label_y) + theme_ggplot

  # Set fill label based on scale method
  fill_label <- switch(scale_method,
                       "none" = "AvgExp",
                       "zscore" = "Z-scored",
                       "minmax" = "Scaled [0,1]")

  p <- p + ggplot2::labs(fill = fill_label)

  return(p)
}




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
#' @param seed Optional seed for reproducibility when `shuffle = TRUE`. Default is `123`.
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
                           return_list = FALSE,
                           point_size = 0.2,
                           angle_x_label = 0,
                           shuffle = FALSE,
                           seed = 123,
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
      set.seed(seed)
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
      p[[i]] <- p[[i]] + ggplot2::facet_wrap(stats::as.formula(paste("~", split_by)))
    }
  }

  if (return_list) {
    return(p)
  } else {
    return(patchwork::wrap_plots(p))
  }
}

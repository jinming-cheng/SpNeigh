
#' Plot expression of genes in cells on the spatial plot
#'
#' @inheritParams ExtractCoords
#' @importFrom rlang .data
#' @param exp_mat A gene expression matrix or dgCMatrix. If `data` is a Seurat object, the `exp_mat` will be automatically obtained.
#' @param genes Specify genes for plotting the expression.
#' @param sub_plot A logical value indicating whether to only plot a subset of cells. Default is FALSE.
#' @param one_cluster Specify a cluster to obtain a subset of cells for `sub_plot`.
#' @param sub_cells Specify a subset of cells for `sub_plot`. If both `one_cluster` and `sub_cells` are provided, the common cells in both the specified cluster and the specified cell subset will be used.
#' @param split_by A column name in the coordinate data.frame (from `data` parameter) used to split the plot.
#' @param return_list Logical. Whether to return a list of ggplot objects. Default is FALSE.
#' @param point_size Point size of cells. Default is 0.2.
#' @param angle_x_label Angle to rotate the x-axis labels.
#' @param theme_ggplot Theme for ggplot.
#' @param shuffle Logical. Whether to shuffle the order of cells before plotting. If `FALSE` (default), cells with larger values are plotted last, so they appear on top.
#' @export
#' @examples
#' # Sample data
#' df <- data.frame(x = c( rnorm(100, mean=1),rnorm(100, mean=5) ),
#'                  y = c( rnorm(100, mean=1),rnorm(100, mean=5) ),
#'                  cell = 1:200,
#'                  cluster = rep(rep(1:2, each = 100)) )
#' exp_mat <- data.frame(gene1 = c(runif(100, min = 4, max = 6),
#'                                 runif(100, min = 0, max = 2) ),
#'                       gene2 = c(runif(100, min = 0, max = 1),
#'                                 runif(100, min = 5, max = 8) ) )
#' exp_mat <- t(exp_mat)
#' colnames(exp_mat) <- df$cell
#'
#' # Plot expression of two genes
#' p <- PlotExpression(data = df, exp_mat = exp_mat,
#'                     genes = c("gene1","gene2"),point_size = 2)
#' p
#'
#' # Plot expression of gene1 in cluster 1
#' p <- PlotExpression(data = df, exp_mat = exp_mat,
#'                     sub_plot = TRUE, one_cluster = 1,
#'                     genes = c("gene1"),point_size = 2)
#' p + ggplot2::facet_wrap(~cluster) +
#'   ggplot2::scale_fill_viridis_c(option = "plasma",limits = c(0, 6)) +
#'   ggplot2::scale_color_viridis_c(option = "plasma",limits = c(0, 6))
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
                           theme_ggplot = my_theme_ggplot()){

  # Obtain the expression data automatically if data is a Seurat object
  if( inherits(data,"Seurat") ){
    exp_mat <- Seurat::GetAssayData(data)
  }else{

    # Check: exp_mat is a matrix or dgCMatrix
    if (!inherits(exp_mat, c("dgCMatrix","matrix") ) ) {
      stop("`exp_mat` must be a numeric matrix or a dgCMatrix with rows corresponding to genes and columns to cells."
      )}

    exp_mat <- exp_mat
  }

  # Check existence of genes
  if (!all(genes %in% rownames(exp_mat)) ) {
    stop("genes are not in the `exp_mat`. Please provide the correct gene names.")
  }

  # Extract coordinates from data
  sp_coords <- ExtractCoords(data = data)

  # Check that cell names in data and exp_mat are identical and in the same order
  cells_data <- as.character(sp_coords$cell)
  cells_expr <- as.character(colnames(exp_mat))
  if (!identical(cells_data, cells_expr)) {
    stop("Cell names in `data` and `exp_mat` must match exactly and be in the same order.")
  }

  # Add gene expression columns
  sp_coords <- cbind.data.frame(sp_coords,
                                t(as.matrix(exp_mat[genes,,drop=FALSE])) )

  # Whether to plot one cluster or all clusters
  if(sub_plot){
    if(!is.null(one_cluster)){
      # Obtain coordinates of the target cluster
      coords_sub <- dplyr::filter(sp_coords,.data$cluster == one_cluster)
    }else{
      coords_sub <- sp_coords
    }

    if(!is.null(sub_cells)){
      # Obtain coordinates of the cell subset
      coords_sub <- dplyr::filter(coords_sub,.data$cell %in% sub_cells)
    }

    data_for_plot <- coords_sub

  }else{
    data_for_plot <- sp_coords
  }

  # Plot expression of genes
  p <- list()
  for (i in genes) {
    if(shuffle){
      # Shuffle cell orders for plotting
      set.seed(123)
      data_for_plot <- data_for_plot[sample(nrow(data_for_plot)), ]

    }else{
      data_for_plot <- dplyr::arrange(data_for_plot, .data[[i]])
    }


    p[[i]] <- ggplot2::ggplot(data_for_plot,
                              ggplot2::aes(x = .data$x, y = .data$y,
                                           color = !!rlang::sym(i),
                                           fill = !!rlang::sym(i) )  ) +
      ggplot2::geom_point(shape = 21, size = point_size) +

      ggplot2::scale_color_viridis_c(option = "plasma") +
      ggplot2::scale_fill_viridis_c(option = "plasma") +
      theme_ggplot +
      ggplot2::theme(axis.text.x =  ggplot2::element_text(angle = angle_x_label, hjust = 1) ) +
      ggplot2::guides(color="none")

    # split plot a by column
    if(!is.null(split_by)){
      # Check existence
      if (!(split_by %in% colnames(data_for_plot)) ) {
        stop("`split_by` should be a column name in the data for plot.")
      }
      p[[i]] <- p[[i]] + ggplot2::facet_wrap(stats::as.formula(paste("~", split_by)))
    }
  }

  if(return_list){
    return(p)
  }else{
    patchwork::wrap_plots(p)
  }

}

#' Differential expression analysis between two groups of cells using limma
#'
#' Performs differential expression analysis between two groups of cells
#' using the `limma` linear modeling framework.
#' Supports optional observation-level weights (e.g., spatial weights)
#' and filters genes by minimum expression threshold across groups.
#'
#' @param exp_mat A normalized gene expression matrix (genes x cells),
#'                either a `matrix` or `dgCMatrix`.
#'                Typically log-normalized counts, e.g., from a Seurat object.
#' @param cells_reference A character vector of cell IDs to use as the
#'                        reference (baseline) group.
#' @param cells_target A character vector of cell IDs to use as the
#'                     target (comparison) group.
#' @param weights Optional numeric vector of observation-level weights.
#'                Must be named with cell IDs and
#'                match the length of `cells_reference + cells_target`.
#' @param adj_p.value Adjusted p-value threshold for reporting
#'                    differentially expressed genes. Default is `0.05`.
#' @param min.pct Minimum proportion of cells expressing the gene in
#'                either group (values between 0 and 1).
#'                Genes not meeting this threshold are excluded before testing.
#'                Default is `0`.
#'
#' @return A data frame with differentially expressed genes,
#'         sorted by absolute log fold change.
#'         Includes columns:
#' \describe{
#'   \item{logFC}{Log2 fold change of expression (target vs. reference)}
#'   \item{AveExpr}{Average expression across both groups}
#'   \item{t, P.Value, adj.P.Val, B}{Statistical results from `limma`}
#'   \item{pct.reference}{Proportion of reference cells expressing the gene}
#'   \item{pct.target}{Proportion of target cells expressing the gene}
#'   \item{gene}{Gene name (rownames from `exp_mat`)}
#' }
#'
#' @export
#'
#' @examples
#' # Load coordinates and log-normalized expression data
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'     package = "SpNeigh"
#' ))
#' logNorm_expr <- readRDS(system.file("extdata", "LogNormExpr.rds",
#'     package = "SpNeigh"
#' ))
#'
#' # Subset cells from cluster 0 and 2
#' cells_ref <- subset(coords, cluster == 0)$cell
#' cells_tar <- subset(coords, cluster == 2)$cell
#'
#' # Run differential expression with minimum expression threshold
#' tab <- RunLimmaDE(
#'     exp_mat = logNorm_expr,
#'     cells_reference = cells_ref,
#'     cells_target = cells_tar,
#'     min.pct = 0.25
#' )
#'
#' head(tab[, c("gene", "logFC", "adj.P.Val", "pct.reference", "pct.target")])
#'
RunLimmaDE <- function(
    exp_mat = NULL,
    cells_reference = NULL,
    cells_target = NULL,
    weights = NULL,
    adj_p.value = 0.05,
    min.pct = 0) {
    # --- Input checks ---
    if (is.null(exp_mat) || !inherits(exp_mat, c("matrix", "dgCMatrix"))) {
        stop(
            "'exp_mat' must be a numeric",
            " matrix or dgCMatrix (genes x cells)."
        )
    }

    if (is.null(cells_reference) || is.null(cells_target)) {
        stop("Both 'cells_reference' and 'cells_target' must be provided.")
    }

    all_cells <- c(cells_reference, cells_target)
    if (!all(all_cells %in% colnames(exp_mat))) {
        stop("Some cell IDs not found in column names of 'exp_mat'.")
    }

    if (!is.null(weights)) {
        if (is.null(names(weights)) || length(weights) != length(all_cells)) {
            stop(
                "'weights' must be a named numeric vector",
                " with the same length as all selected cells."
            )
        }
        weights <- weights[all_cells] # ensure correct order
    }

    # --- Convert to sparse matrix if needed ---
    if (inherits(exp_mat, "matrix")) {
        exp_mat <- Matrix::Matrix(exp_mat, sparse = TRUE)
    }

    # --- Subset expression matrix ---
    expr <- exp_mat[, all_cells, drop = FALSE]

    # --- Compute percentage of expressing cells ---
    df_pct <- data.frame(
        pct.reference = Matrix::rowMeans(expr[, cells_reference] > 0),
        pct.target = Matrix::rowMeans(expr[, cells_target] > 0)
    )
    rownames(df_pct) <- rownames(expr)

    # --- Filter genes ---
    keep_genes <- rowSums(df_pct >= min.pct) >= 1
    expr <- expr[keep_genes, , drop = FALSE]
    df_pct <- df_pct[keep_genes, , drop = FALSE]

    if (nrow(expr) == 0) {
        stop("No genes passed the 'min.pct' filter.")
    }

    # --- Build design matrix ---
    group <- factor(c(
        rep(1, length(cells_reference)),
        rep(2, length(cells_target))
    ))
    design <- stats::model.matrix(~group)

    # --- Fit model using limma ---
    fit <- limma::lmFit(expr, design, weights = weights)
    fit <- limma::eBayes(fit)
    tab <- limma::topTable(fit, coef = 2, number = Inf, p.value = adj_p.value)

    if (nrow(tab) == 0) {
        warning("No genes passed the adjusted p-value threshold.")
        return(tab)
    }

    # --- Add gene-level annotations ---
    tab <- cbind(tab, df_pct[rownames(tab), ])
    tab$gene <- rownames(tab)
    tab <- tab[order(-abs(tab$logFC)), ]

    return(tab)
}



#' Differential expression along spatial distance gradients using splines
#'
#' Performs spatially-aware differential expression (DE) analysis by
#' modeling gene expression as a smooth function of a continuous
#' spatial covariate (e.g., distance to a boundary or centroid).
#' Natural spline basis functions are used to capture non-linear trends
#' in expression relative to spatial distance.
#' This method is suitable for identifying genes whose expression varies
#' continuously across spatial structures.
#'
#' @inheritParams RunLimmaDE
#' @inheritParams SplineDesign
#' @param cell_ids A character vector of cell IDs (column names of `exp_mat`)
#'                 used for DE analysis.
#' @param spatial_distance A named numeric vector containing the spatial
#'                         distance (or weights) for each cell.
#'                         Must be the same length as `cell_ids`.
#'                         Scaled distances are recommended.
#'
#' @return A data frame of differentially expressed genes, including:
#' \describe{
#'   \item{AveExpr, F, P.Value, adj.P.Val}{limma differential expression
#'                                         outputs}
#'   \item{Z1, Z2, Z3}{Spline coefficients (Z1 typically corresponds
#'                     to linear trend)}
#'   \item{gene}{Gene name (from `exp_mat`)}
#'   \item{trend}{"Positive" or "Negative" trend based on the sign of Z1}
#' }
#' The first spline coefficient (`Z1`) captures the main expression trend
#' along the spatial distance.
#'
#'
#' @export
#'
#' @examples
#' # Load example data
#' coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
#'     package = "SpNeigh"
#' ))
#' logNorm_expr <- readRDS(system.file("extdata", "LogNormExpr.rds",
#'     package = "SpNeigh"
#' ))
#'
#' # Identify cluster-specific cells and compute spatial weights
#' cells_c0 <- subset(coords, cluster == 0)$cell
#' bon_c0 <- GetBoundary(data = coords, one_cluster = 0)
#' weights <- ComputeBoundaryWeights(
#'     data = coords, cell_ids = cells_c0,
#'     boundary = bon_c0
#' )
#'
#' # Run spatial DE
#' result <- RunSpatialDE(
#'     exp_mat = logNorm_expr, cell_ids = cells_c0,
#'     spatial_distance = weights
#' )
#' head(result)
#'
RunSpatialDE <- function(
    exp_mat = NULL,
    cell_ids = NULL,
    spatial_distance = NULL,
    adj_p.value = 0.05,
    df = 3) {
    # --- Input checks ---
    if (is.null(exp_mat) || !inherits(exp_mat, c("matrix", "dgCMatrix"))) {
        stop("'exp_mat' must be a non-null numeric matrix (genes x cells).")
    }

    if (is.null(cell_ids)) {
        stop("'cell_ids' must be provided.")
    }

    if (!all(cell_ids %in% colnames(exp_mat))) {
        stop("Some cell IDs not found in column names of 'exp_mat'.")
    }

    all_cells <- cell_ids


    if (is.null(spatial_distance) || is.null(names(spatial_distance)) ||
        length(spatial_distance) != length(all_cells)) {
        stop(
            "'spatial_distance' must be a named vector, and its",
            " length must match the number of selected cells."
        )
    }


    # --- Make spatial_distance numeric and in correct order ---
    t1 <- as.numeric(spatial_distance[all_cells])

    # --- Subset expression matrix ---
    expr <- exp_mat[, all_cells, drop = FALSE]

    # --- Create spline-based design matrix ---
    Z <- SplineDesign(t1, df = df)
    design <- stats::model.matrix(~Z)

    # --- Run limma ---
    fit <- limma::lmFit(expr, design)
    fit <- limma::eBayes(fit)
    tab <- limma::topTable(fit,
        coef = 2:(df + 1), number = Inf,
        p.value = adj_p.value
    )

    if (nrow(tab) == 0) {
        warning("No genes passed the adjusted p-value threshold.")
        return(tab)
    }

    # --- Annotate results ---
    tab$gene <- rownames(tab)
    tab$trend <- ifelse(tab$Z1 > 0, "Positive", "Negative")

    return(tab)
}


#' Generate an orthonormal spline-based design matrix
#'
#' Constructs an orthonormal design matrix from a numeric covariate
#' (e.g., spatial distance) using a natural cubic spline basis.
#' The output matrix can be used in linear modeling to capture smooth,
#' non-linear trends along continuous variables.
#'
#' The first column of the resulting matrix is aligned to show
#' a positive correlation with the input vector and typically captures
#' the main linear or monotonic trend.
#'
#' @param x A numeric vector representing a continuous covariate
#'          (e.g., distance or pseudotime).
#' @param df Integer. Degrees of freedom for the spline basis. Default is 3.
#'
#' @return A numeric matrix with orthonormal columns
#'         (same number of rows as `x`). The columns represent smoothed
#'         trends extracted from the spline basis. The first column is
#'         directionally aligned with the input vector (`x`).
#'
#' @importFrom splines ns
#' @importFrom stats cor
#'
#' @export
#'
#' @examples
#' x <- seq(0, 1, length.out = 100)
#' Z <- SplineDesign(x)
#' cor(Z[, 1], x) # Should be > 0
#'
#' # Use Z in modeling
#' y <- sin(2 * pi * x) + rnorm(100, sd = 0.2)
#' fit <- lm(y ~ Z)
#' summary(fit)
SplineDesign <- function(x, df = 3) {
    if (!is.numeric(x)) {
        stop("'x' must be a numeric vector.")
    }

    # Generate spline basis
    X <- splines::ns(as.numeric(x), df = df)

    # Augment matrix with intercept and original variable
    A <- cbind(1, x, X)

    # Orthonormalize via QR decomposition
    QR <- qr(A)
    r <- QR$rank
    R_rank <- QR$qr[seq(from = 1, to = r), seq(from = 1, to = r)]
    Z <- t(backsolve(R_rank, t(A), transpose = TRUE))

    # Remove intercept column
    Z <- Z[, -1]

    # Ensure first component is positively correlated with input
    if (stats::cor(Z[, 1], x) < 0) {
        Z[, 1] <- -Z[, 1]
    }

    return(Z)
}

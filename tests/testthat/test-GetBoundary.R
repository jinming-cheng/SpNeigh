test_that("Test removeOutliers", {
    set.seed(123)
    x1 <- rnorm(20, mean = 5, sd = 0.5)
    y1 <- rnorm(20, mean = 5, sd = 0.5)
    x2 <- runif(5, min = 10, max = 15)
    y2 <- runif(5, min = 10, max = 15)
    coords <- data.frame(x = c(x1, x2), y = c(y1, y2))

    new_coords <- removeOutliers(coords, k = 5, distance_cutoff = 2)

    expect_true(is.data.frame(new_coords))

    expect_error(removeOutliers(x1))
})


test_that("removeOutliers removes distant point", {
    coords <- data.frame(
        x = c(0, 0, 1, 1, 100),
        y = c(0, 1, 0, 1, 100)
    )

    res <- removeOutliers(coords, k = 2, distance_cutoff = 10)

    expect_equal(nrow(res), 4)
    expect_false(any(res$x == 100))
})


test_that("Test extractCoords", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    coords <- extractCoords(coords)

    expect_true(is.data.frame(coords))

    expect_error(extractCoords(coords[, c(1, 2)]))

    expect_error(extractCoords(111))

    coords_sub <- subset(coords, cluster %in% c("0", "2"))
    coords_sub <- as.matrix(coords_sub[, c("x", "y")])
    metadata_sub <- subset(
        coords[, c("cell", "cluster")],
        cluster %in% c("0", "2")
    )
    logNorm_expr <- readRDS(system.file("extdata", "LogNormExpr.rds",
        package = "SpNeigh"
    ))
    spe <- SpatialExperiment::SpatialExperiment(
        assay = list("logcounts" = logNorm_expr),
        colData = metadata_sub,
        spatialCoords = coords_sub
    )

    coords <- extractCoords(spe)

    expect_true(is.data.frame(coords))

    expect_error(extractCoords(spe, cluster_col = "unknown"))

    seu_sp <- Seurat::CreateSeuratObject(
        assay = "Spatial",
        counts = logNorm_expr,
        meta.data = metadata_sub
    )
    SeuratObject::LayerData(seu_sp,
        assay = "Spatial",
        layer = "data"
    ) <- logNorm_expr
    cents <- SeuratObject::CreateCentroids(coords_sub[, c("x", "y")])
    fov <- SeuratObject::CreateFOV(
        coords = list("centroids" = cents),
        type = c("centroids"),
        assay = "Spatial"
    )
    seu_sp[["fov"]] <- fov
    seu_sp@meta.data$seurat_clusters <- seu_sp$cluster

    coords <- extractCoords(seu_sp)

    expect_error(extractCoords(seu_sp, cluster_col = "unknown"))

    expect_true(is.data.frame(coords))
})


test_that("extractCoords works for data.frame input", {
    df <- data.frame(
        cell = c("c1", "c2"),
        x = c(1, 2),
        y = c(3, 4),
        cluster = c("A", "B")
    )

    res <- extractCoords(df)

    expect_true(is.data.frame(res))
    expect_equal(nrow(res), 2)
    expect_equal(colnames(res), c("cell", "x", "y", "cluster"))
    expect_equal(res$cell, df$cell)
    expect_equal(res$cluster, df$cluster)
})


test_that("extractCoords errors when required columns missing", {
    df <- data.frame(
        cell = c("c1", "c2"),
        x = c(1, 2),
        y = c(3, 4)
    )

    expect_error(
        extractCoords(df, extract_cluster = TRUE),
        "missing required columns"
    )
})


test_that("extractCoords respects extract_cluster = FALSE", {
    df <- data.frame(
        cell = c("c1", "c2"),
        x = c(1, 2),
        y = c(3, 4),
        cluster = c("A", "B")
    )

    res <- extractCoords(df, extract_cluster = FALSE)

    expect_equal(colnames(res), c("cell", "x", "y"))
})

test_that("extractCoords errors on unsupported class", {
    expect_error(extractCoords(1:10))
})


test_that("extractCoords works for SpatialExperiment", {
    skip_if_not_installed("SpatialExperiment")

    mat <- matrix(rpois(20, 5), nrow = 5)
    rownames(mat) <- paste0("gene", 1:5)
    colnames(mat) <- paste0("cell", 1:4)

    spe <- SpatialExperiment::SpatialExperiment(
        assays = list(counts = mat),
        spatialCoords = cbind(x = 1:4, y = 5:8)
    )

    SummarizedExperiment::colData(spe)$cluster <- c("A", "B", "A", "B")

    res <- extractCoords(spe)

    expect_equal(nrow(res), 4)
    expect_true("cluster" %in% colnames(res))
})


test_that("extractCoords works for Seurat objects", {
    skip_if_not_installed("Seurat")

    # minimal Seurat object
    mat <- matrix(rpois(20, 5), nrow = 5)
    rownames(mat) <- paste0("gene", 1:5)
    colnames(mat) <- paste0("cell", 1:4)

    mat <- Matrix::Matrix(mat, sparse = TRUE)
    seu <- Seurat::CreateSeuratObject(counts = mat, assay = "Spatial")
    seu$seurat_clusters <- c("A", "B", "A", "B")

    # fake spatial coordinates
    coords <- data.frame(
        x = 1:4,
        y = 5:8,
        row.names = colnames(seu)
    )

    cents <- SeuratObject::CreateCentroids(coords[, c("x", "y")])

    fov <- SeuratObject::CreateFOV(
        coords = list("centroids" = cents),
        type = c("centroids"),
        assay = "Spatial"
    )

    seu[["fov"]] <- fov

    res <- extractCoords(seu)

    expect_true(all(c("cell", "x", "y", "cluster") %in% colnames(res)))
    expect_equal(nrow(res), 4)
})


test_that("Test getBoundary related functions", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    boundary_points <- getBoundary(
        data = coords,
        one_cluster = 2,
        subregion_method = "kmeans"
    )

    expect_true(is.data.frame(boundary_points))

    boundary_points <- getBoundary(
        data = coords,
        one_cluster = 2,
        multi_region = FALSE
    )

    expect_true(is.data.frame(boundary_points))

    expect_error(getBoundary(data = coords, one_cluster = NULL))

    expect_error(getBoundary(data = coords, one_cluster = "abc"))

    expect_error(getBoundary(
        data = coords,
        one_cluster = 2, distance_cutoff = 0
    ))

    expect_error(getBoundary(
        data = coords,
        one_cluster = 2, distance_cutoff = 8
    ))

    boundary_points <- getBoundary(data = coords, one_cluster = 2)

    expect_true(is.data.frame(boundary_points))

    boundary_polys <- buildBoundaryPoly(boundary_points)

    expect_error(buildBoundaryPoly(boundary = NULL))

    expect_s3_class(boundary_polys, "sf")
    # 3 for tutorial data
    expect_equal(nrow(boundary_polys), 3)
    expect_true(all(sf::st_geometry_type(boundary_polys) == "POLYGON"))

    expect_warning(boundaryPolyToPoints(boundary_polys))

    expect_error(boundaryPolyToPoints(boundary_poly = boundary_points))

    outer <- getOuterBoundary(boundary_points)
    # Area must increase
    expect_true(all(sf::st_area(outer) > sf::st_area(boundary_polys)))

    expect_s3_class(outer, "sf")

    outer <- getOuterBoundary(boundary_polys)

    expect_s3_class(outer, "sf")

    expect_warning(getInnerBoundary(boundary_polys))

    expect_error(getInnerBoundary(boundary_polys, dist = -10))

    inner <- suppressWarnings(getInnerBoundary(boundary_polys))

    expect_error(boundaryPolyToPoints(inner))

    rings <- getRingRegion(boundary = boundary_points, outer_boundary = outer)

    expect_s3_class(rings, "sf")
    expect_equal(nrow(rings), nrow(boundary_polys))
    expect_true(all(sf::st_area(rings) > 0))

    expect_equal(
        as.numeric(sf::st_area(outer) - sf::st_area(boundary_polys)),
        as.numeric(sf::st_area(rings)),
        tolerance = 1e-6
    )

    rings <- getRingRegion(boundary = boundary_polys, outer_boundary = NULL)

    expect_s3_class(rings, "sf")

    boundary_edges <- splitBoundaryPolyByAnchor(boundary_polys[1, ])

    expect_s3_class(boundary_edges, "sf")
    expect_equal(nrow(boundary_edges), 2)
    expect_true(all(sf::st_geometry_type(boundary_edges) %in%
        c("LINESTRING", "POLYGON")))

    expect_error(splitBoundaryPolyByAnchor(boundary_points))

    expect_error(splitBoundaryPolyByAnchor(boundary_polys[1, ],
        pt1 = "a", pt2 = "b"
    ))
})

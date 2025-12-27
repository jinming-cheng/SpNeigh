test_that("Test RemoveOutliers", {
    set.seed(123)
    x1 <- rnorm(20, mean = 5, sd = 0.5)
    y1 <- rnorm(20, mean = 5, sd = 0.5)
    x2 <- runif(5, min = 10, max = 15)
    y2 <- runif(5, min = 10, max = 15)
    coords <- data.frame(x = c(x1, x2), y = c(y1, y2))

    new_coords <- RemoveOutliers(coords, k = 5, distance_cutoff = 2)

    expect_true(is.data.frame(new_coords))

    expect_error(RemoveOutliers(x1))
})


test_that("Test ExtractCoords", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    coords <- ExtractCoords(coords)

    expect_true(is.data.frame(coords))

    expect_error(ExtractCoords(coords[, c(1, 2)]))

    expect_error(ExtractCoords(111))

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

    coords <- ExtractCoords(spe)

    expect_true(is.data.frame(coords))

    expect_error(ExtractCoords(spe, cluster_col = "unknown"))

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

    coords <- ExtractCoords(seu_sp)

    expect_error(ExtractCoords(seu_sp, cluster_col = "unknown"))

    expect_true(is.data.frame(coords))
})


test_that("Test GetBoundary related functions", {
    coords <- readRDS(system.file("extdata", "MouseBrainCoords.rds",
        package = "SpNeigh"
    ))

    boundary_points <- GetBoundary(
        data = coords,
        one_cluster = 2,
        subregion_method = "kmeans"
    )

    expect_true(is.data.frame(boundary_points))

    boundary_points <- GetBoundary(
        data = coords,
        one_cluster = 2,
        multi_region = FALSE
    )

    expect_true(is.data.frame(boundary_points))

    expect_error(GetBoundary(data = coords, one_cluster = NULL))

    expect_error(GetBoundary(data = coords, one_cluster = "abc"))

    expect_error(GetBoundary(
        data = coords,
        one_cluster = 2, distance_cutoff = 0
    ))

    expect_error(GetBoundary(
        data = coords,
        one_cluster = 2, distance_cutoff = 8
    ))

    boundary_points <- GetBoundary(data = coords, one_cluster = 2)

    expect_true(is.data.frame(boundary_points))

    boundary_polys <- BuildBoundaryPoly(boundary_points)

    expect_error(BuildBoundaryPoly(boundary = NULL))

    expect_true(inherits(boundary_polys, "sf"))

    expect_warning(BoundaryPolyToPoints(boundary_polys))

    expect_error(BoundaryPolyToPoints(boundary_poly = boundary_points))

    outer <- GetOuterBoundary(boundary_points)

    expect_true(inherits(outer, "sf"))

    outer <- GetOuterBoundary(boundary_polys)

    expect_true(inherits(outer, "sf"))

    expect_warning(GetInnerBoundary(boundary_polys))

    expect_error(GetInnerBoundary(boundary_polys, dist = -10))

    inner <- suppressWarnings(GetInnerBoundary(boundary_polys))

    expect_error(BoundaryPolyToPoints(inner))

    rings <- GetRingRegion(boundary = boundary_points, outer_boundary = outer)

    expect_true(inherits(rings, "sf"))

    rings <- GetRingRegion(boundary = boundary_polys, outer_boundary = NULL)

    expect_true(inherits(rings, "sf"))

    boundary_edges <- SplitBoundaryPolyByAnchor(boundary_polys[1, ])

    expect_true(inherits(boundary_edges, "sf"))

    expect_error(SplitBoundaryPolyByAnchor(boundary_points))

    expect_error(SplitBoundaryPolyByAnchor(boundary_polys[1, ],
        pt1 = "a", pt2 = "b"
    ))
})

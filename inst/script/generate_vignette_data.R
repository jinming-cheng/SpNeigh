
# 2025.12.27
# R codes for generating vignette data

#######################
#                     #
# Seurat analysis of the mouse brain Xenium data

# download data
# The (Fresh Frozen) Mouse Brain Tiny Xenium data from 10x Genomics website
# was used in this analysis:
#   https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-for-xenium-explorer-demo-1-standard.

# read in data
path <- paste0("Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs")
xenium.obj <- LoadXenium(path)

# remove low quality cells
xenium.obj <- subset(xenium.obj,
                     subset = (nFeature_Xenium >= 10) &
                       (nCount_Xenium >= 15) )

# run Seurat pipeline
xenium.obj <- NormalizeData(xenium.obj, assay = "Xenium")
xenium.obj <- ScaleData(xenium.obj, assay = "Xenium")
xenium.obj <- RunPCA(xenium.obj, npcs = 30, features = rownames(xenium.obj))
xenium.obj <- RunUMAP(xenium.obj, dims = 1:30)
xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:30)
xenium.obj <- FindClusters(xenium.obj, resolution = 0.1)
#                     #
#######################

#######################
#                     #
# Subset Xenium data for vignette

# Extract and save coordinates
coords <- SpNeigh::ExtractCoords(xenium.obj)
saveRDS(coords, file = "MouseBrainCoords.rds", compress = "xz")

# Extract and save log-normalized expression data of clusters 0 and 2
logNormExpr <- Seurat::GetAssayData(xenium.obj)
cells_keep <- colnames(xenium.obj)[(xenium.obj$seurat_clusters %in% c(0,2) )]
logNormExpr <- logNormExpr[,cells_keep]
saveRDS(logNormExpr, file = "LogNormExpr.rds", compress = "xz")

#                     #
#######################

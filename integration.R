library(Seurat)
library(harmony)
library(slingshot)
library(SingleCellExperiment)

#load sox9 seurat objects of the subclusters of interest
#verify if the subclusters are okay
sox95 <- readRDS("~/BINF/yushi scrnaseq/New folder/rds objects/subclusters/sox95.rds")
sox105 <- readRDS("~/BINF/yushi scrnaseq/New folder/rds objects/subclusters/sox105.rds")
sox115 <- readRDS("~/BINF/yushi scrnaseq/New folder/rds objects/subclusters/sox115.rds")

sox95$eday = "E9.5"
sox105$eday = "E10.5"
sox115$eday = "E11.5"

soxCombined = merge(sox95, y = list(sox105, sox115))
soxCombined = NormalizeData(soxCombined)
soxCombined = FindVariableFeatures(soxCombined)
soxCombined = ScaleData(soxCombined)
soxCombined = RunPCA(soxCombined)

soxCombined = RunHarmony(soxCombined, group.by.vars = "eday")
soxCombined = RunUMAP(soxCombined, reduction = "harmony", dims = 1:30)
soxCombined <- FindNeighbors(soxCombined, reduction = "harmony", dims = 1:30)
soxCombined = FindClusters(soxCombined, resolution = 0.5) #change later

DimPlot(soxCombined, reduction = "umap", group.by = "eday")
DimPlot(soxCombined, reduction = "umap", label = TRUE)

#slingshot
#joinLayers coz GetAssayDAta doesn't work for multiple layers of seurat v5
sce = soxCombined %<>% JoinLayers() %>% as.SingleCellExperiment()
sce <- slingshot(sce, clusterLabels = 'seurat_clusters', reducedDim = 'UMAP')

plot(reducedDims(sce)$UMAP, col = as.factor(sce$seurat_clusters), pch = 16, asp = 1,
     main = "Trajectory Analysis with Slingshot")
lines(SlingshotDataSet(sce), lwd = 2, col = 'black')



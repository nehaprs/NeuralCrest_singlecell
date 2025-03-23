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
#sox.markers <- FindAllMarkers(soxCombined, only.pos = TRUE)

soxCombined = JoinLayers(soxCombined)
resolution.range <- seq(from = 0, to = 1, by = 0.1)

# Loop over each resolution
for (res in resolution.range) {
  # Perform clustering with the current resolution
  soxCombined<- FindClusters(soxCombined, resolution = res)
  
  # Find all markers for the clusters at this resolution
  #sox.markers <- FindAllMarkers(soxCombined, only.pos = TRUE)
  
  # Define the file name for saving the markers
  #file_name <- paste0("markers_resolution_", res, ".xlsx")
  
  # Save the markers as an Excel file
  #write_xlsx(sox.markers, file_name)
  
  # Print a message to confirm completion for each resolution
  #print(paste("Markers for resolution", res, "saved to", file_name))
}

saveRDS(soxCombined, "soxCombined.rds")
clustree(soxCombined) #choose res 0.5

DimPlot(soxCombined, reduction = "umap", group.by = "eday")
DimPlot(soxCombined, reduction = "umap", label = TRUE)

rm(sox105, sox115, sox95)


setwd("~/BINF/yushi scrnaseq/time series/harmony_slingshot")
saveRDS(soxCombined, "soxCombined.rds")

#annotate the clusters using combined reference

#combine the references
tomeE9.5 <- readRDS("~/BINF/yushi scrnaseq/E9.5/tomeE9.5.rds")
tomeE10.5 <- readRDS("~/BINF/yushi scrnaseq/E10.5/tome_E10.5.rds")
tomeE11.5 <- readRDS("~/BINF/yushi scrnaseq/E11.5/tomeRef_E11.5.rds") 










#slingshot
#joinLayers coz GetAssayDAta doesn't work for multiple layers of seurat v5
sce = soxCombined %<>% JoinLayers() %>% as.SingleCellExperiment()
sce <- slingshot(sce, clusterLabels = 'seurat_clusters', reducedDim = 'UMAP')

plot(reducedDims(sce)$UMAP, col = as.factor(sce$seurat_clusters), pch = 16, asp = 1,
     main = "Trajectory Analysis with Slingshot")
lines(SlingshotDataSet(sce), lwd = 2, col = 'black')

#annotate the clusters using combined reference

#combine the references
tomeE9.5 <- readRDS("~/BINF/yushi scrnaseq/E9.5/tomeE9.5.rds")
tomeE10.5 <- readRDS("~/BINF/yushi scrnaseq/E10.5/tome_E10.5.rds")
tomeE11.5 <- readRDS("~/BINF/yushi scrnaseq/E11.5/tomeRef_E11.5.rds")

refCombined = merge(tomeE9.5, y = list(tomeE10.5, tomeE11.5))


#save refCombined to run rest of the pipeline on biomix
setwd("~/BINF/yushi scrnaseq/time series")
saveRDS(refCombined,"refCombined.rds")


refCombined = NormalizeData(refCombined)
refCombined = FindVariableFeatures(refCombined)
refCombined = ScaleData(refCombined)
refCombined = RunPCA(soxCombined)

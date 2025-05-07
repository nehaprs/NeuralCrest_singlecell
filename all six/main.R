#integrate three pairs of data first
#find double positives
#then do monocle pseudotime
library(dplyr)
library(Seurat)
library(patchwork)
library(writexl)
library(readxl)
library(clustree)
setwd("~/BINF/yushi scrnaseq/all six")
source("~/GitHub/NeuralCrest_singlecell/all six/filter_for_dp.R")

#load sox9
sproc <- readRDS("~/BINF/yushi scrnaseq/New folder/rds objects/all clusters/sox105.rds")

#load rawdata
sraw = ReadMtx(
mtx = "~/BINF/yushi scrnaseq/E10.5/sox9/1-sox9-10-5_clean_RSEC_MolsPerCell_MEX/matrix.mtx.gz",
cells = "~/BINF/yushi scrnaseq/E10.5/sox9/1-sox9-10-5_clean_RSEC_MolsPerCell_MEX/barcodes.tsv.gz",
features = "~/BINF/yushi scrnaseq/E10.5/sox9/1-sox9-10-5_clean_RSEC_MolsPerCell_MEX/features.tsv.gz")

sraw = CreateSeuratObject(sraw)


dpSox105 = filter_for_dp(sproc, sraw)

#load Pax3
pproc <- readRDS("~/BINF/yushi scrnaseq/New folder/rds objects/all clusters/pax105.rds")

#load rawdata
praw = ReadMtx(
  mtx = "~/BINF/yushi scrnaseq/E10.5/Pax3/3-pax3-10-5_clean_RSEC_MolsPerCell_MEX/matrix.mtx.gz",
  cells = "~/BINF/yushi scrnaseq/E10.5/Pax3/3-pax3-10-5_clean_RSEC_MolsPerCell_MEX/barcodes.tsv.gz",
  features = "~/BINF/yushi scrnaseq/E10.5/Pax3/3-pax3-10-5_clean_RSEC_MolsPerCell_MEX/features.tsv.gz")

praw = CreateSeuratObject(praw)


dpPax105 = filter_for_dp(pproc, praw)


######double positive objects created
##merge dp objects
s.list = list(dpSox105, dpPax105)
# normalize and identify variable features for each dataset independently
s.list = lapply(X= s.list, 
                FUN = function(x){
                  x = NormalizeData(x)
                  x = FindVariableFeatures(x, selection.method = "vst")
                })

features <- SelectIntegrationFeatures(object.list = s.list)

####
s.list <- lapply(X = s.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

s.anchors <- FindIntegrationAnchors(object.list = s.list, anchor.features = features, reduction = "rpca")
# this command creates an 'integrated' data assay
s.combined <- IntegrateData(anchorset = s.anchors)

#run a single integrated analysis
DefaultAssay(s.combined) <- "integrated"

# Run the standard workflow for visualization and clustering


s.combined <- ScaleData(s.combined, verbose = FALSE)
s.combined <- RunPCA(s.combined, npcs = 30, verbose = FALSE)
elbow = ElbowPlot(s.combined) #16
s.combined <- RunUMAP(s.combined, reduction = "pca", dims = 1:30)
s.combined <- FindNeighbors(s.combined, reduction = "pca", dims = 1:30)
s.combined <- FindClusters(s.combined, resolution = 0.5)

p1 <- DimPlot(s.combined, reduction = "umap", group.by = "orig.ident")+ggtitle("Combined Dataset at E10.5 by Origin")
p2 <- DimPlot(s.combined, reduction = "umap", group.by = "predicted.id", label = TRUE,
              repel = TRUE, pt.size = 1.5) + NoLegend() +ggtitle("Combined Dataset at E10.5")
saveRDS(s.combined,"e105Combined.rds")

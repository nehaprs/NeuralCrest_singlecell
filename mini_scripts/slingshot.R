#slingshot

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(slingshot)
library(SingleCellExperiment)
setwd("~/BINF/yushi scrnaseq/time series/harmony_slingshot/paxfull")
s.processed = readRDS("paxFullCombined_procesd.rds")


s.processed[["RNA"]] = JoinLayers(s.processed[["RNA"]])
sce = as.SingleCellExperiment(s.processed)


#set initial clusters: 

clusters = s.processed$cell_state
sce$clusters = clusters

##
eday95_clusters = unique(sce$clusters[s.processed$eday == "E9.5"])

# Run Slingshot using UMAP embeddings 
sce <- slingshot(sce, clusterLabels = 'clusters', reducedDim = 'UMAP', start.clus = eday95_clusters)

#solvingerror: nan values in umap embeddings


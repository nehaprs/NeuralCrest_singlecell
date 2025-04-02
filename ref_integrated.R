library(Seurat)
#library(writexl)
library(dplyr)
#library(ggplot2)
library(harmony)
ref95 <- readRDS("~/BINF/yushi scrnaseq/E9.5/tomeE9.5.rds")
ref105 <- readRDS("~/BINF/yushi scrnaseq/E10.5/tome_E10.5.rds")
ref115 <- readRDS("~/BINF/yushi scrnaseq/E11.5/tomeRef_E11.5.rds")

ref95$eday = "E9.5"
ref105$eday = "E10.5"
ref115$eday = "E11.5"


refCombined = merge(ref95, y = list(ref105, ref115))
rm(ref95, ref105, ref115)
refCombined = NormalizeData(refCombined)
refCombined = FindVariableFeatures(refCombined)
refCombined = ScaleData(refCombined)
refCombined = RunPCA(refCombined)

refCombined = RunHarmony(refCombined, group.by.vars = "eday")
refCombined = RunUMAP(refCombined, reduction = "harmony", dims = 1:30)
refCombined <- FindNeighbors(refCombined, reduction = "harmony", dims = 1:30)
refCombined = FindClusters(refCombined, resolution = 0.5) #change later


refCombined = JoinLayers(refCombined)
#==========================================
#analysis of all the six datasets together
#==========================================


#integrate three pairs of data first
#find double positives
#then do monocle pseudotime
library(dplyr)
library(Seurat)
library(patchwork)
library(writexl)
library(readxl)
library(clustree)
library(monocle3)
library(harmony)
library(slingshot)
library(SingleCellExperiment)
library(SeuratWrappers)
setwd("~/BINF/yushi scrnaseq/all six")


#=====================================================
#1. integrate across cell lines
#=====================================================

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



#=====================================================
#2. integrate across time points
#uses harmony
#edit: 6.16.2025:cell cycle regression
#=====================================================

#code similar to top portion of integration.R
setwd("~/BINF/yushi scrnaseq/all six/threshold0/harmony/CCregressed_monocle")

#load seurat objects
e95Combined <- readRDS("~/BINF/yushi scrnaseq/all six/threshold0/e95/e95Combined.rds")
e105Combined <- readRDS("~/BINF/yushi scrnaseq/all six/threshold0/e105/e105Combined.rds")
e115Combined <- readRDS("~/BINF/yushi scrnaseq/all six/threshold0/e115/e115Combined.rds")

e95Combined$eday = "E9.5"
e105Combined$eday = "E10.5"
e115Combined$eday = "E11.5"

allCombined = merge(e95Combined, y = list(e105Combined, e115Combined))
DefaultAssay(allCombined) = "RNA"
allCombined = NormalizeData(allCombined)
allCombined = JoinLayers(allCombined)
#cell cycle
# A list of HUMAN cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
source('~/GitHub/NeuralCrest_singlecell/functions/ConvertHumanGeneListtoMM.R')

s.genes_mm = ConvertHumanGeneListToMM(s.genes)
g2m.genes_mm = ConvertHumanGeneListToMM(g2m.genes)

#calculate cell cycle phase score

allCombined = CellCycleScoring(object = allCombined, s.features = s.genes_mm, g2m.features = g2m.genes_mm)
allCombined$CC.Difference = allCombined$S.Score - allCombined$G2M.Score
colnames(allCombined@meta.data)


allCombined <- ScaleData(allCombined, assay = "RNA", vars.to.regress = c("nCount_RNA", "CC.Difference"), 
                         features = rownames(allCombined[["RNA"]]))
allCombined = FindVariableFeatures(allCombined, assay = "RNA")

allCombined = RunPCA(allCombined, assay = "RNA")


allCombined = RunHarmony(allCombined, group.by.vars = "eday")
allCombined = RunUMAP(allCombined, reduction = "harmony", dims = 1:20)
allCombined <- FindNeighbors(allCombined, reduction = "harmony", dims = 1:20)
names(allCombined@graphs)

allCombined = FindClusters(allCombined, resolution = 0.5) #change later
#setwd("~/BINF/yushi scrnaseq/all six/threshold0/all_trial")
plot1 = DimPlot(allCombined, reduction = "umap", group.by = "eday")
plot2 = DimPlot(allCombined, reduction = "umap", label = TRUE)
markers = FindAllMarkers(allCombined, only.pos = TRUE)
write_xlsx(markers, "markers_res.0.5v.xlsx")
saveRDS(allCombined,"allCombined.rds")

############
#further testing harmony
DimHeatmap(object = allCombined, reduction = "harmony", cells = 500, dims = 1:4)

###################
#monocle


##################
#7.28.2025
#monocle with 2 partition: Sox9 and PAx3
#constrain the trajectory to follow the biological order
##################

allCombined <- readRDS("C:/Users/neha/Documents/BINF/yushi scrnaseq/all six/threshold0/harmony/monocle/allcombined_with_harmony.rds")

colnames(allCombined@meta.data)
#no cell cycle scoring

#create new cds 

cds = as.cell_data_set(allCombined)

cds <- preprocess_cds(cds, num_dim = 50) #default uses PCA
cds <- reduce_dimension(cds)  # by default this uses UMAP
cds <- cluster_cells(cds)     
#cds <- learn_graph(cds, use_partition = TRUE)




'
#set partitions = orig.ident
colData(cds)$orig.ident <- allCombined$orig.ident
cds@clusters$UMAP$partitions <- factor(colData(cds)$orig.ident, levels = c("Sox9", "Pax3"))

head(cds@clusters$UMAP$partitions)
'




head(colData(cds)$timepoint)
colData(cds)$timepoint <- allCombined$timepoint
#convert to ordered factors to timepoint = directional pseudotime
colData(cds)$timepoint <- factor(colData(cds)$timepoint,
                                 levels = c("E9.5", "E10.5", "E11.5"),
                                 ordered = TRUE)


cds <- learn_graph(cds, use_partition = FALSE)


root_cells <- colnames(cds)[colData(cds)$timepoint == "E9.5"]

# Order cells in pseudotime
cds <- order_cells(cds, root_cells = root_cells)


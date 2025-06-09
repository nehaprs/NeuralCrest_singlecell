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

library(harmony)
library(slingshot)
library(SingleCellExperiment)
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
#=====================================================

#code similar to top portion of integration.R


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
allCombined <- ScaleData(allCombined, assay = "RNA", features = rownames(allCombined[["RNA"]]))
allCombined = FindVariableFeatures(allCombined, assay = "RNA")

allCombined = RunPCA(allCombined, assay = "RNA")


allCombined = RunHarmony(allCombined, group.by.vars = "eday")
allCombined = RunUMAP(allCombined, reduction = "harmony", dims = 1:20)
allCombined <- FindNeighbors(allCombined, reduction = "harmony", dims = 1:20)
names(allCombined@graphs)
allCombined = JoinLayers(allCombined)
allCombined = FindClusters(allCombined, resolution = 0.5) #change later
setwd("~/BINF/yushi scrnaseq/all six/threshold0/all_trial")
plot1 = DimPlot(allCombined, reduction = "umap", group.by = "eday")
plot2 = DimPlot(allCombined, reduction = "umap", label = TRUE)
markers = FindAllMarkers(allCombined, only.pos = TRUE)
write_xlsx(markers, "markers_res.0.5v2.xlsx")
saveRDS(allCombined,"allCombined.rds")

############
#further testing harmony
DimHeatmap(object = allCombined, reduction = "harmony", cells = 500, dims = 1:4)

###################
#monocle

#ignore everything below and try with perplexity.
library(stringr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(monocle3)
library(stringr)
library(SeuratWrappers) 

cds = as.cell_data_set(allCombined)
head(colData(cds))
rownames(fData(cds))[1:10]
fData(cds)$gene_short_name <- rownames(fData(cds))

recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

list.cluster <- allCombined@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- allCombined@reductions$umap@cell.embeddings
cluster.before.traj <-plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F, 
                                 group_label_size = 5) + theme(legend.position = "right")

cds <- learn_graph(cds, use_partition = F)
plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5)


root_cells = WhichCells(allCombined,expression = eday == "E9.5")
e95_root_cells = colnames(cds)[cds$eday == "E9.5"]
#root_clusters = unique(allCombined$eday == "E9.5")
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = e95_root_cells)
#plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F)



plot_cells(cds, 
           color_cells_by="pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE) +
  facet_wrap(~eday)

ggplot(colData(cds), aes(x=pseudotime, fill=eday)) + 
  geom_density(alpha=0.5) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#009E73"))

head(pseudotime(cds), 10)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))


ggplot(data.pseudo, aes(monocle3_pseudotime, seurat_clusters, fill = seurat_clusters)) + geom_boxplot()
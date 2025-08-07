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
setwd("~/BINF/yushi scrnaseq/all six/threshold0/harmony/monocle/choose_Roots")
allCombined <- readRDS("C:/Users/neha/Documents/BINF/yushi scrnaseq/all six/threshold0/harmony/monocle/allcombined_with_harmony.rds")
DimPlot(allCombined, group.by = "predicted.id")

#find nc cells
target_cells = WhichCells(allCombined, expression = predicted.id %in% 
                            c("Neural crest (PNS glia)", "Neural Crest (PNS Neuron)"))
# Add a new metadata column 'nc.cells', marking only the selected cells as TRUE
allCombined$nc.cells = FALSE
allCombined$nc.cells[target_cells] = TRUE
length(target_cells) #there are 1767 nc cells
ncol(allCombined) #there are 16157 cells in allCombined


# Create logical vector for the root cells

root.cells = allCombined$nc.cells == TRUE & allCombined$timepoint == "E9.5"
allCombined$root.cells = root.cells
sum(allCombined$root.cells == TRUE) #275 cells

head(allCombined$nc.cells)

colnames(allCombined@meta.data)
#no cell cycle scoring

#create new cds 

cds = as.cell_data_set(allCombined)

table(colData(cds)$root.cells) #275 root cells

cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds)
cds <- learn_graph(cds,use_partition = FALSE)
root_cells <- rownames(subset(colData(cds), root.cells == TRUE))
cds <- order_cells(cds, root_cells = root_cells)
plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE, label_principal_points = FALSE
      )
plot_cells(cds, color_cells_by = "timepoint", show_trajectory_graph = TRUE, label_principal_points = FALSE,
           label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE, label_roots = FALSE)
pseudotime_df <- data.frame(cell_id = colnames(cds), celltype = cds$predicted.id,
                            pseudotime = pseudotime(cds))

table(cds$ident)


#compute average pseudotime per cluster

pseudotime_df <- data.frame(
  cluster = colData(cds)$seurat_clusters,
  pseudotime = pseudotime(cds)
)

avg_pt <- pseudotime_df %>%
  group_by(cluster) %>%
  summarize(mean_pseudotime = mean(pseudotime, na.rm = TRUE)) %>%
  arrange(mean_pseudotime)



type_eday_df <- data.frame(
  cluster = colData(cds)$seurat_clusters,
  pseudotime = cds$timepoint,
  origin = cds$orig.ident
)

head(type_eday_df)

# 
cluster_summary <- df %>%
  group_by(cluster) %>%
  summarize(
    most_common_pseudotime = names(sort(table(pseudotime), decreasing = TRUE))[1],
    most_common_origin = names(sort(table(origin), decreasing = TRUE))[1]
  ) %>%
  ungroup()


type_eday_summary = type_eday_df %>%
  group_by(cluster) %>%
  summarize(timepoint = names(sort(table(pseudotime), decreasing = TRUE))[1],
            origin = names(sort(table(origin), decreasing = TRUE))[1]
            )%>%
  ungroup()
  
sort(table(type_eday_df$pseudotime), decreasing = TRUE)

write_xlsx(avg_pt, "avg_pseudotime_clusters.xlsx")
write_xlsx(type_eday_summary,"type_eday_cluster.xlsx")
saveRDS(cds, "cds.rds")

pr_graph_test_res <- principal_graph(cds)[["UMAP"]]



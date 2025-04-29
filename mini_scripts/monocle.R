library(stringr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(monocle3)
library(stringr)
library(SeuratWrappers) 


setwd("~/BINF/yushi scrnaseq/time series/harmony_slingshot/paxfull")
s.processed = readRDS("paxFullCombined_procesd.rds")


s.processed[["RNA"]] = JoinLayers(s.processed[["RNA"]])

# Find root cells: E9.5 neural crest
table(s.processed@meta.data$cell_state)
root_cells = WhichCells(s.processed,
                        expression = eday == "E9.5" & str_detect(s.processed@meta.data$predicted.id, regex("Neural crest", ignore_case = TRUE)))


root_clusters = unique(s.processed$cell_state[root_cells])
cds = as.cell_data_set(s.processed)
cds <- cluster_cells(cds, reduction_method = "UMAP")
colData(cds)$cluster <-s.processed$cell_state
cds <- learn_graph(cds, use_partition = FALSE)
cds <- order_cells(cds, root_cells = root_cells)
plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_roots = FALSE,
  label_branch_points = FALSE,
  label_leaves  = FALSE,
  label_cell_groups = FALSE
)

pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

s.processed$pseudotime <- pseudotime(cds)
Idents(s.processed) <-s.processed$pseudotime
FeaturePlot(s.processed, features = "pseudotime", label = F)


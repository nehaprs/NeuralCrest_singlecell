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
Idents(s.processed)= s.processed@meta.data$cell_state
DimPlot(s.processed) +NoLegend()
#pre-process and do partitions
cds = as.cell_data_set(s.processed)

# to get cell metadata
colData(cds)
# to gene metdata
fData(cds)
rownames(fData(cds))[1:10]

# since it misses the gene_short_name column, let's add it
fData(cds)$gene_short_name <- rownames(fData(cds))


# to get counts
counts(cds)


# Assign the cluster info 

list_cluster <- s.processed@active.ident
cds@clusters$UMAP$clusters <- list_cluster


# Assign UMAP coordinate - cell embeddings

cds@int_colData@listData$reducedDims$UMAP <-  s.processed@reductions$umap@cell.embeddings


# plot

cluster.before.trajectory <- plot_cells(cds,
                                        color_cells_by = 'cluster',
                                        label_groups_by_cluster = FALSE) +
  NoLegend()



cluster.names <- plot_cells(cds,
                            color_cells_by = "cell_state",
                            label_groups_by_cluster = FALSE) +
  NoLegend()



# ...3. Learn trajectory graph ------------------------
cds <- learn_graph(cds, use_partition = FALSE)

plot_cells(cds,
           color_cells_by =  "cell_state",
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE)














##################################
#cds <- preprocess_cds(cds, num_dim = 30)
#cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds, reduction_method = "UMAP")

cds@clusters$UMAP$partitions

#find partition with nc cells


plot_cells(cds, color_cells_by = "partition", label_groups_by_cluster = TRUE)


cds <- learn_graph(cds, use_partition =FALSE)
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


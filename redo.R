library(dplyr)
library(Seurat)
library(patchwork)
library(writexl)
library(readxl)
library(clustree)
library(hdf5r)
library(scCATCH)
library(tidyr)  
library(harmony)
library(monocle3)
library(SeuratWrappers)

setwd("~/BINF/yushi scrnaseq/redo")
#sox9
sox95 <- readRDS("~/BINF/yushi scrnaseq/E9.5/Sox9/redo1125/sox9.rds")
sox105 <- readRDS("~/BINF/yushi scrnaseq/E10.5/sox9/seurat output/sox9.rds")
sox115 <- readRDS("~/BINF/yushi scrnaseq/E11.5/sox9/seurat output/sox9.rds")

#make cell names unique across samples
sox95 = RenameCells(sox95, add.cell.id = "s95")
sox105 = RenameCells(sox105, add.cell.id = "s105")
sox115 = RenameCells(sox115, add.cell.id = "s115")

#from clustree, fix a stable cluster
Idents(sox95) = sox95$RNA_snn_res.0.5 # 18 levels
Idents(sox105) = sox105$RNA_snn_res.0.5
Idents(sox115) = sox115$RNA_snn_res.0.6

seu_list = list(sox95, sox105, sox115)

features = SelectIntegrationFeatures(seu_list, nfeatures = 3000) #2000 is default, but 3k is recommended.

seu_list <- lapply(seu_list, \(x) ScaleData(x, features = features))
seu_list <- lapply(seu_list, \(x) RunPCA(x, features = features))

I2 = merge(seu_list[[1]],  c(seu_list[[2]], seu_list[[3]]))
saveRDS(I2,"mergedSox.rds")
head(colnames(I2))
head(I2@meta.data$orig.ident)
cn <- colnames(I2)
first_char <- substring(cn, 1, 3)
table(first_char)

I2$batch = first_char

I2 = NormalizeData(I2)
I2 = FindVariableFeatures(I2, nfeatures = 3000)
I2 = ScaleData(I2)
I2 = RunPCA(I2, npcs = 30)

elbow = ElbowPlot(I2) #9

I2 = FindNeighbors(I2, dims = 1:9)
I2 = RunUMAP(I2, dims = 1:9)
I2 = FindClusters(I2)

DimPlot(I2, group.by = "batch")

#create ori_cluster
I2$orig_cluster = paste0(I2$batch, "_", I2$seurat_clusters)

saveRDS(I2,"mergNoHarm.rds")


#quantify batch effect

library(lisi)
emb <- Embeddings(I2, "pca")[,1:30]
lisi <- compute_lisi(emb, I2@meta.data, "batch")
mean_iLISI <- mean(lisi$iLISI, na.rm=TRUE)
mean_iLISI

#################
#monocle
###############



m = I2
m[["RNA"]] <- JoinLayers(m[["RNA"]])
#Layers(m)
cds = as.cell_data_set(m)

colData(cds)$orig_cluster = m$orig_cluster
colData(cds)$timepoint = m$batch

reducedDims(cds)$UMAP <- Embeddings(m, "umap")
reducedDims(cds)$PCA  <- Embeddings(m, "pca")


cds = cluster_cells(cds, reduction_method = "UMAP")
cds = learn_graph(cds, use_partition = TRUE)

orig_cluster <- colData(cds)$orig_cluster
root_cells = rownames(colData(cds))[orig_cluster == "s95_4"]
cds <- order_cells(cds, root_cells = root_cells)
# sanity plots
bypt = plot_cells(cds, color_cells_by="pseudotime", label_groups_by_cluster=FALSE)

table(m$orig_cluster)

#assign branch lineages of s1_10
library(igraph)

# principal graph on UMAP
g <- principal_graph(cds)$UMAP 

class(g) #makes sure g is an igraph object
head(V(g))
vn <- V(g)$name
comp <- components(g)  




# matrix 'closest' linking nearest graph vertex for each cell
#rows(closest)== cells
#clms(closest)==nearest vertex
#principal_graph_aux: list of auxillary data encapsulating Spatial relationship between cells and the graph
closest <- cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex
closest <- as.matrix(closest[colnames(cds), , drop=FALSE])
head(closest)
# Replace numeric indices of closest with the corresponding vertex names
if (is.numeric(closest[,1])) {
  closest[,1] <- vn[closest[,1]]
}
head(closest)



#Create a named vector mapping each cell barcode to its closest graph vertex
vmap <- setNames(closest[,1], rownames(closest))
head(vmap)



# root vertex:retrieve the nearest vertex (node) for root cells.
root_vs <- vmap[root_cells]
'
table(root_vs)
root_vs
  1   2   3   4   5  10  12  15  16  17  20  22  24  25  27  28  31  36 
  3 267   6   2   1   3  26   8   6  32   2 239  11 131   1   1   4  12

'



#select the most common graph vertex among the root cells
root_v  <- names(sort(table(root_vs), decreasing=TRUE))[1]
#root_v = Y_2


# leaves and partitions
'
definitions:
degree(v)=number of edges incident to v
Interpretation in Monocle3

In the trajectory graph:

Degree = 1 → a leaf node (end of a lineage).

Degree = 2 → an intermediate node along a continuous branch.

Degree ≥ 3 → a branch point (where a lineage splits).

'
deg <- degree(g)
leaf_vs <- names(deg[deg==1])
head(leaf_vs)
# precompute shortest paths from root to each leaf
paths <- lapply(leaf_vs, function(L) shortest_paths(g, from=root_v, to=L)$vpath[[1]])
#convert each graph path into a vector of vertex names
#each element of path_set = character vector of node names along one path
path_sets <- lapply(paths, function(p) names(p))
# assign each vertex to the furthest leaf whose path contains it
# vertices shared by multiple leaf paths are "pre-bifurcation"
assign_leaf <- function(v){
  idx <- which(vapply(path_sets, function(S) v %in% S, logical(1)))
  if(length(idx)==0) return(NA_character_)
  if(length(idx)==1) return(leaf_vs[idx])
  return("pre_bifurcation")
}
vertex_to_leaf <- setNames(vapply(V(g)$name, assign_leaf, character(1)), V(g)$name)




# cell lineage label
cell_lineage <- setNames(vertex_to_leaf[vmap], names(vmap))
colData(cds)$lineage_leaf <- cell_lineage
head(colData(cds)$lineage_leaf)



# restrict to descendants reachable from cluster 10 root
descendant_cells <- names(cell_lineage)[!is.na(cell_lineage)]
desc_cds <- cds[, descendant_cells]

#check what label was used for time point
table(cds@colData$batch)
table(cds@colData$orig_cluster)


lin_flow <- as.data.frame(colData(desc_cds)) |>
  mutate(leaf = lineage_leaf,
         tp = batch,
         is_root = (orig_cluster == "s95_4")) |>
  group_by(tp, leaf) |>
  summarise(n=n(), .groups="drop") |>
  arrange(tp, desc(n))
lin_flow
#lineages in lin_flow are Y_1, Y_43, Y_5
write_xlsx(lin_flow,"lin_flow.xlsx")
ggplot(lin_flow, aes(tp, n, fill = leaf)) + geom_bar(stat = "identity", position = "fill")

library(ggalluvial)
df <- as.data.frame(colData(desc_cds)) |>
  select(batch, orig_cluster, lineage_leaf) |>
  mutate(n=1) |>
  group_by(batch, orig_cluster, lineage_leaf) |>
  summarise(n=n(), .groups="drop")


ggplot(df,
       aes(axis1=batch, axis2=orig_cluster, axis3=lineage_leaf, y=n)) +
  geom_alluvium() + geom_stratum() + geom_text(stat="stratum", infer.label=TRUE) +
  scale_x_discrete(limits=c("batch","cluster","leaf"), expand=c(.1,.1)) +
  ylab("cells")
table(df)

df1 <- df |>
  dplyr::group_by(orig_cluster) |>
  dplyr::summarise(n = sum(n), .groups = "drop") |>
  dplyr::arrange(desc(n))

df1 = df %>%
  group_by(orig_cluster) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  arrange(desc(n))

write_xlsx(df1,"desc_clustersR1.xlsx")

###checking why we see biologically impossible lineages
igraph::components(principal_graph(cds)$UMAP)$no


leaf_tab <- as.data.frame(colData(desc_cds)) |>
  filter(batch== 4) |>
  dplyr::count(lineage_leaf, orig_cluster, name="n") |>
  group_by(lineage_leaf) |>
  mutate(frac = n/sum(n)) |>
  arrange(lineage_leaf, desc(frac))
leaf_tab


names(colData(desc_cds))



end_alloc <- as.data.frame(colData(desc_cds)) |>
  group_by(batch, lineage_leaf) |>
  summarise(cells=n(), .groups="drop")
write.csv(end_alloc, "descendants_endpoint_allocation.csv", row.names=FALSE)



























#plot cells
plot_cells(cds, color_cells_by="lineage_leaf", label_groups_by_cluster=FALSE,
           label_leaves=TRUE, label_branch_points=TRUE)

plot_cells(cds, color_cells_by="orig_cluster", label_groups_by_cluster=FALSE,
           label_leaves=, label_branch_points=FALSE)


#highlight root cells

# 1. Flag root cells
colData(cds)$is_root <- ifelse(colnames(cds) %in% root_cells, "root", "other")

# 2. Plot with custom colors
plot_cells(
  cds,
  color_cells_by = "is_root",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  show_trajectory_graph = TRUE
) + scale_color_manual(values = c("root" = "red", "other" = "grey"))






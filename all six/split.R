#######################
#8.13.2025: split all six into sox and pax objects. pseudotime on that.
######################

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

setwd("~/BINF/yushi scrnaseq/all six/threshold0/split")

s.obj <- readRDS("~/BINF/yushi scrnaseq/all six/threshold0/harmony/monocle/allcombined_with_harmony.rds")

#split into sox and pax cells

sox9 = subset(s.obj, subset = orig.ident == "sox9")
pax3 = subset(s.obj, subset = orig.ident == "pax3")

####################
# Run the standard workflow for visualization and clustering
s.combined = sox9

s.combined <- ScaleData(s.combined, verbose = FALSE)
s.combined <- RunPCA(s.combined, npcs = 30, verbose = FALSE)
elbow = ElbowPlot(s.combined) #16
s.combined <- RunUMAP(s.combined, reduction = "pca", dims = 1:16)
s.combined <- FindNeighbors(s.combined, reduction = "pca", dims = 1:16)
s.combined <- FindClusters(s.combined, resolution = 0.5)

p1 <- DimPlot(s.combined, reduction = "umap", group.by = "seurat_clusters")+ggtitle("Combined Dataset all Sox9 Cells")


#use 
resolution.range <- seq(from = 0, to = 2, by = 0.1)

# Loop over each resolution
for (res in resolution.range) {
  # Perform clustering with the current resolution
  s.combined<- FindClusters(s.combined, resolution = res)
  
  # Find all markers for the clusters at this resolution
  s.combined.markers <- FindAllMarkers(s.combined, only.pos = TRUE)
  
  # Define the file name for saving the markers
  file_name <- paste0("markers_resolution_", res, ".xlsx")
  
  # Save the markers as an Excel file
  write_xlsx(s.combined.markers, file_name)
  
  # Print a message to confirm completion for each resolution
  print(paste("Markers for resolution", res, "saved to", file_name))
}


#list all xlsx files in wd

xlsx_file = list.files(pattern = "\\.xlsx$")

for (file in xlsx_file){
  df = read_xlsx(file)
  dff = df[df$avg_log2FC > 0,]
  dfff = dff[dff$p_val_adj < 0.05,]
  file_new = paste0("filt",file)
  write_xlsx(dfff, file_new)
}

s.combinedclust = clustree(s.combined)

#number of cells in each cluster at desired res
s.combined$seurat_clusters = s.combined$RNA_snn_res.0.6
table(Idents(s.combined))
Idents(s.combined) = s.combined$RNA_snn_res.0.6
table(Idents(s.combined))

saveRDS(s.combined,"sox9Combined.rds")


###################
#monocle
##################

#root cells: E9.5 cells



root.cells =  s.combined$timepoint == "E9.5" 
s.combined$root.cells = root.cells
sum(s.combined$root.cells == TRUE) #3232 E9.5 nc cells

cds = as.cell_data_set(s.combined)
table(colData(cds)$root.cells) #3232 root cells

cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds)
cds <- learn_graph(cds,use_partition = FALSE)
root_cells <- rownames(subset(colData(cds), root.cells == TRUE))
cds <- order_cells(cds, root_cells = root_cells)
plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE, label_principal_points = FALSE,
           label_groups_by_cluster = FALSE, label_leaves = FALSE,  label_roots = FALSE
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
library(dplyr)
library(Seurat)
library(patchwork)
library(writexl)
library(readxl)
library(clustree)
library(tidyr)
library(harmony)
library(slingshot)
library(SingleCellExperiment)


#=====================================
#integration without harmony
#=====================================
setwd("~/BINF/yushi scrnaseq/all six/threshold0/noharmony")

e95Combined <- readRDS("~/BINF/yushi scrnaseq/all six/threshold0/e95/e95Combined.rds")
e105Combined <- readRDS("~/BINF/yushi scrnaseq/all six/threshold0/e105/e105Combined.rds")
e115Combined <- readRDS("~/BINF/yushi scrnaseq/all six/threshold0/e115/e115Combined.rds")

e95Combined$timepoint = "E9.5"
e105Combined$timepoint = "E10.5"
e115Combined$timepoint = "E11.5"

#prefix cell bar codes so they stay unique after merging
e95Combined  = RenameCells(e95Combined,   new.names = paste0("E9.5_",  colnames(e95Combined)))
e105Combined = RenameCells(e105Combined, new.names = paste0("E10.5_", colnames(e105Combined)))
e115Combined = RenameCells(e115Combined, new.names = paste0("E11.5_", colnames(e115Combined)))

#Merge them into one object
allcombined <- merge(
  x = e95Combined,
  y = c(e105Combined, e115Combined),
  add.cell.ids = c("E9.5", "E10.5", "E11.5")  # this prefixes cell IDs under the hood
)

table(allcombined$timepoint)
'
E10.5 E11.5  E9.5 
 4797  6704  4656
'

#preprocess the combined object
DefaultAssay(allcombined) = "RNA"
allcombined = JoinLayers(allcombined, assay = "RNA")
allcombined = NormalizeData(allcombined)

allcombined = FindVariableFeatures(allcombined)
all.genes = rownames(allcombined)

allcombined = ScaleData(allcombined, features = all.genes)
allcombined = RunPCA(allcombined, features = VariableFeatures(allcombined))
elbow = ElbowPlot(allcombined) #10
allcombined = RunUMAP(allcombined, dims = 1:10)
allcombined = FindNeighbors(allcombined, dims = 1:10)

resolution.range <- seq(from = 0, to = 2, by = 0.1)

for (res in resolution.range) {
  # Perform clustering with the current resolution
  allcombined<- FindClusters(allcombined, resolution = res)
  
  # Find all markers for the clusters at this resolution
  #all.markers <- FindAllMarkers(allcombined, only.pos = TRUE)
  
  # Define the file name for saving the markers
  #file_name <- paste0("markers_resolution_", res, ".xlsx")
  
  # Save the markers as an Excel file
  #write_xlsx(all.markers, file_name)
  
  # Print a message to confirm completion for each resolution
  #print(paste("Markers for resolution", res, "saved to", file_name))
}

clus = clustree(allcombined)
#choosw res 2 for now. later try larger/smaller clusters

DimPlot(allcombined, group.by = "seurat_clusters", label = TRUE,
        repel = TRUE) + ggtitle("Clusters, res 0.5")
DimPlot(allcombined, group.by = "timepoint") + ggtitle("Timepoints")


#####################################################################
#try res 0.5

allcombined$seurat_clusters = allcombined$RNA_snn_res.0.5
setwd("~/BINF/yushi scrnaseq/all six/threshold0/noharmony/res 0.5")
######################################################################

#run slingshot

#extract umap embeddings and cluster labels

umap_embedding = Embeddings(allcombined,"umap")
clusterlabels = allcombined$seurat_clusters #res 0.5

#id E9.5 clusters for starting clusters
cluster_time_df = as.data.frame(table(clusterlabels, allcombined$timepoint)) 
colnames(cluster_time_df) <- c("cluster", "timepoint", "count")
write_xlsx(cluster_time_df,"cluster_time.xlsx")

#find fraction of E9.5 
frac_E9.5 = cluster_time_df %>%
  filter(timepoint == "E9.5") %>%
  #rename(E9.5_count = count) %>%
  select(cluster, count) %>%
  left_join(
    cluster_time_df %>% 
      group_by(cluster) %>% 
      summarize(total = sum(count)), 
    by = "cluster"
  ) %>%
  mutate(frac_E9.5 = count / total)

write_xlsx(frac_E9.5, "fraction_e9.5.xlsx")

start_clusters <- as.character(frac_E9.5 %>% filter(frac_E9.5 > 0.8) %>% pull(cluster))
start_clusters
'
"6"  "8"  "10" "15" "16" "18" "19" "22" "28" "35" "36" "37" #for res 2
'

cluster_factor = as.factor(clusterlabels)

#run slingshot
sds <- slingshot(
  umap_embedding,
  clusterLabels = cluster_factor,
  start.clus = start_clusters)
  
  'shrink = TRUE,             # shrinks smooth curves toward the spanning tree
  reweight = TRUE            # reweights outliers during curve fitting
)
'
  
#inspect the lineages found
n_lineages = length(slingLineages(sds))
# 20 lineages found
all_lineages = slingLineages(sds)


cluster_time_stats <- cluster_time_df %>%
  group_by(cluster) %>%
  mutate(total_cells = sum(count)) %>%
  ungroup() %>%
  pivot_wider(names_from = timepoint, values_from = count, values_fill = 0L) %>%
  mutate(
    frac_E9.5  = E9.5  / total_cells,
    frac_E10.5 = E10.5 / total_cells,
    frac_E11.5 = E11.5 / total_cells
  )




#  Function to label each cluster by its dominant timepoint
label_by_time <- function(cluster_id) {
  row <- cluster_time_stats %>% filter(cluster == cluster_id)
  # Find which timepoint fraction is highest
  times <- c("E9.5", "E10.5", "E11.5")
  fracs <- c(row$frac_E9.5, row$frac_E10.5, row$frac_E11.5)
  tp  <- times[which.max(fracs)]
  pct <- round(max(fracs) * 100, 1)
  paste0("Cluster ", cluster_id, " (", tp, " – ", pct, "%)")
}

#  For each lineage, produce a string
for (ln in names(all_lineages)) {
  clusters_in_lineage <- all_lineages[[ln]]
  labels <- vapply(clusters_in_lineage, label_by_time, FUN.VALUE = character(1))
  cat(ln, ": ", paste(labels, collapse = " → "), "\n\n")
}



#visualize

library(RColorBrewer)
library(slingshot)
library(ggplot2)

# (5a) Extract Slingshot curves in UMAP space
curve_data <- slingCurves(sds)  
# slingCurves(sds) is a list of SlingshotCurve objects. Each has a slot $s (the smoothed coordinates).

# (5b) Plot UMAP + curves
umap_df <- as.data.frame(umap_embedding)
umap_df$cluster <- clusterlabels
umap_df$timepoint <- allcombined$timepoint

# Basic UMAP scatter
p <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = timepoint)) +
  geom_point(size = 0.5, alpha = 0.6) +
  theme_classic() +
  ggtitle("UMAP colored by Timepoint, with Slingshot curves")

# Overlay each lineage curve
for (i in seq_along(curve_data)) {
  curve_coords <- curve_data[[i]]$s
  curve_df <- as.data.frame(curve_coords)
  colnames(curve_df) <- c("umap_1", "umap_2")
  p <- p + geom_path(data = curve_df, aes(x = umap_1, y = umap_2),
                     color = "black", size = 1.2, linetype = "solid")
}

print(p)

#troubleshooting
setdiff(start_clusters, levels(cluster_factor))
# A data frame of each cluster center's coordinates
centers_df <- slingClusterLabels(sds)  
# The MST edge list (pairs of cluster IDs connected)
mst_edges  <- slingMST(sds)
print(mst_edges)

#save stuff
setwd("~/BINF/yushi scrnaseq/all six/threshold0/noharmony")
saveRDS(allcombined,"combined_seurat.rds")
saveRDS(sds,"combined_slingshot_obj.rds")
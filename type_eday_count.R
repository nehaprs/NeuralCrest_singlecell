#############################
#for each cluster, counts the % of cells from each time point and each original identity
#############################


meta <- as.data.frame(colData(cds)) %>%
  select(seurat_clusters, timepoint, orig.ident)

cluster_sizes  = as.data.frame(table(as.numeric(meta$seurat_clusters)))
colnames(cluster_sizes) = c("seurat_cluster", "total_cells")

meta <- as.data.frame(colData(cds))

library(dplyr)
library(tidyr)

# Assuming `colData(cds)` has 'cluster', 'timepoint', and 'orig.ident'
meta <- as.data.frame(colData(cds)) %>%
  select(seurat_clusters, timepoint, orig.ident)

# Count number of cells per cluster
cluster_sizes <- meta %>%
  dplyr::count(seurat_clusters, name = "total_cells")

# Count each timepoint per cluster
timepoint_counts <- meta %>%
  dplyr::count(seurat_clusters, timepoint) %>%
  left_join(cluster_sizes, by = "seurat_clusters") %>%
  mutate(percent = 100 * n / total_cells) %>%
  select(seurat_clusters, timepoint, percent) %>%
  pivot_wider(names_from = timepoint, values_from = percent, values_fill = 0)

# Count each orig.ident per cluster
orig_counts <- meta %>%
  dplyr::count(seurat_clusters, orig.ident) %>%
  left_join(cluster_sizes, by = "seurat_clusters") %>%
  mutate(percent = 100 * n / total_cells) %>%
  select(seurat_clusters, orig.ident, percent) %>%
  pivot_wider(names_from = orig.ident, values_from = percent, values_fill = 0)

# Join the two summaries
final_table <- left_join(timepoint_counts, orig_counts, by = "seurat_clusters")

write_xlsx(final_table, "type_eday_list.xlsx")

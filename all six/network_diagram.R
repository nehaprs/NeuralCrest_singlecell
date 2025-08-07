#########################
#network diagram showing parent child relationships
#########################

library(monocle3)
library(dplyr)
library(igraph)
library(ggraph)
library(ggplot2)


pseudotime = pseudotime(cds)
colData(cds)$pseudotime = pseudotime
#Extract pseudotime and cluster metadata
meta <- as.data.frame(colData(cds))
meta$cell <- rownames(meta)
meta <- meta %>%
  select(cell, pseudotime = pseudotime, cluster = seurat_clusters)

#compute average pseudotime for each cluster
cluster_pt <- meta %>%
  group_by(cluster) %>%
  summarise(avg_pseudotime = mean(pseudotime, na.rm = TRUE)) %>%
  arrange(avg_pseudotime) %>%
  mutate(cluster = as.character(cluster))

#Define edges based on pseudotime ordering
edges <- data.frame(
  from = head(cluster_pt$cluster, -1),
  to   = tail(cluster_pt$cluster, -1)
)

# Create graph and plot
g <- graph_from_data_frame(edges, vertices = cluster_pt, directed = TRUE)

p <- ggraph(g, layout = 'sugiyama') +
  geom_edge_link(arrow = arrow(length = unit(4, 'mm')), end_cap = circle(3, 'mm')) +
  geom_node_circle(aes(r = 0.5), fill = 'lightblue') +
  geom_node_text(aes(label = name), vjust = -1.2) +
  theme_void() +
  ggtitle("Cluster Transitions by Pseudotime")

print(p)

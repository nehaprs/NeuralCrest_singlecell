#slingshot

library(stringr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(slingshot)
library(SingleCellExperiment)
setwd("~/BINF/yushi scrnaseq/time series/harmony_slingshot/paxfull")
s.processed = readRDS("paxFullCombined_procesd.rds")


s.processed[["RNA"]] = JoinLayers(s.processed[["RNA"]])



# Find root cells: E9.5 neural crest
table(s.processed@meta.data$cell_state)
root_cells = WhichCells(s.processed,
                        expression = eday == "E9.5" & str_detect(s.processed@meta.data$predicted.id, regex("Neural crest", ignore_case = TRUE)))


root_clusters = unique(s.processed$cell_state[root_cells])

#sce object and slingshot
sce = as.SingleCellExperiment(s.processed)

colData(sce)$cluster = s.processed$cell_state







#remove smaller clusters etc


umap_emb <- reducedDim(sce, "UMAP")
sum(is.na(umap_emb))
sum(is.infinite(umap_emb))

umap_emb = reducedDim(sce,"UMAP")
sum(is.na(umap_emb)) #0
sum(is.infinite(umap_emb)) #0
#found none

#check for none or inf in clusters

table(is.na(sce$cluster))
#FALSE



#error might be because some starting clusters are too small

cluster_sizes = table(sce$cluster)
small_clusters = names(cluster_sizes[cluster_sizes < 3])
small_clusters
write.table(small_clusters,"small_clustersv2.txt")
#remove small clusters
valid_cells = !(sce$cluster %in% small_clusters) 
sce = sce[, valid_cells]
head(table(sce2$cluster))



sce2 = slingshot(sce2,
                clusterLabels = 'cluster',
                reducedDim = "UMAP", 
                start.clus = root_clusters,
                allow.breaks = FALSE,
                approx_points = 18717)
lineages = slingLineages(sce)



#set initial clusters: 

clusters = s.processed$cell_state
sce$clusters = clusters

##
eday95_clusters = unique(sce$clusters[s.processed$eday == "E9.5"])

# Run Slingshot using UMAP embeddings 
sce <- slingshot(sce, clusterLabels = 'clusters', reducedDim = 'UMAP', start.clus = eday95_clusters)

#solvingerror: nan values 
## Check for NA or infinite values in umap embeddings
umap_emb <- reducedDim(sce, "UMAP")
sum(is.na(umap_emb))
sum(is.infinite(umap_emb))

umap_emb = reducedDim(sce,"UMAP")
sum(is.na(umap_emb)) #0
sum(is.infinite(umap_emb)) #0
#found none

#check for none or inf in clusters
 
table(is.na(sce$clusters))
#FALSE
eday95_clusters
head(table(sce$clusters))

#error might be because some starting clusters are too small

cluster_sizes = table(sce$clusters)
small_clusters = names(cluster_sizes[cluster_sizes < 3])
small_clusters
write.table(small_clusters,"small_clusters.txt")



#remove small clusters
valid_cells = !(sce$clusters %in% small_clusters) 
sce2 = sce[, valid_cells]

# Recompute start.clus 
sce2$clusters <- factor(sce2$clusters)
eday95_clusters = unique(sce2$clusters[sce2$eday == "E9.5"])


#rerun slingshot


sce2 = slingshot(sce2, clusterLabels = "clusters",
                 reducedDim = "UMAP",
                 start.clus = eday95_clusters,
                 approx_points = 18717,
                 
                 allow.breaks = TRUE)

summary(slingCurves(sce2))
plot(reducedDims(sce2)$UMAP, col = as.numeric(sce2$clusters), pch = 16)
lines(SlingshotDataSet(sce2), lwd = 2, col = 'black')
lineages <- slingLineages(sce2)
##############################################################
library(RColorBrewer)

# Plot UMAP colored by pseudotime
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)

plot(reducedDims(sce2)$UMAP, col = colors[cut(slingPseudotime(sce2)[,1], breaks=100)],
     pch=16, asp = 1, cex = 0.7, main = "Slingshot Trajectory (Lineage 1)")
lines(SlingshotDataSet(sce2), lwd = 2, col = 'black')


# Example of plotting multiple lineages clearly:
lineages <- slingLineages(sce2)
pseudotime_L1 <- slingPseudotime(sce2)[, "Lineage1"]

library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(11, "Spectral"))(100)

plot(reducedDims(sce2)$UMAP,
     col = colors[cut(pseudotime_L1, breaks = 100)],
     pch = 16, asp = 1, main = "Pseudotime (Lineage1)")
lines(SlingshotDataSet(sce2), lwd = 2, col = 'black', type = 'lineages', lineage = 1)
saveRDS(sce2,"slingshotv1.rds")


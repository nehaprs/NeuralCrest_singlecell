library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
#refCombined <- readRDS("~/BINF/yushi scrnaseq/time series/harmony_slingshot/reference/refCombined.rds")
#s.Combined <- readRDS("~/BINF/yushi scrnaseq/time series/harmony_slingshot/soxfull/soxCombined.rds")

setwd("~/BINF/yushi scrnaseq/time series/harmony_slingshot/paxfull")
s.Combined = readRDS("paxCombined.rds")
s.Combined$cell_state = paste(s.Combined$eday, s.Combined$predicted.id,sep = ":")

tail(s.Combined@meta.data$cell_state)
DimPlot(s.Combined, group.by = "cell_state", label = TRUE, repel = TRUE, label.size = 4, pt.size = 1.5
        ) + NoLegend()

#color by seurat_clusters and label by cell_state

# Get UMAP embeddings
umap_coords = Embeddings(s.Combined, reduction = "umap") %>%
  as.data.frame() %>%
  mutate(
    
    ColorGroup = s.Combined@meta.data$predicted.id,
    LabelGroup = s.Combined@meta.data$cell_state
  )

# Compute the median coordinates for labels
label_coords = umap_coords %>% 
  group_by(LabelGroup) %>%
  summarize(umap_1 = median(umap_1), umap_2 = median(umap_2))

# Plot cells colored by one metadata column and labeled by another
p = ggplot(umap_coords, aes(x = umap_1, y = umap_2, color = ColorGroup)) +
  geom_point(size = 0.5, alpha = 0.8) +
  geom_text_repel(data = label_coords, aes(label = LabelGroup, size = 2),
  color = "black") +
  theme_void() +
    #guides(color = guide_legend(title = "Color Group")) +
  NoLegend() 

p

#change identity to cell_type

Idents(s.Combined)
Idents(s.Combined) = "cell_type"
#setwd("~/BINF/yushi scrnaseq/time series/harmony_slingshot/soxfull")
saveRDS(s.Combined, "paxFullCombined_procesd.rds")

library(Seurat)
#devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(ggplot2)
library(readxl)
library(dplyr)

pax95 <- readRDS("C:/Users/neha/Documents/BINF/yushi scrnaseq/E9.5/Pax3/ref_annot/pax95.rds")

DimPlot(pax95, group.by = "predicted.id", label = TRUE, repel = TRUE, label.size = 4, pt.size = 1.5) + NoLegend() + ggtitle("pax3+ Cells at E9.5")

Idents(pax95) <- "seurat_clusters"  # Set clusters as the identity
library(dplyr)

cluster_annotations <- pax95@meta.data %>%
  group_by(seurat_clusters) %>%
  count(predicted.id, sort = TRUE) %>%
  slice_max(n, n = 1) %>%
  ungroup() %>%
  select(seurat_clusters, predicted_cell_type = predicted.id)

new.cluster.ids <- c( "Skeletal muscle progenitors",
                     "Spinal cord (ventral)",
                     "Pre-epidermal keratinocytes", "Neural crest (PNS glia)",
                     "Mesencephalon/MHB",
                     "Extraembryonic mesoderm",
                     "Mesenchymal stromal cells",
                     "Neural crest (PNS glia)",
                     "Neuromesodermal progenitors","Gut and lung epithelium",
                     "Pre-epidermal keratinocytes",
                     "Neural crest (PNS glia)",
                     "Skeletal muscle progenitors",
                     "Neuron progenitor cells",
                     "First heart field","Endothelium"
)

names(new.cluster.ids) <- levels(pax95)          
pax95 <- RenameIdents(pax95, new.cluster.ids)
DimPlot(pax95, reduction = "umap", label = TRUE, pt.size = 0.6,repel = TRUE) + NoLegend() + ggtitle("Pax3+ Cells at E9.5")

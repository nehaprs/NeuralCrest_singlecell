library(Seurat)
#devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(ggplot2)
library(readxl)
library(dplyr)

#s.query <- readRDS("~/BINF/yushi scrnaseq/E9.5/Sox9/seurat output/round1/sox9_resolution_0.4.rds")
#s.ref <- readRDS("~/BINF/yushi scrnaseq/E9.5/tomeE9.5.rds")

setwd("~/BINF/yushi scrnaseq/E10.5/Pax3/ref_annot")
s.query <- readRDS("~/BINF/yushi scrnaseq/E10.5/Pax3/seurat output/pax3.rds")
s.ref <- readRDS("~/BINF/yushi scrnaseq/E10.5/tome_E10.5.rds")

#function to process seurat object
process = function(obj){
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)
  obj <- FindNeighbors(obj, dims = 1:30)
  obj <- FindClusters(obj)
  RunUMAP(obj, dims = 1:30)
  
  return(obj)
}

s.ref = process(s.ref)

#convert query feature names from gene symbol to ensembl

library(biomaRt)



#Set up the biomaRt connection to the Ensembl database for mouse
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

# Extract the row names (gene symbols) from the Seurat object
mouse_gene_symbols <- rownames(s.query)

# Map mouse gene symbols to Ensembl IDs using biomaRt
gene_mapping <- getBM(
  attributes = c("mgi_symbol", "ensembl_gene_id"),
  filters = "mgi_symbol",
  values = mouse_gene_symbols,
  mart = ensembl
)

# Match Ensembl IDs to the row names in the Seurat object
# Filter out rows without mapping
mapped_genes <- gene_mapping[match(mouse_gene_symbols, gene_mapping$mgi_symbol), ]
rownames(s.query) <- mapped_genes$ensembl_gene_id

# Verify the update
head(rownames(s.query))


#integrate with default common features
#common_features <- intersect(rownames(s.ref), rownames(s.query))
#common_features <- intersect(VariableFeatures(s.ref), VariableFeatures(s.query))

#length(common_features)
# Subset Seurat objects to shared features
#s.ref = subset(s.ref, features = common_features)
#s.query = subset(s.query, features = common_features)
'''
any(duplicated(rownames(s.query)))
duplicated_genes_query <- rownames(s.query)[duplicated(rownames(s.query))]
s.query <- s.query[!is.na(rownames(s.query)), ]
'''

s.anchors <- FindTransferAnchors(
  reference = s.ref,
  query = s.query,
  dims = 1:10,
  features = NULL
)


#s.anchors = FindTransferAnchors(reference = s.ref, query = s.query, dims = 1:15, features = common_features)
pedictions = TransferData(anchorset = s.anchors, refdata = s.ref$cell_type, dims = 1:10)

#Add transferred labels to s2 metadata

s.query = AddMetaData(s.query, metadata = pedictions)
#DimPlot(s2, group.by = "predicted.id", label = TRUE, label.size = 4)
DimPlot(s.query, group.by = "predicted.id", label = TRUE, repel = TRUE, label.size = 4, pt.size = 1.5) + NoLegend() + ggtitle("Pax3+ Cells at E10.5")









#######################
#pax95 <- readRDS("C:/Users/neha/Documents/BINF/yushi scrnaseq/E9.5/Pax3/ref_annot/pax95.rds")


Idents(s.query) <- "seurat_clusters"  # Set clusters as the identity
library(dplyr)

cluster_annotations <- s.query@meta.data %>%
  group_by(seurat_clusters) %>%
  count(predicted.id, sort = TRUE) %>%
  slice_max(n, n = 1) %>%
  ungroup() %>%
  dplyr::select(seurat_clusters, predicted_cell_type = predicted.id)

idnames = cluster_annotations$predicted_cell_type
new.cluster.ids = idnames


names(new.cluster.ids) <- levels(s.query)          
s.query <- RenameIdents(s.query, new.cluster.ids)
DimPlot(s.query, reduction = "umap", label = TRUE, pt.size = 0.6,repel = TRUE) + NoLegend() + ggtitle("Pax3+ Cells at E10.5")



saveRDS(s.query,"pax105.rds")
[]
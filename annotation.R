library(Seurat)

#library(SeuratData)
library(ggplot2)
library(readxl)
library(dplyr)
library(ggrepel)

#variable inputs
#setwd("~/BINF/yushi scrnaseq/E9.5/sox9/ref_annot")
#s.query <- readRDS("~/BINF/yushi scrnaseq/E9.5/Sox9/seurat output/round1/sox9_resolution_0.4.rds")
#s.ref <- readRDS("~/BINF/yushi scrnaseq/E9.5/Pax3/ref_annot/reference9.5.rds")
#setwd("~/BINF/yushi scrnaseq/E9.5/Pax3/ref_annot")
#s.query <- readRDS("~/BINF/yushi scrnaseq/E9.5/Pax3/ref_annot/pax3E9.5_resolution_0.6.rds")
#s.ref <- readRDS("~/BINF/yushi scrnaseq/E11.5/tomeRef_E11.5.rds")

s.query <- readRDS("~/BINF/yushi scrnaseq/time series/harmony_slingshot/soxCombined.rds")
s.ref <- readRDS("~/BINF/yushi scrnaseq/time series/harmony_slingshot/reference/refCombined_procesd.rds")
setwd("~/BINF/yushi scrnaseq/time series/harmony_slingshot/soxAnnot")

########
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
#s.query = cc
library(biomaRt)

##frog
source("~/GitHub/DorsalMigration/frog_nameconvert.R")
s.ref = frog_nameconvert(s.ref)

##mouse
#Set up the biomaRt connection to the Ensembl database for mouse
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
print("a")
# Extract the row names (gene symbols) from the Seurat object
mouse_gene_symbols <- rownames(s.query)

head(rownames(s.query))
# Map mouse gene symbols to Ensembl IDs using biomaRt
gene_mapping <- getBM(
  attributes = c("mgi_symbol", "ensembl_gene_id"),
  filters = "mgi_symbol",
  values = mouse_gene_symbols,
  mart = ensembl
)
print("a1")
#saveRDS(s.query,"soxFltd4bmx.rds")
setwd("~/BINF/yushi scrnaseq/time series/harmony_slingshot/soxAnnot")
#readRDS("soxFltd4bmx.rds")


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


any(duplicated(rownames(s.query)))
duplicated_genes_query <- rownames(s.query)[duplicated(rownames(s.query))]
s.query <- s.query[!is.na(rownames(s.query)), ]
any(duplicated(rownames(s.query)))


tail(rownames(s.ref))

s.anchors <- FindTransferAnchors(
  reference = s.ref,
  query = s.query,
  dims = 1:10,
  features = NULL
)



#s.anchors = FindTransferAnchors(reference = s.ref, query = s.query, dims = 1:15, features = common_features)
pedictions = TransferData(anchorset = s.anchors, refdata = s.ref$cell_state, dims = 1:10)

#Add transferred labels to s2 metadata

s.query = AddMetaData(s.query, metadata = pedictions)
#DimPlot(s2, group.by = "predicted.id", label = TRUE, label.size = 4)
p = DimPlot(s.query, group.by = "predicted.id", label = TRUE, repel = TRUE, label.size = 3, 
            pt.size = 0.7) + NoLegend() + ggtitle("Pax3+ Cells at E9.5")
  
###

# Extract embeddings and metadata:
umap_coords <- Embeddings(s.query, "umap") %>% 
  as.data.frame() %>% 
  mutate( ColorGroup = s.query@meta.data$seurat_clusters,
          LabelGroup = s.query@meta.data$predicted.id)
    
# Compute the median coordinates for labels
label_coords = umap_coords %>% 
  group_by(LabelGroup) %>%
  summarize(umap_1 = median(umap_1), umap_2 = median(umap_2))

# Plot cells colored by one metadata column and labeled by another
p = ggplot(umap_coords, aes(x = umap_1, y = umap_2, color = ColorGroup)) +
  geom_point(size = 0.5, alpha = 0.8) +
  geom_text_repel(data = label_coords, aes(label = LabelGroup, size = 0.5),
                  color = "black", max.overlaps = Inf) +
  theme_void() +
  #guides(color = guide_legend(title = "Color Group")) +
  NoLegend() 

p

################

# Compute median coordinates for labels
label_coords <- umap_coords %>%
  group_by(predicted.id) %>%
  summarize(umap_1 = median(umap_1), umap_2 = median(umap_2))

# Plot without labels first, then manually add labels:
ggplot(umap_coords, aes(umap_1, umap_2, color = predicted.id)) +
  geom_point(size = 0.7, alpha = 0.8) +
  NoLegend() +
  theme_void() +
  geom_text_repel(data = label_coords, aes(label = predicted.id),
                  size = 3, color = "black",
                  max.overlaps = Inf) +  
  
  theme(plot.margin = margin(5,5,5,5, "mm"))
 
ggplot(umap_coords, aes(umap_1, umap_2)) +
  geom_point(size = 0.7, alpha = 0.8) +
  NoLegend() +
  theme_void() +
  geom_text_repel(data = label_coords, aes(label = predicted.id),
                  size = 3, color = "black",
                  max.overlaps = Inf) +  
  
  theme(plot.margin = margin(5,5,5,5, "mm"))


saveRDS(s.query, sox_SC_Combined.rds)

#######################


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
DimPlot(s.query, reduction = "umap", label = TRUE, pt.size = 0.6,repel = TRUE) + NoLegend() + ggtitle("Pax3+ Cells at E9.5")



saveRDS(s.query,"pax95.rds")



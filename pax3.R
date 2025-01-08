#install.packages("dplyr")
#install.packages("Seurat")
#install.packages("patchwork")
#install.packages("writexl")
library(GPTCelltype)
library(openai)
library(dplyr)
library(Seurat)
library(patchwork)
library(writexl)
setwd("~/BINF/yushi scrnaseq/Pax3/pax3_RSEC_MolsPerCell_MEX")
pax3data = ReadMtx(mtx = "matrix.mtx.gz",
                   cells = "barcodes.tsv.gz",
                   features = "features.tsv.gz")

setwd("~/BINF/yushi scrnaseq/Pax3/trial")
pax3 = CreateSeuratObject(counts = pax3data, project = "pax3")
pax3[["percent.mt"]] <- PercentageFeatureSet(pax3, pattern = "^MT-")

vln = VlnPlot(pax3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pax3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pax3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pax3 <- subset(pax3, subset = nFeature_RNA > 1000 & nFeature_RNA < 9500 & percent.mt < 5)
pax3 <- NormalizeData(pax3)


pax3 <- FindVariableFeatures(pax3, selection.method = "vst", nfeatures = 5000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pax3), 10)

#scaling data
all.genes <- rownames(pax3)
pax3 <- ScaleData(pax3, features = all.genes)

#PCA
pax3 <- RunPCA(pax3, features = VariableFeatures(object = pax3))
heat = DimHeatmap(pax3, dims = 1:20, cells = 500, balanced = TRUE)
#variation until 14-15

elbow = ElbowPlot(pax3)
#elbows at 6, 14
#choose 14
pax3 = FindNeighbors(pax3, dims = 1:14)
#larger resolution for more no. of cells
pax3 <- FindClusters(pax3, resolution = 0.6)

#umap
pax3 = RunUMAP(pax3, dims = 1:14)
umapplot = DimPlot(pax3, reduction = "umap", label = TRUE)
#saveRDS(pax3, "pax3.rds")

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pax3.markers <- FindAllMarkers(pax3, only.pos = TRUE)
pax3.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

write_xlsx(pax3.markers,"pax3_markers_respoint6.xlsx")

new.cluster.ids <- c( "Paraxial mesenchyme cells",
                      "Neural progenitor cells", 
                      "Epithelial cells", 
                      "Neural crest cells",  
                      "Epicardium cells/Coronary Vascular Progenitor cells",  
                     "Forebrain neurons", 
                      "Craniofacial mesenchyme cells",  
                      "Schwann cell precursors or neural crest cells",  
                      "Endoderm or Hepatocyte precursors",  
                      "Lower urinary tract epithelial cells",
                      "Primitive Streak Cells",  
                      "Non-specialized/undifferentiated cells ", 
                      "Neural progenitor cells", 
                      "Cardiac progenitor cells ",
                      "Endothelial cells"  
)
           
SaveSeuratRds(pax3, file = "pax.rds")         
names(new.cluster.ids) <- levels(pax3)          
pax3 <- RenameIdents(pax3, new.cluster.ids)
DimPlot(pax3, reduction = "umap", label = TRUE, pt.size = 0.6,repel = TRUE) + NoLegend()     
?DimPlot           
#cluster 2 has enrichment of pax3, snai1, sox8, sox9, foxd3. notably not snai2. finding how much is the log2fc of snai2 in cluster 2

# Extract average expression per cluster
avg_exp <- AverageExpression(pax3, return.seurat = TRUE)

# Extract log-normalized expression values for Snai2
snai2_exp <- FetchData(pax3, vars = "Snai2")

# Add cluster identities
snai2_exp$cluster <- pax3$seurat_clusters

# Calculate average expression in cluster 2 and all other clusters
avg_exp_cluster2 <- mean(snai2_exp$Snai2[snai2_exp$cluster == 2])
avg_exp_other_clusters <- mean(snai2_exp$Snai2[snai2_exp$cluster != 2])

# Calculate log2 fold change
log2FC_snai2 <- log2(avg_exp_cluster2 + 1) - log2(avg_exp_other_clusters + 1)

log2FC_snai2 #-0.1539886

     
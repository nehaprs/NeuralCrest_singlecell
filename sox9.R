#sox9 cells at E9.5
library(dplyr)
library(Seurat)
library(patchwork)
library(writexl)
library(readxl)
library(clustree)
setwd("~/BINF/yushi scrnaseq/E9.5/Sox9/Sox9_RSEC_MolsPerCell_MEX")

sox9data = ReadMtx(mtx = "matrix.mtx.gz",
                   cells = "barcodes.tsv.gz",
                   features = "features.tsv.gz")
setwd("~/BINF/yushi scrnaseq/E9.5/Sox9/seurat output")
sox9 = CreateSeuratObject(counts = sox9data, project = "sox9")
sox9[["percent.mt"]] <- PercentageFeatureSet(sox9, pattern = "^MT-")

vln = VlnPlot(sox9, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(sox9, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sox9, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#more conservative version
setwd("~/BINF/yushi scrnaseq/E9.5/Sox9/seurat output/conservative")
sox9 <- subset(sox9, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & percent.mt < 5)
#was 1000 and 8000 in previous version
sox9 <- NormalizeData(sox9)

sox9  <- FindVariableFeatures(sox9 , selection.method = "vst", nfeatures = 4000)
#nFeat = 5k earlier, whose minimum variance was 1.066
#for nfeeatures = 4k, it is 1.15
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sox9), 10)

#scaling data
all.genes <- rownames(sox9)
sox9 <- ScaleData(sox9, features = all.genes)

#PCA
sox9 <- RunPCA(sox9, features = VariableFeatures(object = sox9))
heat = DimHeatmap(sox9, dims = 1:20, cells = 500, balanced = TRUE)
#11,12, or 17?

elbow = ElbowPlot(sox9)
#elbow at 12, 14?
#chose 14

sox9 = FindNeighbors(sox9, dims = 1:14)
#larger resolution for more no. of cells
sox9 <- FindClusters(sox9, resolution = 0.5)
#res = 0.5, 18 clusters

#umap
sox9 = RunUMAP(sox9, dims = 1:14)
umapplot = DimPlot(sox9, reduction = "umap", label = TRUE)
saveRDS(sox9, "sox9.rds")

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
sox9.markers <- FindAllMarkers(sox9, only.pos = TRUE)
sox9.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

write_xlsx(sox9.markers,"sox9_markers_respoint5.xlsx")

#finding the cut-off of variable features
sox9variables = VariableFeatures(sox9)

variable_feature = list()

standardizedVariance = list()

for(i in 1:29437){
  if (sox9@assays$RNA@meta.data$var.features[i] %in% sox9variables){
    variable_feature[[length(variable_feature) + 1]] <- sox9@assays$RNA@meta.data$var.features[i]
    standardizedVariance[[length(standardizedVariance) + 1]] <- sox9@assays$RNA@meta.data$vf_vst_counts_variance.standardized[i]
   
  }
}

variable_feature_std_var = cbind(variable_feature, standardizedVariance)
variable_feature_std_var[,2] = as.numeric(variable_feature_std_var[,2])
#order them lowest to highest
sorted_variable_features <- variable_feature_std_var[order(variable_feature_std_var$standardizedVariance, decreasing = TRUE), ]

sorted_variable_feature = as.data.frame(variable_feature_std_var[as.vector(order(variable_feature_std_var$standardizedVariance), decreasing = TRUE), ])
write.table(variable_feature_std_var,"variable_feature_std_var")
#highest std variance = 33.2565927497255 lowest std variance = 1.06617979222262


###############################
#Do clustree for sox data
setwd("~/BINF/yushi scrnaseq/Sox9/clustree/highres")
sox <- readRDS("~/BINF/yushi scrnaseq/Sox9/seurat output/round1/sox9.rds")

resolution.range <- seq(from = 0.8, to = 1.5, by = 0.1)
# Loop over each resolution
for (res in resolution.range) {
  # Perform clustering with the current resolution
  sox <- FindClusters(sox, resolution = res)
  
  # Find all markers for the clusters at this resolution
  sox.markers <- FindAllMarkers(sox, only.pos = TRUE)
  
  # Define the file name for saving the markers
  file_name <- paste0("markers_resolution_", res, ".xlsx")
  
  # Save the markers as an Excel file
  write_xlsx(sox.markers, file_name)
  
  # Print a message to confirm completion for each resolution
  print(paste("Markers for resolution", res, "saved to", file_name))
}

#list all xlsx files in wd

xlsx_file = list.files(pattern = "\\.xlsx$")

for (file in xlsx_file){
  df = read_xlsx(file)
  dff = df[df$avg_log2FC > 1,]
  file_new = paste0("filt",file)
  write_xlsx(dff, file_new)
}

soxclust = clustree(sox)
#draw umap for res 0.1
?DimPlot
DimPlot(sox, reduction = "umap", label = TRUE,
        group.by = "RNA_snn_res.0.1")
saveRDS(sox, "sox9highres.rds")
#find expression of marker genes
migratory.markers = c("Foxd3", "Ets1", "Sox8", "Sox9", "Pax7", "Tfap2a",
                      "Pax3", "Sox10", "Lmo4", "Rxrg", "Ltk", "Erbb3")

sox.4 = SetIdent(sox, value = "RNA_snn_res.0.4")
sox.4 = AddModuleScore(sox.4, features = migratory.markers)

#Addmodulescore@features is case-sensitive

FeaturePlot(sox.4, reduction = "umap", pt.size = 0.5,
            features = migratory.markers)

#try clustree with a diff set of resolution

resolution.range2 <- seq(from = 0, to = 1, by = 0.05)
sox <- FindClusters(sox, resolution = resolution.range)
soxclust = clustree(sox)

#save workspace
save.image(file='sox9.RData')

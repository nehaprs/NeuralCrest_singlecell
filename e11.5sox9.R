library(dplyr)
library(Seurat)
library(patchwork)
library(writexl)
library(readxl)
library(clustree)
library(GPTCelltype)
library(openai)

setwd("~/BINF/yushi scrnaseq/E11.5/sox9/1-sox9-11-5_clean_RSEC_MolsPerCell_MEX")
sox9data = ReadMtx(mtx = "matrix.mtx.gz",
                   cells = "barcodes.tsv.gz",
                   features = "features.tsv.gz")
setwd("~/BINF/yushi scrnaseq/E11.5/sox9/seurat output")
sox9 = CreateSeuratObject(counts = sox9data, project = "sox9")
sox9[["percent.mt"]] <- PercentageFeatureSet(sox9, pattern = "^MT-")

vln = VlnPlot(sox9, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(sox9, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sox9, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

sox9 <- subset(sox9, subset = nFeature_RNA > 2000 & nFeature_RNA < 7000 & percent.mt < 5)
sox9 <- NormalizeData(sox9)

sox9  <- FindVariableFeatures(sox9 , selection.method = "vst", nfeatures = 4000)

top10 <- head(VariableFeatures(sox9), 10)

#scaling data
all.genes <- rownames(sox9)
sox9 <- ScaleData(sox9, features = all.genes)

#PCA
sox9 <- RunPCA(sox9, features = VariableFeatures(object = sox9))
heat = DimHeatmap(sox9, dims = 1:20, cells = 500, balanced = TRUE)
#8 or 13?

elbow = ElbowPlot(sox9) #12

sox9 = FindNeighbors(sox9, dims = 1:12)

resolution.range <- seq(from = 0, to = 1, by = 0.1)

# Loop over each resolution
for (res in resolution.range) {
  # Perform clustering with the current resolution
  sox9<- FindClusters(sox9, resolution = res)
  
  # Find all markers for the clusters at this resolution
  sox.markers <- FindAllMarkers(sox9, only.pos = TRUE)
  
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

soxclust = clustree(sox9)
sox9 = RunUMAP(sox9, dims = 1:12)

DimPlot(sox9, reduction = "umap", label = TRUE,
        group.by = "RNA_snn_res.0.1")
migratory.markers = c("Foxd3", "Ets1", "Sox8", "Sox9", "Pax7", "Tfap2a",
                      "Pax3", "Sox10", "Lmo4", "Rxrg", "Ltk", "Erbb3")


sox.1 = SetIdent(sox9, value = "RNA_snn_res.0.1")
sox.1 = AddModuleScore(sox.1, features = migratory.markers)
FeaturePlot(sox.1, reduction = "umap", pt.size = 0.5,
            features = migratory.markers)




chondrocyte.markers = c("Smad3", "Sox9", "Pax3", "Sox5", "Sox6", 
                        "Col2a1")

sox.1 = SetIdent(sox9, value = "RNA_snn_res.0.1")
sox.1 = AddModuleScore(sox.1, features = chondrocyte.markers)
FeaturePlot(sox.1, reduction = "umap", pt.size = 0.5,
            features = chondrocyte.markers)


melanocyte.markers = c( "Sox9", "Pax3", "Sox10", "Pax7",
                        "Mitf", "Mef2c", "Dct", "Tyr")

sox.1 = SetIdent(sox9, value = "RNA_snn_res.0.1")
sox.1 = AddModuleScore(sox.1, features = melanocyte.markers)
FeaturePlot(sox.1, reduction = "umap", pt.size = 0.5,
            features = melanocyte.markers)

sym.neuronal.markers = c( "Sox9", "Pax3", "Sox10", "Smad4", "Insm1",
                      "Ascl1", "Phox2b", "Hand2", "Gata3", "Th", "Dbh")

sox.1 = SetIdent(sox9, value = "RNA_snn_res.0.1")
sox.1 = AddModuleScore(sox.1, features = sym.neuronal.markers)
FeaturePlot(sox.1, reduction = "umap", pt.size = 0.5,
            features = sym.neuronal.markers)
#res.6

DimPlot(sox9, reduction = "umap", label = TRUE,
        group.by = "RNA_snn_res.0.6")
#nc markers likely
sox.6 = SetIdent(sox9, value = "RNA_snn_res.0.6")

#annotation using gpt4
Sys.setenv(OPENAI_API_KEY = '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')

sox_6_markers <- read_excel("filtmarkers_resolution_0.6.xlsx")
sox_6_markers = as.data.frame(sox_6_markers)
res6 <- gptcelltype(sox_6_markers, tissuename = 'mouse embryo at E11.5', 
                    model = 'gpt-4')

new.cluster.ids =c( "Presomitic mesoderm cells", 
                    "Intermediate mesoderm cells/Kidney progenitor cells" ,
                    "Mitotic cells",
                    "Lateral plate mesoderm cells",
                    "Skeletal muscle progenitor cells",
                    "Neural stem/progenitor cells",
                    "Meningeal cells",
                    "Metanephric mesenchyme cells",
                    "Neural progenitor cells (Cortical)",
                    "Chondrocytes",
                    "Epicardium cells",
                    "Surface ectoderm cells/Epithelial cells",
                    "Urothelial cells",
                    "Neural crest cells",
                    "GABAergic neurons",
                    " Erythroid progenitor cells",
                    "Ependymal cells(Choroid plexus)")


sox9 = SetIdent(sox9, value = "RNA_snn_res.0.6")
names(new.cluster.ids) <- levels(sox9)
sox9 <- RenameIdents(sox9, new.cluster.ids)
DimPlot(sox9, reduction = "umap", label = TRUE, pt.size = 0.5,repel = TRUE
) + NoLegend()
saveRDS(sox9, "sox9.rds")

#BiocManager::install("SingleCellExperiment")
#BiocManager::install("celldex")
library(SingleR)
library(SingleCellExperiment)
library(celldex)
sces11 <- as.SingleCellExperiment(sox9)
ref <- MouseRNAseqData()
pred <- SingleR(test = sces11, ref = ref, labels = ref$label.main)
sox9$SingleR.labels <- pred$labels
'
names(sox9$SingleR.labels)
sox9 = SetIdent(sox9, value = "RNA_snn_res.0.6")
names(sox9$SingleR.labels) <- levels(sox9)
sox9 <- RenameIdents(sox9, SingleR.labels)
DimPlot(sox9, reduction = "umap", 
        group.by = "SingleR.labels", label = TRUE, repel = TRUE) + NoLegend()
'

cluster_annotation_table <- table(
  Cluster = Idents(sox9),
  SingleR_Label = sox9$SingleR.labels
)
library(dplyr)

# Convert the table to a data frame for easier manipulation
cluster_annotation_df <- as.data.frame(cluster_annotation_table)

# For each cluster, find the SingleR label with the maximum count
cluster_celltype <- cluster_annotation_df %>%
  group_by(Cluster) %>%
  summarize(CellType = SingleR_Label[which.max(Freq)])

cluster_ids <- as.character(cluster_celltype$Cluster)
cell_types <- as.character(cluster_celltype$CellType)
names(cell_types) <- cluster_ids

cluster_ids <- as.character(cluster_celltype$Cluster)
cell_types <- as.character(cluster_celltype$CellType)
names(cell_types) <- cluster_ids

sox9$ClusterCellType <- plyr::mapvalues(
  x = as.character(Idents(sox9)),
  from = cluster_ids,
  to = cell_types
)


Idents(sox9) <- "ClusterCellType"

DimPlot(
  object = sox9,
  reduction = "umap",
  label = TRUE,
  label.size = 4,
  repel = TRUE
) +
  ggtitle("UMAP Plot of Clusters Labeled by singleR Annotated Cell Types")

##now use mouse gastrulation atlas data for reference

#BiocManager::install("MouseGastrulationData")
#BiocManager::install("scuttle")

library(MouseGastrulationData)
library(scuttle)

# Load the reference dataset
reference <- EmbryoAtlasData()

# Normalize reference data (if not already normalized)
reference <- logNormCounts(reference)

# Convert Seurat object to SingleCellExperiment
query_sce <- as.SingleCellExperiment(sox9)
rm(sces11)
# Normalize query data
query_sce <- logNormCounts(query_sce)

predictions <- SingleR(
  test = query_sce,
  ref = reference,
  labels = reference$celltype.mapped  # Adjust 'celltype.mapped' based on available metadata
)

pred <- SingleR(test = sces11, ref = ref, labels = ref$label.main)
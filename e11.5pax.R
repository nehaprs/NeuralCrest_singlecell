library(dplyr)
library(Seurat)
library(patchwork)
library(writexl)
library(readxl)
library(clustree)
library(GPTCelltype)
library(openai)

setwd("~/BINF/yushi scrnaseq/E11.5/pax3/3-pax3-10-5_clean_RSEC_MolsPerCell_MEX")
pax3data = ReadMtx(mtx = "matrix.mtx.gz",
                   cells = "barcodes.tsv.gz",
                   features = "features.tsv.gz")
setwd("~/BINF/yushi scrnaseq/E11.5/pax3/seurat output")

pax3 = CreateSeuratObject(counts = pax3data, project = "pax3")
pax3[["percent.mt"]] <- PercentageFeatureSet(pax3, pattern = "^MT-")

vln = VlnPlot(pax3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pax3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pax3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pax3 <- subset(pax3, subset = nFeature_RNA > 2500 & nFeature_RNA < 7500 & percent.mt < 5)
pax3 <- NormalizeData(pax3)

pax3  <- FindVariableFeatures(pax3 , selection.method = "vst", nfeatures = 5000)

top10 <- head(VariableFeatures(pax3), 10)

#scaling data
all.genes <- rownames(pax3)
pax3 <- ScaleData(pax3, features = all.genes)

#PCA
pax3 <- RunPCA(pax3, features = VariableFeatures(object = pax3))
heat = DimHeatmap(pax3, dims = 1:20, cells = 500, balanced = TRUE)
#11-13 ?

elbow = ElbowPlot(pax3)
#12

pax3 = FindNeighbors(pax3, dims = 1:12)

resolution.range <- seq(from = 0, to = 1, by = 0.1)

# Loop over each resolution
for (res in resolution.range) {
  # Perform clustering with the current resolution
  pax3<- FindClusters(pax3, resolution = res)
  
  # Find all markers for the clusters at this resolution
  pax.markers <- FindAllMarkers(pax3, only.pos = TRUE)
  
  # Define the file name for saving the markers
  file_name <- paste0("paxmarkers_resolution_", res, ".xlsx")
  
  # Save the markers as an Excel file
  write_xlsx(pax.markers, file_name)
  
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

paxclust = clustree(pax3)
pax3 = RunUMAP(pax3, dims = 1:12)

DimPlot(pax3, reduction = "umap", label = TRUE,
        group.by = "RNA_snn_res.0.1")

#find markers + sox9 +pax3 manually on res.1, as a starting point

chondrocyte.markers = c("Smad3", "Sox9", "Pax3", "Sox5", "Sox6", 
                        "Col2a1")

pax.1 = SetIdent(pax3, value = "RNA_snn_res.0.1")
pax.1 = AddModuleScore(pax.1, features = chondrocyte.markers)
FeaturePlot(pax.1, reduction = "umap", pt.size = 0.5,
            features = chondrocyte.markers)


melanocyte.markers = c( "Sox9", "Pax3", "Sox10", "Pax7",
                        "Mitf", "Mef2c", "Dct", "Tyr")

pax.1 = SetIdent(pax3, value = "RNA_snn_res.0.1")
pax.1 = AddModuleScore(pax.1, features = melanocyte.markers)
FeaturePlot(pax.1, reduction = "umap", pt.size = 0.5,
            features = melanocyte.markers)

neuronal.markers = c( "Sox9", "Pax3", "Sox10", "Smad4", "Insm1",
                      "Ascl1", "Phox2b", "Hand2", "Gata3", "Th", "Dbh")

pax.1 = SetIdent(pax3, value = "RNA_snn_res.0.1")
pax.1 = AddModuleScore(pax.1, features = neuronal.markers)
FeaturePlot(pax.1, reduction = "umap", pt.size = 0.5,
            features = neuronal.markers)

#clusters 0,1,2,4 tend to mean something. 
#choose 0.4, before clustree goes crazy

DimPlot(pax3, reduction = "umap", label = TRUE,
        group.by = "RNA_snn_res.0.4")

pax.3 = SetIdent(pax3, value = "RNA_snn_res.0.4")

#annotation using gpt4
Sys.setenv(OPENAI_API_KEY = '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')

pax_4_markers <- as.data.frame(read_excel("filtpaxmarkers_resolution_0.4.xlsx"))
res4.3 <- gptcelltype(pax_4_markers, tissuename = 'mouse embryo at Embryonic Day 11.5' , 
                    model = 'gpt-4')

new.cluster.ids = c("Mesodermal progenitor cells",
                    "Mesenchymal Stem Cells",
                    "Neural Progenitor Cells",
                    "Neural Progenitor Cells",
                    "Paraxial Mesoderm Cells",
                    "Immature Neurons",
                    "Neuroendocrine Cells",
                    "Neural Crest Cells" ,
                    "Adrenergic Neurons",
                     "Epidermal Cells",
                    "Neural Plate Border Specifiers",
                    "Primitive Erythroid Cells"
)

pax3 = SetIdent(pax3, value = "RNA_snn_res.0.4")
names(new.cluster.ids) <- levels(pax3)
pax3 <- RenameIdents(pax3, new.cluster.ids)
DimPlot(pax3, reduction = "umap", label = TRUE, pt.size = 0.7,repel = TRUE
) + NoLegend()
saveRDS(pax3, "pax3.rds")

#load RDS and find nc markers
pax3 <- readRDS("~/BINF/yushi scrnaseq/E11.5/pax3/seurat output/pax3.rds")

migratorync.markers = c("Foxd3", "Ets1", "Sox8", "Sox9", "Pax7", "Tfap2a",
                      "Pax3", "Sox10", "Lmo4", "Rxrg", "Ltk", "Erbb3")

pax.1 = SetIdent(pax3, value = "RNA_snn_res.0.1")
pax.1 = AddModuleScore(pax.1, features =migratorync.markers)
FeaturePlot(pax.1, reduction = "umap", pt.size = 0.5,
            features = migratorync.markers)

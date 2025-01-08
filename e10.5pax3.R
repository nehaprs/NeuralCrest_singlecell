library(dplyr)
library(Seurat)
library(patchwork)
library(writexl)
library(readxl)
library(clustree)
library(GPTCelltype)
library(openai)

setwd("~/BINF/yushi scrnaseq/E10.5/Pax3/3-pax3-10-5_clean_RSEC_MolsPerCell_MEX")
pax3data = ReadMtx(mtx = "matrix.mtx.gz",
                   cells = "barcodes.tsv.gz",
                   features = "features.tsv.gz")
setwd("~/BINF/yushi scrnaseq/E10.5/Pax3/seurat output")
pax3 = CreateSeuratObject(counts = pax3data, project = "pax3")
pax3[["percent.mt"]] <- PercentageFeatureSet(pax3, pattern = "^MT-")

vln = VlnPlot(pax3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pax3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pax3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pax3 <- subset(pax3, subset = nFeature_RNA > 4000 & nFeature_RNA < 8500 & percent.mt < 5)
pax3 <- NormalizeData(pax3)

pax3  <- FindVariableFeatures(pax3 , selection.method = "vst", nfeatures = 4000)

top10 <- head(VariableFeatures(pax3), 10)

#scaling data
all.genes <- rownames(pax3)
pax3 <- ScaleData(pax3, features = all.genes)

#PCA
pax3 <- RunPCA(pax3, features = VariableFeatures(object = pax3))
heat = DimHeatmap(pax3, dims = 1:20, cells = 500, balanced = TRUE)
#12 ?

elbow = ElbowPlot(pax3)
#12 urappichu


pax3 = FindNeighbors(pax3, dims = 1:12)
###code changes to only one resolution
#resolution.range <- seq(from = 0, to = 1, by = 0.1)
resolution.range <- 0.2
# Loop over each resolution
for (res in resolution.range) {
  # Perform clustering with the current resolution
  pax3<- FindClusters(pax3, resolution = res)
  
  # Find all markers for the clusters at this resolution
  #pax.markers <- FindAllMarkers(pax3, only.pos = TRUE)
  
  # Define the file name for saving the markers
  #file_name <- paste0("paxmarkers_resolution_", res, ".xlsx")
  
  # Save the markers as an Excel file
  #write_xlsx(pax.markers, file_name)
  
  # Print a message to confirm completion for each resolution
  #print(paste("Markers for resolution", res, "saved to", file_name))
}

#list all xlsx files in wd

xlsx_file = list.files(pattern = "\\.xlsx$")

for (file in xlsx_file){
  df = read_xlsx(file)
  dff = df[df$avg_log2FC > 1,]
  file_new = paste0("filt",file)
  write_xlsx(dff, file_new)
}

# Loop over each resolution
for (res in resolution.range) {
  # Perform clustering with the current resolution
  pax3<- FindClusters(pax3, resolution = res)}
paxclust = clustree(pax3)
pax3 = RunUMAP(pax3, dims = 1:12)



#from res 0.1, find nc like clusters. Then find the subclustrs at the optimal res

DimPlot(pax3, reduction = "umap", label = TRUE,
        group.by = "RNA_snn_res.0.1")
migratory.markers = c("Foxd3", "Ets1", "Sox8", "Sox9", "Pax7", "Tfap2a",
                      "Pax3", "Sox10", "Lmo4", "Rxrg", "Ltk", "Erbb3")

pax.1 = SetIdent(pax3, value = "RNA_snn_res.0.1")
pax.1 = AddModuleScore(pax.1, features = migratory.markers)
FeaturePlot(pax.1, reduction = "umap", pt.size = 0.5,
            features = migratory.markers)
#focus on parent cluster 4. may be to some extend also 3
#res .2 chosen

DimPlot(pax3, reduction = "umap", label = TRUE,
        group.by = "RNA_snn_res.0.2")
#nc markers likely 3. Also considering 6 and 7
pax.3 = SetIdent(pax3, value = "RNA_snn_res.0.2")

#annotation using gpt4
Sys.setenv(OPENAI_API_KEY = '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')

pax_2_markers <- as.data.frame(read_excel("filtpaxmarkers_resolution_0.2.xlsx"))
res2 <- gptcelltype(pax_2_markers, tissuename = 'mouse embryo', 
                    model = 'gpt-4')
new.cluster.ids = c("Mesenchymal Stem Cells",
                    "Neural Progenitor Cells",
                    "Craniofacial Mesenchyme Cells",
                    "Neural Crest Cells",
                    "Neurons",
                    "Embryonic Epithelial Cells",
                    "Epithelial Cells",
                    "Epidermal Cells",
                    "Paraxial Mesoderm Cells",
                    "Erythrocytes",
                    "T Lymphocytes" ,
                    "Hematopoietic Stem Cells" 
  
)

pax3 = SetIdent(pax3, value = "RNA_snn_res.0.2")
names(new.cluster.ids) <- levels(pax3)
pax3 <- RenameIdents(pax3, new.cluster.ids)
DimPlot(pax3, reduction = "umap", label = TRUE, pt.size = 0.5,repel = TRUE
) + NoLegend()
saveRDS(pax3, "pax3.rds")
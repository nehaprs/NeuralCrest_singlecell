#integrate three pairs of data first
#find double positives
#then do monocle pseudotime
library(dplyr)
library(Seurat)
library(patchwork)
library(writexl)
library(readxl)
library(clustree)

setwd("~/BINF/yushi scrnaseq/all six/e105")
#load the surat objects
sox105 <- readRDS("~/BINF/yushi scrnaseq/New folder/rds objects/all clusters/sox105.rds")
pax105 <- readRDS("~/BINF/yushi scrnaseq/New folder/rds objects/all clusters/pax105.rds")

s.list = list(sox105, pax105)
# normalize and identify variable features for each dataset independently
s.list = lapply(X= s.list, 
                FUN = function(x){
                  x = NormalizeData(x)
                  x = FindVariableFeatures(x, selection.method = "vst")
                })


features <- SelectIntegrationFeatures(object.list = s.list)

####
s.list <- lapply(X = s.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

s.anchors <- FindIntegrationAnchors(object.list = s.list, anchor.features = features, reduction = "rpca")
# this command creates an 'integrated' data assay
s.combined <- IntegrateData(anchorset = s.anchors)


#run a single integrated analysis
DefaultAssay(s.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
s.combined <- ScaleData(s.combined, verbose = FALSE)
s.combined <- RunPCA(s.combined, npcs = 30, verbose = FALSE)
s.combined <- RunUMAP(s.combined, reduction = "pca", dims = 1:30)
s.combined <- FindNeighbors(s.combined, reduction = "pca", dims = 1:30)
s.combined <- FindClusters(s.combined, resolution = 0.5)

p1 <- DimPlot(s.combined, reduction = "umap", group.by = "orig.ident")+ggtitle("Combined Dataset at E10.5 by Origin")
p2 <- DimPlot(s.combined, reduction = "umap", group.by = "predicted.id", label = TRUE,
              repel = TRUE, pt.size = 1) + NoLegend() +ggtitle("Combined Dataset at E10.5")

##find cells with minimus sox9 and pax3
##subcluster
head(rownames(s.combined))
#genes are in genename format

#find minimum sox9 and pax3 expressions
sox9Cells = subset(s.combined, orig.ident == "sox9")
sox9_expr = FetchData(sox9Cells, vars = "Sox9", assay = "RNA")
any(sox9_expr < 0)
# scatterplot  SOX9 expression
sox9_expr$cell <- rownames(sox9_expr)
sox9_expr$index <- seq_len(nrow(sox9_expr))

# Scatter plot

ggplot(sox9_expr, aes(x = index, y = Sox9)) +
  geom_point(alpha = 0.6, color = "darkgreen") +
  labs(title = "Sox9 Expression in Cells with orig.ident == 'sox9' at E10.5",
       x = "Cell Index",
       y = "Sox9 Expression (log-normalized)") +
  theme_minimal()

#find the smallest non-zero value
non_zero_values <- sox9_expr$Sox9[sox9_expr$Sox9 > 0]
min_non_zero <- min(non_zero_values)
min(sox9_expr$Sox9)
#min non-zero value at E9.5 = 0.211
#E10.5 = 0.2557998

unique(s.combined$orig.ident)

pax3Cells = subset(s.combined, orig.ident == "pax3")
pax_expr = FetchData(pax3Cells, vars ="Pax3", assay = "RNA")
# scatterplot expression
pax_expr$cell = rownames(pax_expr)
pax_expr$index = seq_len(nrow(pax_expr))

# Scatter plot

ggplot(pax_expr, aes(x = index, y = Pax3)) +
  geom_point(alpha = 0.6, color = "darkblue") +
  labs(title = "Pax3 Expression in Cells with orig.ident == 'pax3' at E10.5",
       x = "Cell Index",
       y = "Pax3 Expression (log-normalized)") +
  theme_minimal()

#find the smallest non-zero value

non_zero_values = pax_expr$Pax3[pax_expr$Pax3 > 0]
min_non_zero <- min(non_zero_values)
min(pax_expr$Pax3)
#2.979014e-05
sort(non_zero_values)[1:10]
#lowest 10 values: 2.979014e-05 1.359133e-04 1.137102e-03 1.458762e-03 1.460009e-03 2.597825e-03 2.891456e-03 3.048795e-03 3.064116e-03
# 3.095681e-03

#choose 1.137102e-03?


soxx = FetchData(sox95, vars = "Sox9")
nzv = soxxx$Sox9[soxx$Sox9 == 0]
#2729 cells with zero sox9

rm(soxx, nzv)

### min sox9 = 0.211
### min pax3 = try 2.979014e-05 or  1.137102e-03


#find pax3 expr in sox9 cells
paxinsox = FetchData(sox9Cells, vars = "Pax3")
non_zero_values = pax_expr$Pax3[pax_expr$Pax3 > 0]

filtpaxinsox = paxinsox$Pax3[paxinsox$Pax3 > 2.979014e-05]
#3824 cells have non-zeo values. same number above both thresholds


#find cells with value above the threshold for both pax and sox

# Subset the double positives
dbpos = subset(s.combined, subset = Pax3 > 1.137102e-03 & Sox9 > 0.211)
#there are 4895 double positive cells
# Count how many cells in the subset have each orig.ident value

table(dbpos$orig.ident)

'''
pax3 sox9 
1663 3232
'''

p1 <- DimPlot(dbpos, reduction = "umap", group.by = "orig.ident")+ggtitle("Double Positives at E9.5 by Origin")
p2 <- DimPlot(dbpos, reduction = "umap", group.by = "predicted.id", label = TRUE,
              repel = TRUE, pt.size = 1) + NoLegend() +ggtitle("Double Positives at E9.5")

#preprocess dbpos


DefaultAssay(dbpos) = "RNA"
dbpos <- NormalizeData(dbpos)

dbpos  <- FindVariableFeatures(dbpos , selection.method = "vst")

top10 <- head(VariableFeatures(dbpos), 10)

#scaling data
all.genes <- rownames(dbpos)
dbpos <- ScaleData(dbpos, features = all.genes)

#PCA
dbpos <- RunPCA(dbpos, features = VariableFeatures(object = dbpos))
heat = DimHeatmap(dbpos, dims = 1:20, cells = 500, balanced = TRUE)

elbow = ElbowPlot(dbpos) #15

dbpos = FindNeighbors(dbpos, dims = 1:15)
dbpos <- RunUMAP(dbpos, reduction = "pca", dims = 1:30)

dbpos <- FindClusters(dbpos, resolution = 0.9)

p1 <- DimPlot(dbpos, reduction = "umap", group.by = "orig.ident")+ggtitle("double positives at E9.5 by Origin")
p2 <- DimPlot(dbpos, reduction = "umap", group.by = "predicted.id", label = TRUE,
              repel = TRUE, pt.size = 1) + NoLegend() +ggtitle("double positives at E9.5")

DimPlot(dbpos, reduction = "umap")
saveRDS(dbpos,"dbpos95.rds")

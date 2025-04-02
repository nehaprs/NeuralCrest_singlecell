#annotate integrated object
#one query with 3 edays, 3 references.
#archived. method doesn't work
library(Seurat)
library(writexl)
library(dplyr)
library(ggplot2)
setwd("~/BINF/yushi scrnaseq/time series/harmony_slingshot")
query <- readRDS("~/BINF/yushi scrnaseq/time series/harmony_slingshot/soxCombined.rds")

any(duplicated(rownames(query)))
#FALSE
'''
duplicated_genes_query <- rownames(query)[duplicated(rownames(query))]
query <- query[!is.na(rownames(query)), ]
any(duplicated(rownames(query)))
'''


ref95 <- readRDS("~/BINF/yushi scrnaseq/E9.5/tomeE9.5.rds")
ref105 <- readRDS("~/BINF/yushi scrnaseq/E10.5/tome_E10.5.rds")
ref115 <- readRDS("~/BINF/yushi scrnaseq/E11.5/tomeRef_E11.5.rds")

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
ref95 = process(ref95)
ref105 = process(ref105)
ref115 = process(ref115)

table(query$eday)

#convert query rownames from geneID to ensembl

# First, split the query object by eday:
query_list <- SplitObject(query, split.by = "eday")
# Create an empty list to store annotated objects
annotated_query <- list()


#process

# Define a named list of references corresponding to each eday
ref_list = list("E9.5" = ref95,
                "E10.5" = ref105,
                "E11.5" = ref115)





# Loop over each subset and annotate separately

for(eday in names(query_list)){
  
  curr_query = query_list[[eday]]
  curr_ref = ref_list[[eday]]
  
  s.anchors <- FindTransferAnchors(
    reference = curr_ref,
    query = curr_query,
    dims = 1:10,
    features = NULL
  )
  
  
  predictions = TransferData(anchorset = s.anchors, refdata = curr_ref$cell_type, dims = 1:10)
  
  curr_query <- AddMetaData(curr_query, metadata = predictions)
  
  # Save into the list
  annotated_query[[eday]] <- curr_query
  
  cat("Completed annotation for", eday, "\n")
}



# Merge annotated query subsets back into a single Seurat object
annotated_query_merged <- merge(annotated_query[[1]], 
                                y = annotated_query[2:length(annotated_query)], 
                                merge.data = TRUE)

# Verify results:
table(annotated_query_merged$eday, annotated_query_merged$predicted.id)




ensemble_toId(annotated_query_merged) #convert to geneID
annotated_query_merged <-process(annotated_query_merged)
annotated_query_merged <- RunUMAP(annotated_query_merged, dims = 1:30, reduction = "harmony")

DimPlot(annotated_query_merged, group.by = "predicted.id",  label = TRUE, repel = TRUE, label.size = 4, pt.size = 1.5) + NoLegend() + ggtitle("merged Sox9+ Cells")








saveRDS(annotated_query_merged,"soxCombinedAnntd.rds")


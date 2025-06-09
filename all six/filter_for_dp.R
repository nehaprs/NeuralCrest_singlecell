#==========================================================================
#filters for cells with double positives, ie, Sox9 and Pax3 expression > 0
#===========================================================================

library(Seurat)
'
add the following information in sproc:
for every cell in sproc, the value of the expression levels of genes Sox9 and Pax3 of the same cell in sraw.
'

filter_for_dp <- function(sproc, sraw) {
  
  # Extract raw expression for Sox9 and Pax3
  raw_expr <- FetchData(sraw, vars = c("Sox9", "Pax3"), slot = "counts")
  
  # Keep only cells that are in both objects
  common_cells <- intersect(Cells(sproc), rownames(raw_expr))
  raw_expr <- raw_expr[common_cells, ]
  sproc <- subset(sproc, cells = common_cells)
  
  # Rename to avoid conflicts and add to metadata
  colnames(raw_expr) <- paste0("raw_", colnames(raw_expr))
  sproc <- AddMetaData(sproc, metadata = raw_expr)
  
  # Subset: keep only cells where both genes have non-zero expression
  sproc <- subset(sproc, subset = raw_Sox9 > 0 & raw_Pax3 > 0)
  
  
  
  # Return the filtered Seurat object
  return(sproc)
}

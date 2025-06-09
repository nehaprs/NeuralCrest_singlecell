library(dplyr)
library(Seurat)
'
a function with the inputs being a seurat object , and an integrer n, 
as the following: output is the same seurat object with its clusters annotated as the following. 
the cells already have a corresponding predicted id in the seurat object. 
for every cluster, if  > n % cells have the same predicted.id associated with it, rename it as that predicted.id.
if none of the predicted.id occurs for more than n% of the cells, 
the cluster is left as it is, (which is likely the number assigned by suerat  clustering.)

'
annotate_clusters_by_pred = function(seurat_obj, n)
  
  # seurat_obj:      a Seurat object containing a "predicted.id" column in meta.data
  # n:               numeric threshold (percent). If > n% of cells in a cluster share the same
  #                  predicted.id, that cluster is renamed to that predicted.id.
  #
  # Returns:         the same Seurat object, with Idents() updated. 
  #                  The original cluster labels are saved in meta.data$orig_cluster.
  
  
  #initial checks
{ 
    if(!"predicted.id" %in% colnames(seurat_obj@meta.data)){
    stop("Error: No column 'predicted.id'")
    }
  
  #copy original ids as a new metadata column, so that it's saved for later
 cluster_orig = Idents(seurat_obj)
 seurat_obj$orig_cluster <- as.character(cluster_orig)
 
 meta = seurat_obj@meta.data
 

 # Force predicted IDs and original clusters to character type
  meta$predicted.id = as.character(meta$predicted.id)
  meta$orig_cluster = as.character(meta$orig_cluster)
  
  # Create a named vector of current cluster assignments
  # (this will be updated cluster-by-cluster)
  current_idents = as.character(cluster_orig)
  names(current_idents) = names(cluster_orig) #cell names
 

  ##For each original cluster, check the fraction of cells per predicted.id
  
  unique_clusters = unique(meta$orig_cluster)
  for (cl in unique_clusters) {
    # Get all cells belonging to this original cluster
    cells_in_cl = rownames(meta)[meta$orig_cluster == cl]
    # Extract their predicted.id values
    preds = meta[cells_in_cl, "predicted.id"]
    
    # Build a frequency table of predicted.id within this cluster
    tbl <- table(preds)
    # Convert to percentages
    pct <- tbl / length(cells_in_cl) * 100
    
    # Find the predicted.id with the highest percentage
    max_pred <- names(pct)[which.max(pct)]
    max_pct  <- max(pct)
    
    # If the top predicted.id exceeds threshold n, reassign cluster to that predicted.id
    if (max_pct > n) {
      # Assign all cells in this cluster to `max_pred`
      current_idents[cells_in_cl] <- max_pred
    } else {
      # Leave them as the original cluster label (cl) – no change needed,
      # since current_idents[cells_in_cl] was already cl by initialization.
      next
    }
  }
  
  
  #update seurat identities
  
  seurat_obj = SetIdent(object = seurat_obj,
                       value = factor(current_idents, levels = unique(current_idents)))
  
return(seurat_obj)
}
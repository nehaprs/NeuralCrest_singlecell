#=========================
#find min expression of Sox9 in Sox9 enriched cells, and min expression of Pax3 in Pax3 enriched cells in 
# works on raw data 's seurat object
#==========================

#incomplete. min expression of Sox9 in Sox9+ object is still 0

library(Seurat)

min_expression = function(sraw){
  
  ## Extract raw expression for Sox9 
  raw_expr <- FetchData(sraw, vars = c("Sox9"), layer = "counts")
}
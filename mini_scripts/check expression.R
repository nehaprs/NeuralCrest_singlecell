#troubleshooting: pax3 at E10.5 unnatual values for pax3 expression


setwd("~/BINF/yushi scrnaseq/E10.5/Pax3/3-pax3-10-5_clean_RSEC_MolsPerCell_MEX")
pax105data= ReadMtx(mtx = "matrix.mtx.gz",
                   cells = "barcodes.tsv.gz",
                   features = "features.tsv.gz")
summary(pax105data)
pax105so = CreateSeuratObject(counts = pax105data)
summary(pax105so)

#cehck what expression looks like before normalization
pax_expr = FetchData(pax105so, vars ="Pax3")
# scatterplot expression
pax_expr$cell = rownames(pax_expr)
pax_expr$index = seq_len(nrow(pax_expr))
min(pax_expr$Pax3[pax_expr$Pax3 > 0] )
# Scatter plot

ggplot(pax_expr, aes(x = index, y = Pax3)) +
  geom_point(alpha = 0.6, color = "darkblue") +
  labs(title = "Pax3 Expression in Cells with orig.ident == 'pax9' at E10.5",
       x = "Cell Index",
       y = "Pax9 Expression (log-normalized)") +
  theme_minimal()

pax105so[["percent.mt"]] <- PercentageFeatureSet(pax105so, pattern = "^MT-")

vln = VlnPlot(pax105so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pax105so <- subset(pax105so, subset = nFeature_RNA > 4000 & nFeature_RNA < 8500 & percent.mt < 5)
pax105so <- NormalizeData(pax105so)

#check what expression looks like after normalization
paxnr_expr = FetchData(pax105so, vars ="Pax3")
# scatterplot expression
paxnr_expr$cell = rownames(paxnr_expr)
paxnr_expr$index = seq_len(nrow(paxnr_expr))


ggplot(paxnr_expr, aes(x = index, y = Pax3)) +
  geom_point(alpha = 0.6, color = "darkblue") +
  labs(title = "Pax3 Expression in Cells with orig.ident == 'pax9' at E10.5 after normalization",
       x = "Cell Index",
       y = "Pax9 Expression (log-normalized)") +
  theme_minimal()

##### the original nrl doesn't bring it below 0
## check for the saved RDS object
pax105 <- readRDS("~/BINF/yushi scrnaseq/New folder/rds objects/all clusters/pax105.rds")

#check what expression looks like after normalization
pax2_expr = FetchData(pax105, vars ="Pax3")
# scatterplot expression
pax2_expr$cell = rownames(pax2_expr)
pax2_expr$index = seq_len(nrow(pax2_expr))


ggplot(pax2_expr, aes(x = index, y = Pax3)) +
  geom_point(alpha = 0.6, color = "darkblue") +
  labs(title = "Pax3 Expression in Cells with orig.ident == 'pax9' at E10.5 after normalization, saved file",
       x = "Cell Index",
       y = "Pax9 Expression (log-normalized)") +
  theme_minimal()

#exactly the same plot as done from the scratch, no negative values.

##does it change after normalizing in the list
paxlist = s.list[[2]]

#check what expression looks like after normalization
pax2_expr = FetchData(paxlist, vars ="Pax3")
any(pax2_expr < 0)
#all positive
# scatterplot expression
pax2_expr$cell = rownames(pax2_expr)
pax2_expr$index = seq_len(nrow(pax2_expr))


ggplot(pax2_expr, aes(x = index, y = Pax3)) +
  geom_point(alpha = 0.6, color = "darkblue") +
  labs(title = "Pax3 Expression in Cells with orig.ident == 'pax9' at E10.5 after normalization, saved file",
       x = "Cell Index",
       y = "Pax9 Expression (log-normalized)") +
  theme_minimal()

#processing as list also doesn't make a difference

##check expression in s.combined 
DefaultAssay(s.combined)
pax3Cells = subset(s.combined, orig.ident == "pax3")
pax3_expr = FetchData(s.combined, vars ="Pax3")
any(pax3_expr < 0) 
soxinpax = FetchData(s.combined, vars = "Sox9")
min(soxinpax)
#############################################################

source("~/GitHub/NeuralCrest_singlecell/all six/plot_normalized_expression.R")

 

#convert query feature names from gene symbol to ensembl
#s.query = cc
library(biomaRt)



#Set up the biomaRt connection to the Ensembl database for mouse
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

# Extract the row names (gene symbols) from the Seurat object
mouse_gene_symbols <- rownames(query)
head(rownames(query))
# Map mouse gene symbols to Ensembl IDs using biomaRt
gene_mapping <- getBM(
  attributes = c("mgi_symbol", "ensembl_gene_id"),
  filters = "mgi_symbol",
  values = mouse_gene_symbols,
  mart = ensembl
)

# Match Ensembl IDs to the row names in the Seurat object
# Filter out rows without mapping
mapped_genes <- gene_mapping[match(mouse_gene_symbols, gene_mapping$mgi_symbol), ]
rownames(query) <- mapped_genes$ensembl_gene_id
head(rownames(query))
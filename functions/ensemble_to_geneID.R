library(biomaRt)
#change rownames of seurat object from ensemble ID to gene symbol
#Set up the biomaRt connection to the Ensembl database for mouse


ensemble_toId = function(s.data){
#ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
  ensembl <- useEnsembl(biomart = "genes", dataset = "xtropicalis_gene_ensembl")
# Extract the row names from the Seurat object
ensembl_ids <- rownames(s.data)

# Map  Ensembl IDs to mouse gene symbols to  using biomaRt
gene_mapping <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)

# Match gene symbols to the row names in the Seurat object
# Filter out rows without mapping
mapped_symbols <- gene_mapping[match(ensembl_ids, gene_mapping$ensembl_gene_id), ]

# Replace rownames with gene symbols (keep original Ensembl ID if no match found)
new_row_names <- ifelse(!is.na(mapped_symbols$mgi_symbol), mapped_symbols$mgi_symbol, ensembl_ids)

# Assign new rownames to the Seurat object
rownames(s.data) <- new_row_names
return(s.data)
}

setwd("~/BINF/yushi scrnaseq/New folder/rds objects")
ensemble_dir = file.path(getwd(),"ensemble_id")

# List all .rds files
rds_files = list.files(ensemble_dir, pattern = "rds$")

#loop through RDS files
for(file in rds_files){
  filename = paste0("~/BINF/yushi scrnaseq/New folder/rds objects/ensemble_id/", file)
  obj = readRDS(filename)
  obj_processed = ensemble_toId(obj)
  output_file = file.path(getwd(), basename(file))
  saveRDS(obj_processed, file = output_file)
  
  cat("Processed and saved:", output_file, "\n")
}

#####################################################

#process references
setwd("~/BINF/yushi scrnaseq/New folder/rds objects/geneID")
ref95 <- readRDS("~/BINF/yushi scrnaseq/E9.5/tomeE9.5.rds")
ref105 <- readRDS("~/BINF/yushi scrnaseq/E10.5/tome_E10.5.rds")
ref115 <- readRDS("~/BINF/yushi scrnaseq/E11.5/tomeRef_E11.5.rds")

#change row names to gene id
ref95 = ensemble_toId(ref95)
head(rownames(ref95))
#Renaming features in v3/v4 assays is not supported 

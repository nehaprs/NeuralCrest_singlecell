ConvertHumanGeneListToMM <- function(x){
  # Load human ensembl attributes
  human = biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org")
  # Load mouse ensembl attributes
  mouse = biomaRt::useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl",  host="https://dec2021.archive.ensembl.org")
  # Link both datasets and retrieve mouse genes from the human genes
  genes.list = biomaRt::getLDS(attributes = c("hgnc_symbol"), 
                               filters = "hgnc_symbol", values = x , 
                               mart = human, attributesL = c("mgi_symbol"),
                               martL = mouse, uniqueRows = F)
  # Get unique names of genes (in case gene names are duplicated)
  mouse.gene.list <- unique(genes.list[, 2])
  return(mouse.gene.list)
}
#counting cells

#count cells per cluster for each of the 6 mouse data

setwd("~/BINF/yushi scrnaseq/New folder/cell_counts")

library(dplyr)
library(Seurat)
library(patchwork)
library(writexl)
library(readxl)

s.data =  readRDS("~/BINF/yushi scrnaseq/E11.5/sox9/seurat output/sox9.rds")

cl.counts = as.data.frame(table(Idents(s.data)))
colnames(cl.counts) = c("cell type", "number of cells")
write_xlsx(cl.counts,"sox9E115_counts.xlsx")


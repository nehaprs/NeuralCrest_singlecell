#counting cells

#count cells per cluster for each of the 6 mouse data

setwd("~/BINF/yushi scrnaseq/New folder/cell_counts")

library(dplyr)
library(Seurat)
library(patchwork)
library(writexl)
library(readxl)

#s.data <- readRDS("~/BINF/yushi scrnaseq/E9.5/Sox9/ref_annot/sox95.rds")

cl.counts = as.data.frame(table(Idents(s.data)))
colnames(cl.counts) = c("cell type", "number of cells")
write_xlsx(cl.counts,"sox9E95_countsv2.xlsx")


#count the cell annotations (not the cluster annotation)
setwd("~/BINF/yushi scrnaseq/New folder/cell_counts/cell counts")

s.data <- readRDS("~/BINF/yushi scrnaseq/E11.5/sox9/ref_annot/sox115.rds")

cell.count = as.data.frame(table(s.data$predicted.id))
colnames(cell.count) = c("cell type", "number of cells")
write_xlsx(cell.count, "sox9e115cellcount.xlsx")

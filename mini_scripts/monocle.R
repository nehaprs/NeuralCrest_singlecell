library(stringr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(monocle3)
library(stringr)
library(SeuratWrappers) 


setwd("~/BINF/yushi scrnaseq/time series/harmony_slingshot/paxfull")
s.processed = readRDS("paxFullCombined_procesd.rds")


s.processed[["RNA"]] = JoinLayers(s.processed[["RNA"]])

# Find root cells: E9.5 neural crest
table(s.processed@meta.data$cell_state)
root_cells = WhichCells(s.processed,
                        expression = eday == "E9.5" & str_detect(s.processed@meta.data$predicted.id, regex("Neural crest", ignore_case = TRUE)))


root_clusters = unique(s.processed$cell_state[root_cells])


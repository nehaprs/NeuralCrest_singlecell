#counting cells

#count cells per cluster for each of the 6 mouse data

setwd("~/BINF/yushi scrnaseq/New folder/cell_counts")

library(dplyr)
library(Seurat)
library(patchwork)
library(writexl)
library(readxl)
library(openxlsx)

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

#finding confidence level of cell-level annotation
setwd("~/BINF/yushi scrnaseq/New folder/cell_counts/cell counts/cell count confidence")
s.data <- readRDS("~/BINF/yushi scrnaseq/E9.5/Pax3/ref_annot/pax95.rds") 

predicted_scores = s.data@meta.data[, c("predicted.id", "prediction.score.max")]


confidence_summary = predicted_scores %>%
  group_by(predicted.id) %>%
  summarize(Avg_Pred_Score = mean(prediction.score.max, na.rm = TRUE),
            Median_Pred_score = median(prediction.score.max, na.rm = TRUE),
            Cell_Count = n())

write_xlsx(as.data.frame(confidence_summary),"pax3e95_predscore.xlsx")

#do this for all standalone rds objects
input_dir = "~/BINF/yushi scrnaseq/New folder/rds objects"
rds_files = list.files(path = input_dir, pattern = "rds$", full.names = TRUE)

for(file in rds_files){
  
  file_name = tools::file_path_sans_ext(basename(file))
  s.data = readRDS(file)
  predicted_scores = s.data@meta.data[, c("predicted.id", "prediction.score.max")]
  
  
  confidence_summary = predicted_scores %>%
    group_by(predicted.id) %>%
    summarize(Avg_Pred_Score = mean(prediction.score.max, na.rm = TRUE),
              Median_Pred_score = median(prediction.score.max, na.rm = TRUE),
              Cell_Count = n(),
              .groups = "drop")
  
  output_file = paste0(file_name,"_predscore_cells.xlsx")
  write_xlsx(as.data.frame(confidence_summary),output_file)
}
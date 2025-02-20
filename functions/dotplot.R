#plot NC markers for subcluster annotation
setwd("~/BINF/yushi scrnaseq/New folder/genes_plots")
library(dplyr)
library(Seurat)
library(ggplot2)
library(Matrix)
library(grid)
s.data <- readRDS("~/BINF/yushi scrnaseq/New folder/rds objects/subclusters/pax115.rds")
name = "Pax115"

marker_list = list(premigratory_NC = c("Zic1","Zic3", "Zic5",  "Msx1", "Mafb", "Gdf7", "Sox9",
                                       "Snai1"),
                   Delaminating_NC = c("Sox9", "Snai1","Dlx5", "Pdgfra","Pak3", "Hapln"),
                   Migratory_NC = c("Sox10", "Foxd3", "Ets1", "Sox9"),
                   Sensory_NC = c("Neurog2", "Neurod1", "Isl1", "Stmn2"),
                   Autonomic_NC = c("Phox2b", "Ascl1", "Hand2"),
                   Mesenchymal_NC = c("Prrx1", "Twist1", "Dlx2","Sox9"),
                   Melanoblasts = c("Mitf", "Pmel", "Dct"),
                   Vagal_Cardiac_NC = c("Hand1",  "Hand2",  "Phox2b",  "Ret"),
                   Trunk_NC = c("Neurog2", "Sox10", "Phox2b"),
                   Cranial_NC =c("Twist1", "Prrx2", "Sox9", "Dlx2"),
                   
                   
                   Sensory_Neurons = c("Neurog2", "Pou4f1", "Neurod1", "Neurog1", "Isl1", "Stmn2", "Dcx"),
                   Autonomic_Neurons = c("Phox2b", "Ascl1", "Hand2", "Dbh", "Th", "Ret"),
                   Glial_Cells = c("Sox10", "Mpz", "Plp1", "Fabp7", "Zfp488"), 
                   Melanocytes = c("Mitf", "Pmel", "Dct"),
                   Craniofacial_Mesenchyme = c("Prrx1", "Twist1", "Dlx2", "Sox9", "Runx2", "Col2a1" ),
                   Chondrocytes = c("Sox9", "Col2a1", "Acan" ),
                   Osteoblasts = c("Runx2", "Sp7", "Col1a1" ),
                   Smooth_Muscle_Cells = c("Hand2", "Acta2", "Myh11" ),
                   Cardiac_Mesenchyme = c("Hand1", "Hand2", "Msx2", "Dlx6", "Gata6" ),
                   Adrenal_Chromaffin_Cells = c("Phox2b", "Th", "Dbh", "Ascl1"),
                   Enteric_Neurons = c("Phox2b", "Ret", "Sox10"),
                   Thyroid_Parafollicular_C_Cells = c("Gata3", "Ascl1", "Calca"),
                   Sympathetic_Neurons = c("Phox2b", "Th", "Ascl1", "Hand2"),
                   Parasympathetic_Neurons = c("Phox2b", "Ret", "Hand2"),
                   Corneal_Endothelium_Stroma = c("Pax6", "Pitx2", "Foxc1"),
                   Dental_Mesenchyme = c("Barx1", "Msx1", "Dlx2")
                   )



for(group in names(marker_list)) {
  # Extract the vector of genes for the current group
  genes <- marker_list[[group]]
  
  # Create custom y-axis labels with cluster names and cell counts
  cluster_counts <- table(Idents(s.data))
  new_labels <- paste0(names(cluster_counts), " (", cluster_counts, ")")
  names(new_labels) <- names(cluster_counts)
  
  # Generate the DotPlot without modifying the Seurat object
  p <- DotPlot(s.data, features = genes) +
    scale_y_discrete(labels = new_labels) +
    ggtitle(paste(group, "markers for", name))+
    scale_x_discrete(drop = FALSE)
  
  # Force the x-axis factor to include all genes in the desired order.
  # This sets the factor levels for the 'features.plot' column to 'genes'
  p$data$features.plot <- factor(p$data$features.plot, levels = genes)
  
  # Adjust theme settings to reduce spacing between axis labels
  p <- p + theme(
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1, margin = margin(t = 0)),
    axis.text.y = element_text(size = 16, margin = margin(r = 0)),
    panel.spacing = unit(0.1, "lines"),
    plot.margin = margin(5, 5, 5, 5),
    text = element_text(size = 16),
    aspect.ratio = 1/1
  )
  
  # Construct a filename based on the group name and the object name, and save the plot
  file_name <- paste0(group, "_", name, ".png")
  ggsave(filename = file_name, plot = p, width = 10, height = 5, dpi = 300, bg = "white")
}

#plot NC markers for subcluster annotation
setwd("~/BINF/yushi scrnaseq/New folder/genes_plots")

#s.data <- readRDS("~/BINF/yushi scrnaseq/New folder/rds objects/sox95.rds")
#table(s.data@active.ident)
#keep_clusters = c("Mesenchymal stromal cells", "Pre-epidermal keratinocytes", "Blood progenitors", "Neural crest (PNS neurons)", "Primitive erythroid cells")
# Extract the name from the file path 
#name = tools::file_path_sans_ext(basename("~/BINF/yushi scrnaseq/New folder/rds objects/sox95.rds"))


setwd("~/BINF/yushi scrnaseq/New folder/genes_plots")

#get s.data from makers_list.txt
s.data <- readRDS("~/BINF/yushi scrnaseq/New folder/rds objects/sox115.rds")
name = "sox115"
table(s.data@active.ident)
keep_clusters = c("Limb mesenchyme progenitors", "Neuron progenitor cells", "Inhibitory interneurons", "Neural crest (PNS glia)","Pre-epidermal keratinocytes","Primitive erythroid cells", "Epidermis")
#these clusters are later made and saved in subclustering.R

##

s.data2 = subset(s.data, idents = keep_clusters)


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


# Loop over each group name in marker_list
for(group in names(marker_list)) {
  # Extract the vector of genes for the current group
  genes <- marker_list[[group]]
  
  # Here, you can use your DotPlot code. For example:
  p <- DotPlot(s.data2, features = genes) +
    ggtitle(paste(group, "markers for", name))
  
  # Construct a filename based on the group name and your RDS object's name
  file_name <- paste0(group, "_", name, ".png")
  
  # Save the plot with a white background
  ggsave(filename = file_name, plot = p, width = 10, height = 8, dpi = 300, bg = "white")
}



















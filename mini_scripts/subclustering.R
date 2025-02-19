#selct clusters of interest 

setwd("~/BINF/yushi scrnaseq/New folder/rds objects/subclusters")

s.data <- readRDS("~/BINF/yushi scrnaseq/New folder/rds objects/all clusters/sox95.rds")
#table(s.data@active.ident)
keep_clusters = c("Mesenchymal stromal cells", "Pre-epidermal keratinocytes", "Blood progenitors", "Neural crest (PNS neurons)", "Primitive erythroid cells")
s.data2 = subset(s.data, idents = keep_clusters)
saveRDS(s.data2, "sox95.rds")
rm(s.data, s.data2, keep_clusters)

s.data <- readRDS("~/BINF/yushi scrnaseq/New folder/rds objects/all clusters/pax95.rds")
table(s.data@active.ident)
keep_clusters = c("Pre-epidermal keratinocytes", "Neural crest (PNS glia)", "Mesenchymal stromal cells", "Neural crest (PNS neurons)")
s.data2 = subset(s.data, idents = keep_clusters)
saveRDS(s.data2, "pax95.rds")
rm(s.data, s.data2, keep_clusters)

s.data <- readRDS("~/BINF/yushi scrnaseq/New folder/rds objects/all clusters/pax105.rds")
name = "pax105"
keep_clusters = c("Mesenchymal Stem Cells","Neural Progenitor Cells","Craniofacial Mesenchyme Cells","Neural Crest Cells","Neurons","Erythrocytes","T Lymphocytes","Hematopoietic Stem Cells")
s.data2 = subset(s.data, idents = keep_clusters)
saveRDS(s.data2, "pax105.rds")
rm(s.data, s.data2, keep_clusters)



s.data <- readRDS("~/BINF/yushi scrnaseq/New folder/rds objects/all clusters/sox105.rds")
name = "sox105"
keep_clusters = c("Mesenchymal stromal cells", "Neural crest (PNS glia)", "Pre-epidermal keratinocytes","Neuron progenitor cells")
s.data2 = subset(s.data, idents = keep_clusters)
saveRDS(s.data2, "sox105.rds")
rm(s.data, s.data2, keep_clusters)



s.data <- readRDS("~/BINF/yushi scrnaseq/New folder/rds objects/all clusters/sox115.rds")
name = "sox115"
keep_clusters = c("Limb mesenchyme progenitors", "Neuron progenitor cells", "Inhibitory interneurons", "Neural crest (PNS glia)","Pre-epidermal keratinocytes","Primitive erythroid cells", "Epidermis")
s.data2 = subset(s.data, idents = keep_clusters)
saveRDS(s.data2, "sox115.rds")
rm(s.data, s.data2, keep_clusters)



s.data <- readRDS("~/BINF/yushi scrnaseq/New folder/rds objects/all clusters/pax115.rds")
name = "pax115"
keep_clusters = c("Limb mesenchyme progenitors", "Myocytes","Neuron progenitor cells","Neural crest (PNS glia)","Intermediate progenitor cells","Pre-epidermal keratinocytes","Neural crest (PNS neurons)","Primitive erythroid cells")
s.data2 = subset(s.data, idents = keep_clusters)
saveRDS(s.data2, "pax115.rds")
rm(s.data, s.data2, keep_clusters)


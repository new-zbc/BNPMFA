library(SingleCellExperiment)
library(aricode)
library(spARI)

sampleID = "BZ14"

load(paste0("application/STARmap/", sampleID, "/data/", sampleID, ".RData"))
sce$label = as.character(sce$label)
sce$label[sce$label == "1"] = "L1"
sce$label[sce$label == "2"] = "L2/3"
sce$label[sce$label == "3"] = "L5"
sce$label[sce$label == "4"] = "L6"

load(paste0("application/STARmap/", sampleID, "/cluster_prediction/BNPMFA.RData"))
sce$BNPMFA = pred_label

load(paste0("application/STARmap/", sampleID, "/cluster_prediction/BayesSpace.RData"))
sce$BayesSpace = pred_label

load(paste0("application/STARmap/", sampleID, "/cluster_prediction/SCMEB.RData"))
sce$SCMEB = c(pred_label)

load(paste0("application/STARmap/", sampleID, "/cluster_prediction/DRSC.RData"))
sce$DRSC = pred_label

load(paste0("application/STARmap/", sampleID, "/cluster_prediction/Louvain.RData"))
sce$Louvain = pred_label

load(paste0("application/STARmap/", sampleID, "/cluster_prediction/SpaGCN.RData"))
sce$SpaGCN = pred_label +1

load(paste0("application/STARmap/", sampleID, "/cluster_prediction/BASS.RData"))
sce$BASS = pred_label

# for method adept
load(paste0("application/STARmap/", sampleID, "/cluster_prediction/ADEPT.RData"))
sce$ADEPT = pred_label 


load(paste0("application/STARmap/", sampleID, "/cluster_prediction/GraphST.RData"))
sce$GraphST = pred_label


load(paste0("application/STARmap/", sampleID, "/cluster_prediction/STAGATE.RData"))
sce$STAGATE = pred_label


load(paste0("application/STARmap/", sampleID, "/cluster_prediction/BANKSY.RData"))
sce$BANKSY = pred_label 


library(ggplot2)
source("reproduce/plotFunctions.R")
p1 = plot_STARmap(sce, sampleID = sampleID, Method = "Ground Truth")
p2 = plot_STARmap(sce, Method = "BNPMFA")
p3 = plot_STARmap(sce, Method = "BayesSpace")
p4 = plot_STARmap(sce, Method = "SCMEB")
p5 = plot_STARmap(sce, Method = "DRSC")
p6 = plot_STARmap(sce, Method = "Louvain")
p7 = plot_STARmap(sce, Method = "SpaGCN")
p8 = plot_STARmap(sce, Method = "BASS")
p9 = plot_STARmap(sce,  Method = "STAGATE")
p10 = plot_STARmap(sce, Method = "GraphST")
p11 = plot_STARmap(sce, Method = "ADEPT")
p12 = plot_STARmap(sce, Method = "BANKSY")


library(cowplot)
p <- plot_grid(p1, p2,  p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, byrow = T, nrow = 4, ncol = 3)
ggsave(p, filename =  "reproduce/img/Figure4d.pdf", width = 12, height =8, units = "in", bg = "white")




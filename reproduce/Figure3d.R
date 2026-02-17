library(flexclust)
library(SingleCellExperiment)
library(ggplot2)
library(aricode)
library(spARI)
source("R/plotFunctions.R")


sampleID = 151672

load(paste0("application/DLPFCdata/", sampleID, "/data/", sampleID, "_counts.RData"))

if(sampleID %in% c(151669, 151670, 151671, 151672)){
  sce1$label = as.character(sce1$label)
}

    load(paste0("application/DLPFCdata/", sampleID, "/cluster_prediction/", "BNPMFA.RData"))
    sce1$BNPMFA = pred_label

    load(paste0("application/DLPFCdata/", sampleID, "/cluster_prediction/", "BayesSpace.RData"))
    sce1$BayesSpace = pred_label
    
    load(paste0("application/DLPFCdata/", sampleID, "/cluster_prediction/", "SCMEB.RData"))
    sce1$SCMEB = c(pred_label)
    
    load(paste0("application/DLPFCdata/", sampleID, "/cluster_prediction/", "DRSC.RData"))
    sce1$DRSC = pred_label
    
    load(paste0("application/DLPFCdata/", sampleID, "/cluster_prediction/", "Louvain.RData"))
    sce1$Louvain = as.numeric(pred_label) 
    
    load(paste0("application/DLPFCdata/", sampleID, "/cluster_prediction/", "BASS.RData"))
    sce1$BASS = as.numeric(pred_label) 

    load(paste0("application/DLPFCdata/", sampleID, "/cluster_prediction/", "SpaGCN.RData"))
    #pred_label = pred_label[!is.na()]
    sce1$SpaGCN = pred_label +1
    
    load(paste0("application/DLPFCdata/", sampleID, "/cluster_prediction/", "ADEPT.RData"))
    sce1$ADEPT = pred_label 
    
    load(paste0("application/DLPFCdata/", sampleID, "/cluster_prediction/", "GraphST.RData"))
    sce1$GraphST = pred_label

    load(paste0("application/DLPFCdata/", sampleID, "/cluster_prediction/", "STAGATE.RData"))
    sce1$STAGATE = pred_label

    load(paste0("application/DLPFCdata/", sampleID, "/cluster_prediction/", "BANKSY.RData"))
    sce1$BANKSY = pred_label 
    
    p1 = plot_ST_Visium(sce1, sampleID= sampleID, platform = "Visium", Method = "Ground Truth")
    p2 = plot_ST_Visium(sce1, platform = "Visium", Method = "BNPMFA")
    p3 = plot_ST_Visium(sce1, platform = "Visium", Method = "BayesSpace")
    p4 = plot_ST_Visium(sce1, platform = "Visium", Method = "SCMEB")
    p5 = plot_ST_Visium(sce1, platform = "Visium", Method = "DRSC")
    p6 = plot_ST_Visium(sce1, platform = "Visium", Method = "Louvain")
    p7 = plot_ST_Visium(sce1, platform = "Visium", Method = "SpaGCN")
    p8 = plot_ST_Visium(sce1, platform = "Visium", Method = "BASS")

    p9 = plot_ST_Visium(sce1, platform = "Visium", Method = "STAGATE")
    p10 = plot_ST_Visium(sce1, platform = "Visium", Method = "GraphST")
    p11 = plot_ST_Visium(sce1, platform = "Visium", Method = "ADEPT")
    p12 = plot_ST_Visium(sce1, platform = "Visium", Method = "BANKSY")
    
    library(cowplot)
    p <- plot_grid(p1, p2,  p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, byrow = T, nrow = 3, ncol = 4)

    if(!dir.exists(paste0("reproduce/img"))){
      dir.create(paste0("reproduce/img"))
    }

    ggsave(p, filename = paste0("reproduce/img/Figure3d.jpg"), width = 12, height = 9, units = "in", bg = "white")
    


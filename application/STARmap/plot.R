library(SingleCellExperiment)
library(aricode)
library(spARI)
ARI_mat = matrix(0, nrow = 3, ncol = 11)
AMI_mat = matrix(0, nrow = 3, ncol = 11)
NMI_mat = matrix(0, nrow = 3, ncol = 11)
spARI_mat = matrix(0, nrow = 3, ncol = 11)





sampleID = "BZ5"

load(paste0("STARmap/", "data/", sampleID, ".RData"))
sce$label = as.character(sce$label)
sce$label[sce$label == "1"] = "L1"
sce$label[sce$label == "2"] = "L2/3"
sce$label[sce$label == "3"] = "L5"
sce$label[sce$label == "4"] = "L6"

load(paste0("STARmap/output/", sampleID, "/cluster_prediction/DRMFM.RData"))
sce$BNPMFA = pred_label

load(paste0("STARmap/output/", sampleID, "/cluster_prediction/BayesSpace.RData"))
sce$BayesSpace = pred_label

load(paste0("STARmap/output/", sampleID, "/cluster_prediction/SCMEB.RData"))
sce$SCMEB = c(pred_label)

load(paste0("STARmap/output/", sampleID, "/cluster_prediction/DRSC.RData"))
sce$DRSC = pred_label

load(paste0("STARmap/output/", sampleID, "/cluster_prediction/Louvain.RData"))
sce$Louvain = pred_label

load(paste0("STARmap/output/", sampleID, "/cluster_prediction/SpaGCN.RData"))
sce$SpaGCN = pred_label +1

load(paste0("STARmap/output/", sampleID, "/cluster_prediction/BASS.RData"))
sce$BASS = pred_label

# for method adept
adept = read.csv(paste0("STARmap/output/", sampleID, "/ADEPT/", "5.csv"), row.names=1, header=TRUE)
sce$ADEPT = adept[,1] 

pred_label = sce$ADEPT
save(pred_label, file = paste0("STARmap/output/", sampleID, "/cluster_prediction/ADEPT.RData"))

graphst = read.csv(paste0("STARmap/output/", sampleID, "/GraphST/", "1.csv"), row.names=1, header=TRUE)
sce$GraphST = graphst[,1] + 1
pred_label = sce$GraphST
save(pred_label, file = paste0("STARmap/output/", sampleID, "/cluster_prediction/GraphST.RData"))


stagate = read.csv(paste0("STARmap/output/", sampleID, "/STAGATE/", "1.csv"), row.names=1, header=TRUE)
sce$STAGATE = stagate[,1] + 1
pred_label = sce$STAGATE
save(pred_label, file = paste0("STARmap/output/", sampleID, "/cluster_prediction/STAGATE.RData"))


load(paste0("STARmap/output/", sampleID, "/BANKSY/", "1.RData"))
sce$BANKSY = pred_label 
pred_label = sce$BANKSY
save(pred_label, file = paste0("STARmap/output/", sampleID, "/cluster_prediction/BANKSY.RData"))


for(i in c("BNPMFA", "BayesSpace", "SCMEB", "DRSC", "Louvain", "SpaGCN", "BASS", "STAGATE", "GraphST", "ADEPT", "BANKSY")){
  ARI_mat[1, which(c("BNPMFA", "BayesSpace", "SCMEB", "DRSC", "Louvain", "SpaGCN", "BASS", "STAGATE", "GraphST", "ADEPT", "BANKSY") == i)] = ARI(sce$label, sce[[i]])
  AMI_mat[1, which(c("BNPMFA", "BayesSpace", "SCMEB", "DRSC", "Louvain", "SpaGCN", "BASS", "STAGATE", "GraphST", "ADEPT", "BANKSY") == i)] = AMI(sce$label, sce[[i]])
  NMI_mat[1, which(c("BNPMFA", "BayesSpace", "SCMEB", "DRSC", "Louvain", "SpaGCN", "BASS", "STAGATE", "GraphST", "ADEPT", "BANKSY") == i)] = NMI(sce$label, sce[[i]])
  spARI_mat[1, which(c("BNPMFA", "BayesSpace", "SCMEB", "DRSC", "Louvain", "SpaGCN", "BASS", "STAGATE", "GraphST", "ADEPT", "BANKSY") == i)] = spARI(sce$label, sce[[i]], coords = as.data.frame(colData(sce)[, c("row", "col")]))[2]
}


library(ggplot2)
source("plot/plotFunctions.R")
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

# p1 = plot_STARmap_spARI(sce, Method = "Ground Truth")
# p2 = plot_STARmap_spARI(sce, Method = "BNPMFA")
# p3 = plot_STARmap_spARI(sce, Method = "BayesSpace")
# p4 = plot_STARmap_spARI(sce, Method = "SCMEB")
# p5 = plot_STARmap_spARI(sce, Method = "DRSC")
# p6 = plot_STARmap_spARI(sce, Method = "Louvain")
# p7 = plot_STARmap_spARI(sce, Method = "SpaGCN")
# p8 = plot_STARmap_spARI(sce, Method = "BASS")


library(cowplot)
p <- plot_grid(p1, p2,  p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, byrow = T, nrow = 4, ncol = 3)

if(!dir.exists(paste0("STARmap/img/", sampleID))){
  dir.create(paste0("STARmap/img/", sampleID))
}

ggsave(p, filename = paste0("STARmap/img/", sampleID, "_cluster.jpg"), width =12, height = 8, units = "in", bg = "white")
ggsave(p, filename =  paste0("STARmap/img/", sampleID, "_cluster.pdf"), width = 12, height =8, units = "in", bg = "white")









library(SingleCellExperiment)
sampleID = "BZ14"

load(paste0("STARmap/", "data/", sampleID, ".RData"))
sce$label = as.character(sce$label)
sce$label[sce$label == "1"] = "L1"
sce$label[sce$label == "2"] = "L2/3"
sce$label[sce$label == "3"] = "L5"
sce$label[sce$label == "4"] = "L6"

load(paste0("STARmap/output/", sampleID, "/cluster_prediction/DRMFM.RData"))
sce$BNPMFA = pred_label

load(paste0("STARmap/output/", sampleID, "/cluster_prediction/BayesSpace.RData"))
sce$BayesSpace = pred_label

load(paste0("STARmap/output/", sampleID, "/cluster_prediction/SCMEB.RData"))
sce$SCMEB = c(pred_label)

load(paste0("STARmap/output/", sampleID, "/cluster_prediction/DRSC.RData"))
sce$DRSC = pred_label

load(paste0("STARmap/output/", sampleID, "/cluster_prediction/Louvain.RData"))
sce$Louvain = pred_label

load(paste0("STARmap/output/", sampleID, "/cluster_prediction/SpaGCN.RData"))
sce$SpaGCN = pred_label +1

load(paste0("STARmap/output/", sampleID, "/cluster_prediction/BASS.RData"))
sce$BASS = pred_label

# for method adept
adept = read.csv(paste0("STARmap/output/", sampleID, "/ADEPT/", "1.csv"), row.names=1, header=TRUE)
sce$ADEPT = adept[,1] 
pred_label = sce$ADEPT
save(pred_label, file = paste0("STARmap/output/", sampleID, "/cluster_prediction/ADEPT.RData"))


graphst = read.csv(paste0("STARmap/output/", sampleID, "/GraphST/", "1.csv"), row.names=1, header=TRUE)
sce$GraphST = graphst[,1] + 1
pred_label = sce$GraphST
save(pred_label, file = paste0("STARmap/output/", sampleID, "/cluster_prediction/GraphST.RData"))


stagate = read.csv(paste0("STARmap/output/", sampleID, "/STAGATE/", "2.csv"), row.names=1, header=TRUE)
sce$STAGATE = stagate[,1] + 1
pred_label = sce$STAGATE
save(pred_label, file = paste0("STARmap/output/", sampleID, "/cluster_prediction/STAGATE.RData"))


load(paste0("STARmap/output/", sampleID, "/BANKSY/", "1.RData"))
#pred_label = pred_label[!is.na()]
sce$BANKSY = pred_label 
pred_label = sce$BANKSY
save(pred_label, file = paste0("STARmap/output/", sampleID, "/cluster_prediction/BANKSY.RData"))



library(ggplot2)
source("plot/plotFunctions.R")
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

# p1 = plot_STARmap_spARI(sce, Method = "Ground Truth")
# p2 = plot_STARmap_spARI(sce, Method = "BNPMFA")
# p3 = plot_STARmap_spARI(sce, Method = "BayesSpace")
# p4 = plot_STARmap_spARI(sce, Method = "SCMEB")
# p5 = plot_STARmap_spARI(sce, Method = "DRSC")
# p6 = plot_STARmap_spARI(sce, Method = "Louvain")
# p7 = plot_STARmap_spARI(sce, Method = "SpaGCN")
# p8 = plot_STARmap_spARI(sce, Method = "BASS")

for(i in c("BNPMFA", "BayesSpace", "SCMEB", "DRSC", "Louvain", "SpaGCN", "BASS", "STAGATE", "GraphST", "ADEPT", "BANKSY")){
  ARI_mat[2, which(c("BNPMFA", "BayesSpace", "SCMEB", "DRSC", "Louvain", "SpaGCN", "BASS", "STAGATE", "GraphST", "ADEPT", "BANKSY") == i)] = ARI(sce$label, sce[[i]])
  AMI_mat[2, which(c("BNPMFA", "BayesSpace", "SCMEB", "DRSC", "Louvain", "SpaGCN", "BASS", "STAGATE", "GraphST", "ADEPT", "BANKSY") == i)] = AMI(sce$label, sce[[i]])
  NMI_mat[2, which(c("BNPMFA", "BayesSpace", "SCMEB", "DRSC", "Louvain", "SpaGCN", "BASS", "STAGATE", "GraphST", "ADEPT", "BANKSY") == i)] = NMI(sce$label, sce[[i]])
  spARI_mat[2, which(c("BNPMFA", "BayesSpace", "SCMEB", "DRSC", "Louvain", "SpaGCN", "BASS", "STAGATE", "GraphST", "ADEPT", "BANKSY") == i)] = spARI(sce$label, sce[[i]], coords = as.data.frame(colData(sce)[, c("row", "col")]))[2]
}




library(cowplot)
p <- plot_grid(p1, p2,  p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, byrow = T, nrow = 4, ncol = 3)

if(!dir.exists(paste0("STARmap/img/", sampleID))){
  dir.create(paste0("STARmap/img/", sampleID))
}

ggsave(p, filename = paste0("STARmap/img/", sampleID, "_cluster.jpg"), width =12, height = 8, units = "in", bg = "white")
ggsave(p, filename =  paste0("STARmap/img/", sampleID, "_cluster.pdf"), width = 12, height =8, units = "in", bg = "white")







library(SingleCellExperiment)
sampleID = "BZ9"

load(paste0("STARmap/", "data/", sampleID, ".RData"))
sce$label = as.character(sce$label)
sce$label[sce$label == "1"] = "L1"
sce$label[sce$label == "2"] = "L2/3"
sce$label[sce$label == "3"] = "L5"
sce$label[sce$label == "4"] = "L6"

# load(paste0("STARmap/output/", sampleID, "/cluster_prediction/f_3.RData"))
# #sce$BNPMFA = result$pred_label
# sce$BNPMFA = result$MCMCList$group_iter[, 100] +1

load(paste0("STARmap/cov_1/", sampleID, "/seed_2/f_3.5.RData"))
sce$BNPMFA = result$pred_label
pred_label = result$pred_label
save(pred_label, file = paste0("STARmap/output/", sampleID, "/cluster_prediction/BNPMFA.RData"))


load(paste0("STARmap/output/", sampleID, "/cluster_prediction/BayesSpace.RData"))
sce$BayesSpace = pred_label

load(paste0("STARmap/output/", sampleID, "/cluster_prediction/SCMEB.RData"))
sce$SCMEB = c(pred_label)

load(paste0("STARmap/output/", sampleID, "/cluster_prediction/DRSC.RData"))
sce$DRSC = pred_label

load(paste0("STARmap/output/", sampleID, "/cluster_prediction/Louvain.RData"))
sce$Louvain = pred_label

load(paste0("STARmap/output/", sampleID, "/cluster_prediction/SpaGCN.RData"))
sce$SpaGCN = pred_label +1

load(paste0("STARmap/output/", sampleID, "/cluster_prediction/BASS.RData"))
sce$BASS = pred_label

# for method adept
adept = read.csv(paste0("STARmap/output/", sampleID, "/ADEPT/", "1.csv"), row.names=1, header=TRUE)
sce$ADEPT = adept[,1] 
pred_label = sce$ADEPT
save(pred_label, file = paste0("STARmap/output/", sampleID, "/cluster_prediction/ADEPT.RData"))


graphst = read.csv(paste0("STARmap/output/", sampleID, "/GraphST/", "1.csv"), row.names=1, header=TRUE)
sce$GraphST = graphst[,1] + 1
pred_label = sce$GraphST
save(pred_label, file = paste0("STARmap/output/", sampleID, "/cluster_prediction/GraphST.RData"))


stagate = read.csv(paste0("STARmap/output/", sampleID, "/STAGATE/", "2.csv"), row.names=1, header=TRUE)
sce$STAGATE = stagate[,1] + 1
pred_label = sce$STAGATE
save(pred_label, file = paste0("STARmap/output/", sampleID, "/cluster_prediction/STAGATE.RData"))


load(paste0("STARmap/output/", sampleID, "/BANKSY/", "1.RData"))
#pred_label = pred_label[!is.na()]
sce$BANKSY = pred_label 
pred_label = sce$BANKSY
save(pred_label, file = paste0("STARmap/output/", sampleID, "/cluster_prediction/BANKSY.RData"))



library(ggplot2)
source("plot/plotFunctions.R")
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

# p1 = plot_STARmap_spARI(sce, Method = "Ground Truth")
# p2 = plot_STARmap_spARI(sce, Method = "BNPMFA")
# p3 = plot_STARmap_spARI(sce, Method = "BayesSpace")
# p4 = plot_STARmap_spARI(sce, Method = "SCMEB")
# p5 = plot_STARmap_spARI(sce, Method = "DRSC")
# p6 = plot_STARmap_spARI(sce, Method = "Louvain")
# p7 = plot_STARmap_spARI(sce, Method = "SpaGCN")
# p8 = plot_STARmap_spARI(sce, Method = "BASS")


for(i in c("BNPMFA", "BayesSpace", "SCMEB", "DRSC", "Louvain", "SpaGCN", "BASS", "STAGATE", "GraphST", "ADEPT", "BANKSY")){
  ARI_mat[3, which(c("BNPMFA", "BayesSpace", "SCMEB", "DRSC", "Louvain", "SpaGCN", "BASS", "STAGATE", "GraphST", "ADEPT", "BANKSY") == i)] = ARI(sce$label, sce[[i]])
  AMI_mat[3, which(c("BNPMFA", "BayesSpace", "SCMEB", "DRSC", "Louvain", "SpaGCN", "BASS", "STAGATE", "GraphST", "ADEPT", "BANKSY") == i)] = AMI(sce$label, sce[[i]])
  NMI_mat[3, which(c("BNPMFA", "BayesSpace", "SCMEB", "DRSC", "Louvain", "SpaGCN", "BASS", "STAGATE", "GraphST", "ADEPT", "BANKSY") == i)] = NMI(sce$label, sce[[i]])
  spARI_mat[3, which(c("BNPMFA", "BayesSpace", "SCMEB", "DRSC", "Louvain", "SpaGCN", "BASS", "STAGATE", "GraphST", "ADEPT", "BANKSY") == i)] = spARI(sce$label, sce[[i]], coords = as.data.frame(colData(sce)[, c("row", "col")]))[2]
}




library(cowplot)
p <- plot_grid(p1, p2,  p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, byrow = T, nrow = 4, ncol = 3)

if(!dir.exists(paste0("STARmap/img/", sampleID))){
  dir.create(paste0("STARmap/img/", sampleID))
}

ggsave(p, filename = paste0("STARmap/img/", sampleID, "_cluster.jpg"), width =12, height = 8, units = "in", bg = "white")
ggsave(p, filename =  paste0("STARmap/img/", sampleID, "_cluster.pdf"), width = 12, height = 8, units = "in", bg = "white")









#############################################
#
#
## plot ARI
#
#
###############################################

# ARI = matrix(0, nrow = 3, ncol = 7)
# ARI[1, ] = c(0.767, 0.214, 0.880, 0.248, 0.186, 0.299, 0.152)
# ARI[2, ] = c(0.611, 0.195, 0.564, 0.277, 0.209, 0.309, 0.185)
# ARI[3, ] = c(0.820, 0.216, 0.761, 0.278, 0.258, 0.347, 0.171)

ARI = as.data.frame(ARI_mat)
colnames(ARI) = c("BNPMFA", "BayesSpace", "SCMEB", "DRSC", "Louvain", "SpaGCN", "BASS", "STAGATE", "GraphST", "ADEPT", "BANKSY")
rownames(ARI) = c("BZ5", "BZ14", "BZ9")
write.csv(ARI, file = "STARmap/output/ARI_summary.csv")
library(reshape2)
ARI_long = melt(ARI)
colnames(ARI_long) = c("Method", "ARI")
ARI_long$Method = factor(ARI_long$Method, levels = c("BNPMFA", "BASS", "STAGATE", "BANKSY", "ADEPT", "GraphST" , "BayesSpace", "SpaGCN", "SCMEB", "DRSC", "Louvain" ))

ARI_mean = aggregate(ARI_long[,2], list(ARI_long$Method), mean)
colnames(ARI_mean) = c("Method", "ARI")
# library(ggplot2)
# plot_box = ggplot(data = ARI_long, aes(x=Method, y=ARI, fill = Method)) + #geom_violin(trim = F) +
#   geom_boxplot(width = 0.2, outlier.shape = NA) +
#   geom_hline(yintercept = mean(ARI[,1]), color = "red", linetype = "dashed", linewidth=1) + 
#   theme_bw() +
#   #ggtitle(paste0("Mean Neighborhood Distance")) +
#   theme(
#     #plot.title = element_text(face="bold", hjust = 0.5, size = 10),
#     axis.title.y = element_text(size = 12,face = "bold"),
#     axis.title.x = element_blank(),
#     legend.title = element_text(size = 10,face = "bold"),
#     axis.text = element_text(size = 10, face = "bold", angle = 45, hjust=1, vjust=1),
#     #axis.text = element_text(size = 8, face = "bold"),
#     legend.position="none") + ylim(c(0,1))



library(ggplot2)
plot_box = ggplot(data = ARI_mean, aes(x=Method, y=ARI, fill = Method)) + #geom_violin(trim = F) +
  geom_bar(stat = "identity") +
  #geom_hline(yintercept = mean(ARI[,1]), color = "red", linetype = "dashed", linewidth=1) + 
  theme_bw() +
  #ggtitle(paste0("Mean Neighborhood Distance")) +
  theme(
    #plot.title = element_text(face="bold", hjust = 0.5, size = 10),
    axis.title.y = element_text(size = 12,face = "bold"),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 10,face = "bold"),
    axis.text = element_text(size = 10, face = "bold", angle = 45, hjust=1, vjust=1),
    #axis.text = element_text(size = 8, face = "bold"),
    legend.position="none") + ylim(c(0,1))

ggsave(plot_box, filename = paste0("STARmap/img/", "ARI.jpg"), width = 4, height = 4, units = "in")
ggsave(plot_box, filename =  paste0("STARmap/img/",  "ARI.pdf"), width = 4, height = 4, units = "in")



#############################################
#
#
## plot AMI
#
#
###############################################

AMI = as.data.frame(AMI_mat)
colnames(AMI) = c("BNPMFA", "BayesSpace", "SCMEB", "DRSC", "Louvain", "SpaGCN", "BASS", "STAGATE", "GraphST", "ADEPT", "BANKSY")
rownames(AMI) = c("BZ5", "BZ14", "BZ9")
write.csv(AMI, file = "STARmap/output/AMI_summary.csv")
library(reshape2)
AMI_long = melt(AMI)
colnames(AMI_long) = c("Method", "AMI")
AMI_long$Method = factor(AMI_long$Method, levels = c("BNPMFA", "BASS", "STAGATE", "BANKSY", "ADEPT", "GraphST" , "BayesSpace", "SpaGCN", "SCMEB", "DRSC", "Louvain" ))

AMI_mean = aggregate(AMI_long[,2], list(AMI_long$Method), mean)
colnames(AMI_mean) = c("Method", "AMI")

library(ggplot2)
plot_box1 = ggplot(data = AMI_mean, aes(x=Method, y=AMI, fill = Method)) + #geom_violin(trim = F) +
  geom_bar(stat = "identity") +
  #geom_hline(yintercept = mean(ARI[,1]), color = "red", linetype = "dashed", linewidth=1) + 
  theme_bw() +
  #ggtitle(paste0("Mean Neighborhood Distance")) +
  theme(
    #plot.title = element_text(face="bold", hjust = 0.5, size = 10),
    axis.title.y = element_text(size = 12,face = "bold"),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 10,face = "bold"),
    axis.text = element_text(size = 10, face = "bold", angle = 45, hjust=1, vjust=1),
    #axis.text = element_text(size = 8, face = "bold"),
    legend.position="none") + ylim(c(0,1))

ggsave(plot_box1, filename = paste0("STARmap/img/", "AMI.jpg"), width = 4, height = 4, units = "in")
ggsave(plot_box1, filename =  paste0("STARmap/img/",  "AMI.pdf"), width = 4, height = 4, units = "in")



#############################################
#
#
## plot NMI
#
#
###############################################

NMI = as.data.frame(NMI_mat)
colnames(NMI) = c("BNPMFA", "BayesSpace", "SCMEB", "DRSC", "Louvain", "SpaGCN", "BASS", "STAGATE", "GraphST", "ADEPT", "BANKSY")
rownames(NMI) = c("BZ5", "BZ14", "BZ9")
write.csv(NMI, file = "STARmap/output/NMI_summary.csv")
library(reshape2)
NMI_long = melt(NMI)
colnames(NMI_long) = c("Method", "NMI")
NMI_long$Method = factor(NMI_long$Method, levels = c("BNPMFA", "BASS", "STAGATE", "BANKSY", "ADEPT", "GraphST" , "BayesSpace", "SpaGCN", "SCMEB", "DRSC", "Louvain" ))

NMI_mean = aggregate(NMI_long[,2], list(NMI_long$Method), mean)
colnames(NMI_mean) = c("Method", "NMI")

library(ggplot2)
plot_box2 = ggplot(data = NMI_mean, aes(x=Method, y=NMI, fill = Method)) + #geom_violin(trim = F) +
  geom_bar(stat = "identity") +
  #geom_hline(yintercept = mean(ARI[,1]), color = "red", linetype = "dashed", linewidth=1) + 
  theme_bw() +
  #ggtitle(paste0("Mean Neighborhood Distance")) +
  theme(
    #plot.title = element_text(face="bold", hjust = 0.5, size = 10),
    axis.title.y = element_text(size = 12,face = "bold"),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 10,face = "bold"),
    axis.text = element_text(size = 10, face = "bold", angle = 45, hjust=1, vjust=1),
    #axis.text = element_text(size = 8, face = "bold"),
    legend.position="none") + ylim(c(0,1))

ggsave(plot_box2, filename = paste0("STARmap/img/", "NMI.jpg"), width = 4, height = 4, units = "in")
ggsave(plot_box2, filename =  paste0("STARmap/img/",  "NMI.pdf"), width = 4, height = 4, units = "in")


library(cowplot)
p <- plot_grid(plot_box2, plot_box1, byrow = T, nrow = 1, ncol = 2)
ggsave(p, filename = paste0("STARmap/img/", "AMI_NMI.jpg"), width = 8, height = 4, units = "in")
ggsave(p, filename =  paste0("STARmap/img/",  "AMI_NMI.pdf"), width = 8, height = 4, units = "in")

#############################################
#
#
## plot spARI
#
#
###############################################

spARI = as.data.frame(spARI_mat)
ARI = as.data.frame(spARI)
colnames(ARI) = c("BNPMFA", "BayesSpace", "SCMEB", "DRSC", "Louvain", "SpaGCN", "BASS", "STAGATE", "GraphST", "ADEPT", "BANKSY")
library(reshape2)
ARI_long = melt(ARI)
colnames(ARI_long) = c("Method", "spARI")
ARI_long$Method = factor(ARI_long$Method, levels = c("BNPMFA", "BASS", "STAGATE", "BANKSY", "ADEPT", "GraphST" , "BayesSpace", "SpaGCN", "SCMEB", "DRSC", "Louvain" ))

ARI_mean = aggregate(ARI_long[,2], list(ARI_long$Method), mean)
colnames(ARI_mean) = c("Method", "spARI")
# library(ggplot2)
# plot_box = ggplot(data = ARI_long, aes(x=Method, y=ARI, fill = Method)) + #geom_violin(trim = F) +
#   geom_boxplot(width = 0.2, outlier.shape = NA) +
#   geom_hline(yintercept = mean(ARI[,1]), color = "red", linetype = "dashed", linewidth=1) + 
#   theme_bw() +
#   #ggtitle(paste0("Mean Neighborhood Distance")) +
#   theme(
#     #plot.title = element_text(face="bold", hjust = 0.5, size = 10),
#     axis.title.y = element_text(size = 12,face = "bold"),
#     axis.title.x = element_blank(),
#     legend.title = element_text(size = 10,face = "bold"),
#     axis.text = element_text(size = 10, face = "bold", angle = 45, hjust=1, vjust=1),
#     #axis.text = element_text(size = 8, face = "bold"),
#     legend.position="none") + ylim(c(0,1))



library(ggplot2)
plot_box = ggplot(data = ARI_mean, aes(x=Method, y=spARI, fill = Method)) + #geom_violin(trim = F) +
  geom_bar(stat = "identity") +
  #geom_hline(yintercept = mean(ARI[,1]), color = "red", linetype = "dashed", linewidth=1) + 
  theme_bw() +
  #ggtitle(paste0("Mean Neighborhood Distance")) +
  theme(
    #plot.title = element_text(face="bold", hjust = 0.5, size = 10),
    axis.title.y = element_text(size = 12,face = "bold"),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 10,face = "bold"),
    axis.text = element_text(size = 10, face = "bold", angle = 45, hjust=1, vjust=1),
    #axis.text = element_text(size = 8, face = "bold"),
    legend.position="none") + ylim(c(0,1))

ggsave(plot_box, filename = paste0("STARmap/img/", "spARI.jpg"), width = 4, height = 4, units = "in")
ggsave(plot_box, filename =  paste0("STARmap/img/",  "spARI.pdf"), width = 4, height = 4, units = "in")



# #############################################
# #
# #
# ## Differential gene analysis
# #
# #
# ###############################################
# library(SingleCellExperiment)
# sampleID = "BZ14"

# load(paste0("STARmap/", "data/", sampleID, ".RData"))
# sce$label = as.character(sce$label)
# sce$label[sce$label == "1"] = "L1"
# sce$label[sce$label == "2"] = "L2/3"
# sce$label[sce$label == "3"] = "L5"
# sce$label[sce$label == "4"] = "L6"

# load(paste0("STARmap/", sampleID, "/others/DRMFM.RData"))
# sce$BNPMFA = pred_label

# #construc Seurat object
# library(Seurat)
# sce = scater::logNormCounts(sce)

# seu = as.Seurat(sce, counts = "counts", data = "logcounts", project = "sce_to_seurat")


# seu = ScaleData(seu)
# Idents(seu) = as.factor(sce$BNPMFA)
# #rownames(seu@assays$originalexp@scale.data) = seu@assays$originalexp@meta.features$gene_name
# #seu@assays$originalexp@data@Dimnames[[1]] = seu@assays$originalexp@meta.features$gene_name
# #seu@assays$originalexp@counts@Dimnames[[1]] = seu@assays$originalexp@meta.features$gene_name

# #names(seu@active.ident) = seu@assays$originalexp@counts@Dimnames[[2]]


# seu.markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# library(dplyr)
# top_gene <- seu.markers %>%
#   group_by(cluster) %>% top_n(n=20, avg_log2FC)

# write.csv(top_gene, file = paste0("STARmap/", sampleID, "/differential_gene.csv"))

# library(ggplot2)
# p = DoHeatmap(seu, features = top_gene$gene,
#               group.bar = T, slot = "scale.data", size = 6,
#               label = T, draw.lines = T, combine = T) + 
#   theme(axis.text.y = element_text(face='bold'),
#         legend.text = element_text(size = 8),
#         legend.title = element_text( size = 8, face='bold', angle = 90)) +
#   guides(colour = guide_colorbar("title")) +
#   #scale_fill_gradientn(colors = c("blue", "white", "red"))+
#   guides(fill = guide_colorbar(title = "Relative expression", title.position = "left"))

# p

# ggsave(p, filename = paste0("STARmap/", sampleID, "/heatmap.jpg"), width = 6, height = 8, units = "in")
# ggsave(p, filename =  paste0("STARmap/", sampleID, "/heatmap.pdf"), width = 6, height = 8, units = "in")

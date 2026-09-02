library(Rcpp)
library(SingleCellExperiment)
library(RcppArmadillo)
library(flexclust)


sampleID = 151509


filename = paste0("application/DLPFCdata/", sampleID, "/data/", sampleID, "_counts.RData")
load(file = filename)

source("R/utils.R")

Adj = find_neighbors(sce1, "Visium", "lattice")
neighbors = find_neighbor_index(Adj, "Visium")

set.seed(0)
sce = spatialPreprocess(sce1, n.PCs = 15, n.HVGs = n_gene)
index = which(rowData(sce)$is.HVG)

load(paste0("application/DLPFCdata/", sampleID, "/low_embeddings.rda"))


pca_dist = rep(0, dim(pca_est)[1])
for(i in 1:dim(pca_est)[1]){
  neighbor_index = neighbors[i, which(neighbors[i, ] != 0)]
  if(length(neighbor_index) > 1){
    pca_dist[i] = mean(rowSums((pca_est[neighbor_index, ] - pca_est[i, ])^2))
  } 
}



BNPMFA_dist = rep(0, dim(bnpmfa_est)[1])
for(i in 1:dim(bnpmfa_est)[1]){
  neighbor_index = neighbors[i, which(neighbors[i, ] != 0)]
  if(length(neighbor_index) > 1){
    BNPMFA_dist[i] = mean(rowSums((bnpmfa_est[neighbor_index, ] - bnpmfa_est[i, ])^2))
  }
}



neighbors_spa = neighbors  
spadist = rep(0, dim(spa_est)[1])
for(i in 1:dim(spa_est)[1]){
  neighbor_index = neighbors_spa[i, which(neighbors_spa[i, ] != 0)]
  if(length(neighbor_index) > 1){
    spadist[i] = mean(rowSums((spa_est[neighbor_index, ] - spa_est[i, ])^2))
  }
}



df_dist = data.frame(Method = c(rep("PCA", length(pca_dist)), 
                                rep("BNPMFA", length(BNPMFA_dist)), 
                                rep("SpatialPCA", length(spadist))),
                     MND = c(pca_dist, BNPMFA_dist, spadist))

df_dist$Method = factor(df_dist$Method, levels = c("BNPMFA", "PCA", "SpatialPCA"))

# wilcox.test(BNPMFA_dist, pca_dist, paired = T)
# wilcox.test(BNPMFA_dist, spadist, paired = T)
# wilcox.test(spadist, pca_dist, paired = T)
# median(BNPMFA_dist - pca_dist)
# median(spadist - pca_dist)


baseline1 = 11
baseline2 = baseline1 + 5
baseline3 = baseline1 + 2.5
library(latex2exp)
library(ggplot2)
df_dist$MND[which(df_dist$MND > 10)] = NA
mean(is.na(df_dist$MND))
plot_box1 = ggplot(data = df_dist, aes(x=Method, y=MND, color = Method)) + geom_violin(trim = T) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  theme_bw() +
  ggtitle(paste0("Mean Neighborhood Distance")) +
  theme(
    plot.title = element_text(face="bold", hjust = 0.5, size = 10),
    axis.title.y = element_text(size = 10,face = "bold"),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 10,face = "bold"),
    #axis.text = element_text(size = 8, face = "bold", angle = 45, hjust=1, vjust=1),
    axis.text = element_text(size = 8, face = "bold"),
    legend.position="none") + ylim(c(-1, 17)) + 
  annotate("segment", x=1, xend=2, y=baseline1, yend=baseline1)+
  annotate("segment", x=1, xend=1, y=baseline1-0.5, yend=baseline1)+
  annotate("segment", x=2, xend=2, y=baseline1-0.5, yend=baseline1)+
  annotate("text", x=1.5, y=baseline1+0.8, label="p < 0.0001", size=3) +
  annotate("text", x=1.5, y=baseline1-0.8, label="diff = -0.34", size=3) +
  
  annotate("segment", x=1, xend=3, y=baseline2, yend=baseline2)+
  annotate("segment", x=1, xend=1, y=baseline2-0.5, yend=baseline2)+
  annotate("segment", x=3, xend=3, y=baseline2-0.5, yend=baseline2)+
  annotate("text", x=2, y=baseline2+0.8, label="p = 0.1893", size=3) +
  
  annotate("segment", x=2, xend=3, y=baseline3, yend=baseline3)+
  annotate("segment", x=2, xend=2, y=baseline3-0.5, yend=baseline3)+
  annotate("segment", x=3, xend=3, y=baseline3-0.5, yend=baseline3)+
  annotate("text", x=2.5, y=baseline3+0.8, label="p < 0.0001", size=3)+
  annotate("text", x=2.5, y=baseline3-0.8, label="diff = -0.46", size=3)









label = sce1$label

pca_dist = rep(0, dim(pca_est)[1])
for(i in 1:length(pca_dist)){
  label_i = label[i]
  index = which(label != label_i)
  pca_dist[i] = mean(rowSums((pca_est[index, ] - pca_est[i, ])^2))
}

BNPMFA_dist = rep(0, dim(bnpmfa_est)[1])
for(i in 1:dim(bnpmfa_est)[1]){
  label_i = label[i]
  index = which(label != label_i)
  BNPMFA_dist[i] = mean(rowSums((bnpmfa_est[index, ] - bnpmfa_est[i, ])^2))
}



spadist = rep(0, dim(spa_est)[1])
for(i in 1:dim(spa_est)[1]){
  label_i = label[i]
  index = which(label != label_i)
  spadist[i] = mean(rowSums((spa_est[index, ] - spa_est[i, ])^2))
}




df_dist = data.frame(Method = c(rep("PCA", length(pca_dist)),
                                rep("BNPMFA", length(BNPMFA_dist)), 
                                rep("SpatialPCA", length(spadist))),
                     MOD = c(pca_dist, BNPMFA_dist, spadist))

df_dist$Method = factor(df_dist$Method, levels = c("BNPMFA", "PCA", "SpatialPCA"))

wilcox.test(BNPMFA_dist, pca_dist, paired = T)
wilcox.test(BNPMFA_dist, spadist, paired = T)
wilcox.test(spadist, pca_dist, paired = T)
median(BNPMFA_dist - pca_dist)
median(spadist - pca_dist)


baseline1 = 11
baseline2 = baseline1 + 5
baseline3 = baseline1 + 2.5
library(latex2exp)
library(ggplot2)
df_dist$MOD[which(df_dist$MOD > 10)] = NA
mean(is.na(df_dist$MOD))
plot_box2 = ggplot(data = df_dist, aes(x=Method, y=MOD, color = Method)) + geom_violin(trim = T) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  theme_bw() +
  ggtitle(paste0("Mean Out-of-domain Distance")) +
  theme(
    plot.title = element_text(face="bold", hjust = 0.5, size = 10),
    axis.title.y = element_text(size = 10,face = "bold"),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 10,face = "bold"),
    #axis.text = element_text(size = 8, face = "bold", angle = 45, hjust=1, vjust=1),
    axis.text = element_text(size = 8, face = "bold"),
    legend.position="none") + ylim(c(0, 17)) + 
  annotate("segment", x=1, xend=2, y=baseline1, yend=baseline1)+
  annotate("segment", x=1, xend=1, y=baseline1-0.5, yend=baseline1)+
  annotate("segment", x=2, xend=2, y=baseline1-0.5, yend=baseline1)+
  annotate("text", x=1.5, y=baseline1+0.8, label="p < 0.0001", size=3) +
  annotate("text", x=1.5, y=baseline1-0.8, label="diff = 0.11", size=3) +
  
  annotate("segment", x=1, xend=3, y=baseline2, yend=baseline2)+
  annotate("segment", x=1, xend=1, y=baseline2-0.5, yend=baseline2)+
  annotate("segment", x=3, xend=3, y=baseline2-0.5, yend=baseline2)+
  annotate("text", x=2, y=baseline2+0.8, label="p < 0.0001", size=3) +
  
  annotate("segment", x=2, xend=3, y=baseline3, yend=baseline3)+
  annotate("segment", x=2, xend=2, y=baseline3-0.5, yend=baseline3)+
  annotate("segment", x=3, xend=3, y=baseline3-0.5, yend=baseline3)+
  annotate("text", x=2.5, y=baseline3+0.8, label="p < 0.0001", size=3)+
  annotate("text", x=2.5, y=baseline3-0.8, label="diff = -0.19", size=3)




p_dist_151509 = cowplot::plot_grid(plot_box1,  plot_box2, ncol = 2)
#ggsave(p_dist_151509, filename = paste0("DLPFCdata/", sampleID, "/img/distance_compare.pdf"), width = 6, height = 3.2, units = "in")
ggsave(p_dist_151509, filename = paste0("reproduce/img/Figure3b.jpg"), width = 6, height = 3.2, units = "in")


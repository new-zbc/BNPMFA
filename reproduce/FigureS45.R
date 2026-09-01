library(Rcpp)
library(SingleCellExperiment)
library(RcppArmadillo)
library(flexclust)
color_pal = c("#1F77B4FF", "#FF7F0EFF", "#D62728FF", "#2CA02CFF", "#9467BDFF",
              "#8C564BFF", "#E377C2FF", "#7F7F7FFF", "#BCBD22FF", "#17BECFFF",
              "#AEC7E8FF", "#FFBB78FF", "#98DF8AFF", "#FF9896FF", "#C5B0D5FF",
              "#C49C94FF", "#F7B6D2FF", "#C7C7C7FF", "#DBDB8DFF", "#9EDAE5FF"
)


load("application/VisiumHD/data/sce.RData")
load("application/VisiumHD/data/var.RData")
rowData(sce) <- vars
rownames(sce)<- vars$gene_ids
colnames(sce) <- paste0("cell_", 1:dim(sce)[2])

library(ggplot2)
source("reproduce/plotFunctions.R")

load("application/VisiumHD/BNPMFA.RData")
sce$label = pred_label
sce$label_tumor = "Tumor"
sce$label_tumor[sce$label != 7] = "Non-tumor"

plot_dim = function(position, label, title = "PCA Embedding", size = 0.1){
  
  library(ggplot2)
  df = data.frame(x = position[, 1], y = position[, 2], label = label)
  ggplot(data = df) + geom_point(aes(x=x, y=y, color = label), size = size) + 
    theme_bw() + ggtitle(title) + #xlab("Dim 1") + ylab("Dim 2")+ 
    xlim(c(-3,3)) + ylim(c(-3,3)) + scale_color_manual(values = color_pal[1:length(unique(label))])+
    theme(
      legend.position = "none",
      plot.title = element_text(face="bold", hjust = 0.5, size = 12),
      axis.title = element_text(size = 10,face = "bold"),
      legend.title = element_text(size = 10,face = "bold"),
      axis.text = element_text(size = 8, face = "bold")
    ) + xlab("Dimension 1") + ylab("Dimension 2") +
    stat_ellipse(aes(x = position[, 1], y = position[, 2],color = label),
                 #alpha = 0.2,
                 type = "norm",
                 show.legend = FALSE,
                 level = 0.85)
}

load("application/VisiumHD/embedding.RData")

load("application/VisiumHD/spatialPCA_embeddings.rda")
source("R/utils.R")
spa_est = t(spaPCA[1:2, ])
outlier = setdiff(colnames(sce), colnames(spaPCA)) 
outlier_index = which(colnames(sce) %in% outlier)
### normalize
spa_est[, 1] = (spa_est[, 1] - mean(spa_est[, 1])) / sd(spa_est[, 1])
spa_est[, 2] = (spa_est[, 2] - mean(spa_est[, 2])) / sd(spa_est[, 2])

if(length(outlier_index) > 0){
  sce_spa = sce[, -outlier_index]
  Adj = find_neighbors(sce_spa, "ST", "lattice")
  neighbors_spa = find_neighbor_index(Adj, "ST")
}else{
  neighbors_spa = neighbors  
}



p1 = plot_dim(pca_est[,1:2], sce$label_tumor, title = "PCA Embedding") 
p2 = plot_dim(Yhat[,1:2], sce$label_tumor, title = "BNPMFA Embedding")
p3 = plot_dim(spa_est[,1:2], sce_spa$label_tumor, title = "SpatialPCA Embedding")
p = cowplot::plot_grid(p2, p1, p3,  ncol = 3)


pca_dist = rep(0, dim(pca_est)[1])
for(i in 1:dim(pca_est)[1]){
  neighbor_index = neighbors[i, which(neighbors[i, ] != 0)]
  if(length(neighbor_index) > 1){
    pca_dist[i] = mean(rowSums((pca_est[neighbor_index, ] - pca_est[i, ])^2))
  }
}


BNPMFA_dist = rep(0, dim(Yhat)[1])
for(i in 1:dim(Yhat)[1]){
  neighbor_index = neighbors[i, which(neighbors[i, ] != 0)]
  if(length(neighbor_index) > 1){
    BNPMFA_dist[i] = mean(rowSums((Yhat[neighbor_index, ] - Yhat[i, ])^2))
  }
}


spa_dist = rep(0, dim(spa_est)[1])
for(i in 1:dim(spa_est)[1]){
  neighbor_index = neighbors_spa[i, which(neighbors_spa[i, ] != 0)]
  if(length(neighbor_index) > 1){
    spa_dist[i] = mean(rowSums((spa_est[neighbor_index, ] - spa_est[i, ])^2))
  }
}


df_dist = data.frame(Method = c(rep("PCA", length(pca_dist)), 
                                rep("BNPMFA", length(BNPMFA_dist)), 
                                rep("SpatialPCA", length(spa_dist))),
                     MND = c(pca_dist, BNPMFA_dist, spa_dist))

df_dist$Method = factor(df_dist$Method, levels = c("BNPMFA", "PCA", "SpatialPCA"))

library(latex2exp)
library(ggplot2)
df_dist$MND[which(df_dist$MND > 10)] = NA
baseline1 = 11
baseline2 = baseline1 + 5
baseline3 = baseline1 + 2.5
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
  annotate("text", x=1.5, y=baseline1-0.8, label="diff = -0.24", size=3) +
  
  annotate("segment", x=1, xend=3, y=baseline2, yend=baseline2)+
  annotate("segment", x=1, xend=1, y=baseline2-0.5, yend=baseline2)+
  annotate("segment", x=3, xend=3, y=baseline2-0.5, yend=baseline2)+
  annotate("text", x=2, y=baseline2+0.8, label="p < 0.0001", size=3) +
  
  annotate("segment", x=2, xend=3, y=baseline3, yend=baseline3)+
  annotate("segment", x=2, xend=2, y=baseline3-0.5, yend=baseline3)+
  annotate("segment", x=3, xend=3, y=baseline3-0.5, yend=baseline3)+
  annotate("text", x=2.5, y=baseline3+0.8, label="p < 0.0001", size=3) +
  annotate("text", x=2.5, y=baseline3-0.8, label="diff = -0.41", size=3)




label = sce$label_tumor
pca_dist = rep(0, dim(pca_est)[1])
for(i in 1:length(pca_dist)){
  label_i = label[i]
  index = which(label != label_i)
  pca_dist[i] = mean(rowSums((pca_est[index, ] - pca_est[i, ])^2))
}

BNPMFA_dist = rep(0, dim(Yhat)[1])
for(i in 1:dim(Yhat)[1]){
  label_i = label[i]
  index = which(label != label_i)
  BNPMFA_dist[i] = mean(rowSums((Yhat[index, ] - Yhat[i, ])^2))
}

label_spa = sce_spa$label_tumor
spa_dist = rep(0, dim(spa_est)[1])
for(i in 1:dim(spa_est)[1]){
  label_i = label_spa[i]
  index = which(label_spa != label_i)
  spa_dist[i] = mean(rowSums((spa_est[index, ] - spa_est[i, ])^2))
}


df_dist = data.frame(Method = c(rep("PCA", length(pca_dist)), 
                                rep("BNPMFA", length(BNPMFA_dist)), 
                                rep("SpatialPCA", length(spa_dist))),
                     MOD = c(pca_dist, BNPMFA_dist, spa_dist))

df_dist$Method = factor(df_dist$Method, levels = c("BNPMFA", "PCA", "SpatialPCA"))
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
  annotate("text", x=1.5, y=baseline1-0.8, label="diff = 0.27", size=3) +
  
  annotate("segment", x=1, xend=3, y=baseline2, yend=baseline2)+
  annotate("segment", x=1, xend=1, y=baseline2-0.5, yend=baseline2)+
  annotate("segment", x=3, xend=3, y=baseline2-0.5, yend=baseline2)+
  annotate("text", x=2, y=baseline2+0.8, label="p < 0.0001", size=3) +
  
  annotate("segment", x=2, xend=3, y=baseline3, yend=baseline3)+
  annotate("segment", x=2, xend=2, y=baseline3-0.5, yend=baseline3)+
  annotate("segment", x=3, xend=3, y=baseline3-0.5, yend=baseline3)+
  annotate("text", x=2.5, y=baseline3+0.8, label="p < 0.0001", size=3)+
  annotate("text", x=2.5, y=baseline3-0.8, label="diff = -0.18", size=3)



p = cowplot::plot_grid(p2, p1, p3, plot_box1, plot_box2,  ncol = 3)
ggsave(p, filename = "reproduce/img/FigureS45.pdf", width = 9, height = 6, units = "in", bg = "white")


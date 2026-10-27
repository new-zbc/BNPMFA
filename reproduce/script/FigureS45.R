library(Rcpp)
library(SingleCellExperiment)
library(RcppArmadillo)
library(flexclust)
library(mclust)


color_pal = c("#1F77B4FF", "#FF7F0EFF", "#D62728FF", "#2CA02CFF", "#9467BDFF",
              "#8C564BFF", "#E377C2FF", "#7F7F7FFF", "#BCBD22FF", "#17BECFFF",
              "#AEC7E8FF", "#FFBB78FF", "#98DF8AFF", "#FF9896FF", "#C5B0D5FF",
              "#C49C94FF", "#F7B6D2FF", "#C7C7C7FF", "#DBDB8DFF", "#9EDAE5FF"
)


plot_dim = function(position, label, title = "PCA Embedding", size = 1){
  
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



sampleID = "BZ5"
load(paste0("application/STARmap/", sampleID, "/data/", sampleID, ".RData"))
sce$label = as.character(sce$label)
sce$label[sce$label == "1"] = "L1"
sce$label[sce$label == "2"] = "L2/3"
sce$label[sce$label == "3"] = "L5"
sce$label[sce$label == "4"] = "L6"

load(paste0("application/STARmap/", sampleID, "/", sampleID, "_embedding.RData"))
p2 = plot_dim(pca_est, sce$label, title = "PCA Embedding") 
p3 = plot_dim(spaPCA_est, sce$label[-outlier_index], title = "SpatialPCA Embedding") 
p1 = plot_dim(bnpmfa_est, sce$label, title = "BNPMFA Embedding") 



source("R/utils.R")
temp <- data.frame(id = 1:dim(sce)[2], x = colData(sce)$col, y = colData(sce)$row)

tt = voronoi_adjacency(data = temp, id~x+y, scale=1, PLOT=FALSE)
Adj = tt$G   ## a n-by-n Adjacency Matrix, 1 if two samples are neighbors; 0 otherwise
neighbors = tt$P   ## a n-by-m matrix, m: maximum # of neighbors changing from data to data, 



BNPMFA_dist = rep(0, dim(bnpmfa_est)[1])
BNPMFA_back = rep(0, dim(bnpmfa_est)[1])
for(i in 1:dim(bnpmfa_est)[1]){
  neighbor_index = neighbors[i, which(neighbors[i, ] != 0)]
  if(length(neighbor_index) > 1){
    BNPMFA_dist[i] = mean(rowSums((bnpmfa_est[neighbor_index, ] - bnpmfa_est[i, ])^2))
    BNPMFA_back[i] = mean(rowSums((bnpmfa_est[-neighbor_index, ] - bnpmfa_est[i, ])^2))
  }
}


PCA_dist = rep(0, dim(pca_est)[1])
PCA_back = rep(0, dim(pca_est)[1])
for(i in 1:dim(pca_est)[1]){
  neighbor_index = neighbors[i, which(neighbors[i, ] != 0)]
  if(length(neighbor_index) > 1){
    PCA_dist[i] = mean(rowSums((pca_est[neighbor_index, ] - pca_est[i, ])^2))
    PCA_back[i] = mean(rowSums((pca_est[-neighbor_index, ] - pca_est[i, ])^2))
  }
}


if(length(outlier_index) > 0){
  temp <- data.frame(id = 1:dim(sce)[2], x = colData(sce)$col, y = colData(sce)$row)
  temp <- temp[-outlier_index, ]
  tt = voronoi_adjacency(data = temp, id~x+y, scale=1, PLOT=FALSE)
  Adj = tt$G   ## a n-by-n Adjacency Matrix, 1 if two samples are neighbors; 0 otherwise
  neighbors_spa = tt$P 
}else{
  neighbors_spa = neighbors  
}


spaPCA_dist = rep(0, dim(spaPCA_est)[1])
spaPCA_back = rep(0, dim(spaPCA_est)[1])
for(i in 1:dim(spaPCA_est)[1]){
  neighbor_index = neighbors_spa[i, which(neighbors_spa[i, ] != 0)]
  if(length(neighbor_index) > 1){
    spaPCA_dist[i] = mean(rowSums((spaPCA_est[neighbor_index, ] - spaPCA_est[i, ])^2))
    spaPCA_back[i] = mean(rowSums((spaPCA_est[-neighbor_index, ] - spaPCA_est[i, ])^2))
  }
}


df_dist = data.frame(Method = c(rep("PCA", length(PCA_dist)),
                                rep("BNPMFA", length(BNPMFA_dist)), 
                                rep("SpatialPCA", length(spaPCA_dist))),
                     MND = c(PCA_dist, BNPMFA_dist, spaPCA_dist))

df_dist$Method = factor(df_dist$Method, levels = c("BNPMFA", "PCA", "SpatialPCA"))

wilcox.test(BNPMFA_dist, PCA_dist, paired = T)
wilcox.test(BNPMFA_dist[-outlier_index], spaPCA_dist, paired = T)
wilcox.test(spaPCA_dist, PCA_dist[-outlier_index], paired = T)

median(BNPMFA_dist - PCA_dist)
median(spaPCA_dist - PCA_dist[-outlier_index])

#library(latex2exp)
library(ggplot2)
plot_box = ggplot(data = df_dist, aes(x=Method, y=MND, color = Method)) + geom_violin(trim = T) +
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
    legend.position="none")+
  ylim(c(-1, 20)) 

plot_BZ5 = cowplot::plot_grid(p1, p2, p3, plot_box,  ncol = 4)






sampleID = "BZ9"
load(paste0("application/STARmap/", sampleID, "/data/", sampleID, ".RData"))
sce$label = as.character(sce$label)
sce$label[sce$label == "1"] = "L1"
sce$label[sce$label == "2"] = "L2/3"
sce$label[sce$label == "3"] = "L5"
sce$label[sce$label == "4"] = "L6"

load(paste0("application/STARmap/", sampleID, "/", sampleID, "_embedding.RData"))
p2 = plot_dim(pca_est, sce$label, title = "PCA Embedding") 
p3 = plot_dim(spaPCA_est, sce$label[-outlier_index], title = "SpatialPCA Embedding") 
p1 = plot_dim(bnpmfa_est, sce$label, title = "BNPMFA Embedding") 



source("R/utils.R")
temp <- data.frame(id = 1:dim(sce)[2], x = colData(sce)$col, y = colData(sce)$row)

tt = voronoi_adjacency(data = temp, id~x+y, scale=1, PLOT=FALSE)
Adj = tt$G   ## a n-by-n Adjacency Matrix, 1 if two samples are neighbors; 0 otherwise
neighbors = tt$P   ## a n-by-m matrix, m: maximum # of neighbors changing from data to data, 



BNPMFA_dist = rep(0, dim(bnpmfa_est)[1])
BNPMFA_back = rep(0, dim(bnpmfa_est)[1])
for(i in 1:dim(bnpmfa_est)[1]){
  neighbor_index = neighbors[i, which(neighbors[i, ] != 0)]
  if(length(neighbor_index) > 1){
    BNPMFA_dist[i] = mean(rowSums((bnpmfa_est[neighbor_index, ] - bnpmfa_est[i, ])^2))
    BNPMFA_back[i] = mean(rowSums((bnpmfa_est[-neighbor_index, ] - bnpmfa_est[i, ])^2))
  }
}


PCA_dist = rep(0, dim(pca_est)[1])
PCA_back = rep(0, dim(pca_est)[1])
for(i in 1:dim(pca_est)[1]){
  neighbor_index = neighbors[i, which(neighbors[i, ] != 0)]
  if(length(neighbor_index) > 1){
    PCA_dist[i] = mean(rowSums((pca_est[neighbor_index, ] - pca_est[i, ])^2))
    PCA_back[i] = mean(rowSums((pca_est[-neighbor_index, ] - pca_est[i, ])^2))
  }
}


if(length(outlier_index) > 0){
  temp <- data.frame(id = 1:dim(sce)[2], x = colData(sce)$col, y = colData(sce)$row)
  temp <- temp[-outlier_index, ]
  tt = voronoi_adjacency(data = temp, id~x+y, scale=1, PLOT=FALSE)
  Adj = tt$G   ## a n-by-n Adjacency Matrix, 1 if two samples are neighbors; 0 otherwise
  neighbors_spa = tt$P 
}else{
  neighbors_spa = neighbors  
}


spaPCA_dist = rep(0, dim(spaPCA_est)[1])
spaPCA_back = rep(0, dim(spaPCA_est)[1])
for(i in 1:dim(spaPCA_est)[1]){
  neighbor_index = neighbors_spa[i, which(neighbors_spa[i, ] != 0)]
  if(length(neighbor_index) > 1){
    spaPCA_dist[i] = mean(rowSums((spaPCA_est[neighbor_index, ] - spaPCA_est[i, ])^2))
    spaPCA_back[i] = mean(rowSums((spaPCA_est[-neighbor_index, ] - spaPCA_est[i, ])^2))
  }
}


df_dist = data.frame(Method = c(rep("PCA", length(PCA_dist)),
                                rep("BNPMFA", length(BNPMFA_dist)), 
                                rep("SpatialPCA", length(spaPCA_dist))),
                     MND = c(PCA_dist, BNPMFA_dist, spaPCA_dist))

df_dist$Method = factor(df_dist$Method, levels = c("BNPMFA", "PCA", "SpatialPCA"))

wilcox.test(BNPMFA_dist, PCA_dist, paired = T)
wilcox.test(BNPMFA_dist[-outlier_index], spaPCA_dist, paired = T)
wilcox.test(spaPCA_dist, PCA_dist[-outlier_index], paired = T)

median(BNPMFA_dist - PCA_dist)
median(spaPCA_dist - PCA_dist[-outlier_index])

#library(latex2exp)
library(ggplot2)
plot_box = ggplot(data = df_dist, aes(x=Method, y=MND, color = Method)) + geom_violin(trim = T) +
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
    legend.position="none")+
  ylim(c(-1, 20)) 

plot_BZ9 = cowplot::plot_grid(p1, p2, p3, plot_box,  ncol = 4)




sampleID = "BZ14"
load(paste0("application/STARmap/", sampleID, "/data/", sampleID, ".RData"))
sce$label = as.character(sce$label)
sce$label[sce$label == "1"] = "L1"
sce$label[sce$label == "2"] = "L2/3"
sce$label[sce$label == "3"] = "L5"
sce$label[sce$label == "4"] = "L6"

load(paste0("application/STARmap/", sampleID, "/", sampleID, "_embedding.RData"))
p2 = plot_dim(pca_est, sce$label, title = "PCA Embedding") 
p3 = plot_dim(spaPCA_est, sce$label, title = "SpatialPCA Embedding") 
p1 = plot_dim(bnpmfa_est, sce$label, title = "BNPMFA Embedding") 



source("R/utils.R")
temp <- data.frame(id = 1:dim(sce)[2], x = colData(sce)$col, y = colData(sce)$row)

tt = voronoi_adjacency(data = temp, id~x+y, scale=1, PLOT=FALSE)
Adj = tt$G   ## a n-by-n Adjacency Matrix, 1 if two samples are neighbors; 0 otherwise
neighbors = tt$P   ## a n-by-m matrix, m: maximum # of neighbors changing from data to data, 



BNPMFA_dist = rep(0, dim(bnpmfa_est)[1])
BNPMFA_back = rep(0, dim(bnpmfa_est)[1])
for(i in 1:dim(bnpmfa_est)[1]){
  neighbor_index = neighbors[i, which(neighbors[i, ] != 0)]
  if(length(neighbor_index) > 1){
    BNPMFA_dist[i] = mean(rowSums((bnpmfa_est[neighbor_index, ] - bnpmfa_est[i, ])^2))
    BNPMFA_back[i] = mean(rowSums((bnpmfa_est[-neighbor_index, ] - bnpmfa_est[i, ])^2))
  }
}


PCA_dist = rep(0, dim(pca_est)[1])
PCA_back = rep(0, dim(pca_est)[1])
for(i in 1:dim(pca_est)[1]){
  neighbor_index = neighbors[i, which(neighbors[i, ] != 0)]
  if(length(neighbor_index) > 1){
    PCA_dist[i] = mean(rowSums((pca_est[neighbor_index, ] - pca_est[i, ])^2))
    PCA_back[i] = mean(rowSums((pca_est[-neighbor_index, ] - pca_est[i, ])^2))
  }
}



neighbors_spa = neighbors  



spaPCA_dist = rep(0, dim(spaPCA_est)[1])
spaPCA_back = rep(0, dim(spaPCA_est)[1])
for(i in 1:dim(spaPCA_est)[1]){
  neighbor_index = neighbors_spa[i, which(neighbors_spa[i, ] != 0)]
  if(length(neighbor_index) > 1){
    spaPCA_dist[i] = mean(rowSums((spaPCA_est[neighbor_index, ] - spaPCA_est[i, ])^2))
    spaPCA_back[i] = mean(rowSums((spaPCA_est[-neighbor_index, ] - spaPCA_est[i, ])^2))
  }
}


df_dist = data.frame(Method = c(rep("PCA", length(PCA_dist)),
                                rep("BNPMFA", length(BNPMFA_dist)), 
                                rep("SpatialPCA", length(spaPCA_dist))),
                     MND = c(PCA_dist, BNPMFA_dist, spaPCA_dist))

df_dist$Method = factor(df_dist$Method, levels = c("BNPMFA", "PCA", "SpatialPCA"))

wilcox.test(BNPMFA_dist, PCA_dist, paired = T)
wilcox.test(BNPMFA_dist, spaPCA_dist, paired = T)
wilcox.test(spaPCA_dist, PCA_dist, paired = T)

median(BNPMFA_dist - PCA_dist)
median(spaPCA_dist - PCA_dist)

#library(latex2exp)
library(ggplot2)
plot_box = ggplot(data = df_dist, aes(x=Method, y=MND, color = Method)) + geom_violin(trim = T) +
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
    legend.position="none")+
  ylim(c(-1, 20)) 

plot_BZ14 = cowplot::plot_grid(p1, p2, p3, plot_box,  ncol = 4)

p = cowplot::plot_grid(plot_BZ5, plot_BZ9, plot_BZ14, ncol = 1, labels = c("(a)", "(b)", "(c)"))
ggsave(p, filename = "reproduce/img/FigureS45.pdf", width = 12, height = 9.6, units = "in")



library(Rcpp)
library(SingleCellExperiment)
library(RcppArmadillo)
library(flexclust)
source("reproduce/plotFunctions.R")
plot_dim = function(position, label, title = "PCA Embedding", size = 0.5){
  
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


color_pal = c("#1F77B4FF", "#FF7F0EFF", "#D62728FF", "#2CA02CFF", "#9467BDFF",
              "#8C564BFF", "#E377C2FF", "#7F7F7FFF", "#BCBD22FF", "#17BECFFF",
              "#AEC7E8FF", "#FFBB78FF", "#98DF8AFF", "#FF9896FF", "#C5B0D5FF",
              "#C49C94FF", "#F7B6D2FF", "#C7C7C7FF", "#DBDB8DFF", "#9EDAE5FF"
)

folder = "application/DLPFCdata"

sampleID = 151507
load(paste0(folder, "/", sampleID, "/data/", sampleID, "_counts.RData"))
load(paste0(folder, "/", sampleID, "/", "embedding.RData"))
p_domain = plot_ST_Visium(sce1, sampleID= sampleID, platform = "Visium", Method = "Ground Truth")
p1 = plot_dim(pca_est, sce1$label, title = "PCA Embedding") 
p2 = plot_dim(spa_est, sce1$label, title = "spatialPCA Embedding") 
p3 = plot_dim(Yhat, sce1$label, title = "BNPMFA Embedding") 

p507 = cowplot::plot_grid(p_domain, p3, p1, p2,  ncol = 4)


sampleID = 151508
load(paste0(folder, "/", sampleID, "/data/", sampleID, "_counts.RData"))
load(paste0(folder, "/", sampleID, "/", "embedding.RData"))
p_domain = plot_ST_Visium(sce1, sampleID= sampleID, platform = "Visium", Method = "Ground Truth")
p1 = plot_dim(pca_est, sce1$label, title = "PCA Embedding") 
p2 = plot_dim(spa_est, sce1$label, title = "spatialPCA Embedding") 
p3 = plot_dim(Yhat, sce1$label, title = "BNPMFA Embedding") 

p508 = cowplot::plot_grid(p_domain, p3, p1, p2,  ncol = 4)


sampleID = 151509
load(paste0(folder, "/", sampleID, "/data/", sampleID, "_counts.RData"))
load(paste0(folder, "/", sampleID, "/", "embedding.RData"))
p_domain = plot_ST_Visium(sce1, sampleID= sampleID, platform = "Visium", Method = "Ground Truth")
p1 = plot_dim(pca_est, sce1$label, title = "PCA Embedding") 
p2 = plot_dim(spa_est, sce1$label, title = "spatialPCA Embedding") 
p3 = plot_dim(Yhat, sce1$label, title = "BNPMFA Embedding") 

p509 = cowplot::plot_grid(p_domain, p3, p1, p2,  ncol = 4)

sampleID = 151510
load(paste0(folder, "/", sampleID, "/data/", sampleID, "_counts.RData"))
load(paste0(folder, "/", sampleID, "/", "embedding.RData"))
p_domain = plot_ST_Visium(sce1, sampleID= sampleID, platform = "Visium", Method = "Ground Truth")
p1 = plot_dim(pca_est, sce1$label, title = "PCA Embedding") 
p2 = plot_dim(spa_est, sce1$label, title = "spatialPCA Embedding") 
p3 = plot_dim(Yhat, sce1$label, title = "BNPMFA Embedding") 

p510 = cowplot::plot_grid(p_domain, p3, p1, p2,  ncol = 4)

p = cowplot::plot_grid(p507, p508, p509, p510,  nrow = 4, labels = c("(a)", "(b)", "(c)", "(d)"))
ggsave(p, filename = paste0("reproduce/img/FigureS23.jpg"), width = 12, height = 12, units = "in")


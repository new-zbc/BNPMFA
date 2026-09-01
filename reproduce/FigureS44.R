library(SingleCellExperiment)
compute_eigens <- function(X){
  n <- nrow(X)
  p <- ncol(X)
  corMat <- cor(X)
  evalues <- eigen(corMat)$values
  return(evalues)
}

plot_eigens <- function(sampleID, X, threshold = 0){
  library(ggplot2)
  qmax = 15
  out = compute_eigens(X)[1:qmax]
  df = data.frame(iters=1:qmax, var = out, threshold = 1/((1:qmax)*sum(1/1:(dim(X)[1]))))
  df_new = df
  for(i in 1:(qmax-1)){
    df_new[i, 2] = df[i, 2] - df[i+1, 2]
  }
  df_new = df_new[1:(qmax-1), ]
  index = which(df_new$var - df_new$threshold < 0.02)[1] -1
  p = ggplot(df_new, aes(x = iters, y = var)) + 
    geom_line() + geom_point() + theme_bw() + xlab("Latent dimensions")+
    geom_line(aes(x = iters, y = threshold), linetype = "dashed", color = "red") + 
    ylab("Eigenvalue difference") + labs(title = sampleID)+ 
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 10,face = "bold"),
          legend.title = element_text(size = 10,face = "bold"),
          axis.text = element_text(size = 8, face = "bold"),
          plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 
  p + annotate("segment", x = df_new$iters[index], xend = df_new$iters[index], y = df_new$var[index]+0.2*max(df_new$var), 
               yend = df_new$var[index]+0.1, colour = "blue", size = 0.5, arrow = arrow(angle = 30, length = unit(0.2, "cm")))
}


source("R/utils.R")

sampleID = "BZ5"
filename = paste0("application/STARmap/",sampleID, "/data/", sampleID, ".RData")
load(file = filename)
sce = spatialPreprocess(sce, n.PCs = 5, n.HVGs = dim(sce)[1])
X = as.matrix(t(assay(sce, "logcounts")))
p1 = plot_eigens(sampleID, X)

sampleID = "BZ9"
filename = paste0("application/STARmap/",sampleID, "/data/", sampleID, ".RData")
load(file = filename)
sce = spatialPreprocess(sce, n.PCs = 5, n.HVGs = dim(sce)[1])
X = as.matrix(t(assay(sce, "logcounts")))
p2 = plot_eigens(sampleID, X)


sampleID = "BZ14"
filename = paste0("application/STARmap/",sampleID, "/data/", sampleID, ".RData")
load(file = filename)
sce = spatialPreprocess(sce, n.PCs = 5, n.HVGs = dim(sce)[1])
X = as.matrix(t(assay(sce, "logcounts")))
p3 = plot_eigens(sampleID, X)





library(ggplot2)
data_file = "BZ5"
df = read.csv(paste0("application/STARmap/",sampleID, "/", "latent_dims_explore.csv"))
p1.1 = ggplot(df, aes(x = eps, y = ARI)) + 
  geom_line()+geom_point() + theme_bw() + xlab("Latent dimensions")+
  ylab("ARI") + labs(title = data_file)+ ylim(c(-0.1,1))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 


library(ggplot2)
data_file = "BZ9"
df = read.csv(paste0("application/STARmap/",sampleID, "/", "latent_dims_explore.csv"))
p1.2 = ggplot(df, aes(x = eps, y = ARI)) + 
  geom_line()+geom_point() + theme_bw() + xlab("Latent dimensions")+
  ylab("ARI") + labs(title = data_file)+ ylim(c(-0.1,1))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 


library(ggplot2)
data_file = "BZ14"
df = read.csv(paste0("application/STARmap/",sampleID, "/", "latent_dims_explore.csv"))
p1.3 = ggplot(df, aes(x = eps, y = ARI)) + 
  geom_line()+geom_point() + theme_bw() + xlab("Latent dimensions")+
  ylab("ARI") + labs(title = data_file)+ ylim(c(-0.1,1))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 



library(cowplot)
p <- plot_grid(p1, p2,  p3, p1.1, p1.2, p1.3, byrow = T, nrow = 2, ncol = 3, labels = c("(a)", "", "",
                                                                                        "(b)", "", ""))

ggsave(p, filename = "reproduce/img/FigureS44.pdf", width = 10, height = 6, bg = "white")


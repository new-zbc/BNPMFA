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
  qmax = 21
  out = compute_eigens(X)[1:qmax]
  df = data.frame(iters=1:qmax, var = out, threshold = 1/((1:qmax)*sum(1/1:(dim(X)[1]))))
  df_new = df
  for(i in 1:(qmax-1)){
    df_new[i, 2] = df[i, 2] - df[i+1, 2]
  }
  df_new = df_new[1:(qmax-1), ]
  index = which(df_new$var - df_new$threshold < 0.01)[1] -1
  p = ggplot(df_new, aes(x = iters, y = var)) + 
    geom_line() + geom_point() + theme_bw() + xlab("Latent dimensions")+
    geom_line(aes(x = iters, y = threshold), linetype = "dashed", color = "red") + 
    ylab("Eigenvalue difference") + labs(title = sampleID)+ 
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 10,face = "bold"),
          legend.title = element_text(size = 10,face = "bold"),
          axis.text = element_text(size = 8, face = "bold"),
          plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 
  p + annotate("segment", x = df_new$iters[index], xend = df_new$iters[index], y = df_new$var[index]+2, 
               yend = df_new$var[index]+0.3, colour = "blue", size = 0.5, arrow = arrow(angle = 30, length = unit(0.2, "cm")))
}


source("R/utils.R")

sampleID = 151507
filename = paste0("application/DLPFCdata/", sampleID, "/data/", sampleID, "_counts.RData")
load(file = filename)
sce2 = spatialPreprocess(sce1, n.PCs = 15, n.HVGs = 2000)
index = which(rowData(sce2)$is.HVG)
X = as.matrix(t(assay(sce2[index, ], "logcounts")))
p1 = plot_eigens(sampleID, X)


sampleID = 151508
filename = paste0("application/DLPFCdata/", sampleID, "/data/", sampleID, "_counts.RData")
load(file = filename)
sce2 = spatialPreprocess(sce1, n.PCs = 15, n.HVGs = 2000)
index = which(rowData(sce2)$is.HVG)
X = as.matrix(t(assay(sce2[index, ], "logcounts")))
p2 = plot_eigens(sampleID, X)


sampleID = 151509
filename = paste0("application/DLPFCdata/", sampleID, "/data/", sampleID, "_counts.RData")
load(file = filename)
sce2 = spatialPreprocess(sce1, n.PCs = 15, n.HVGs = 2000)
index = which(rowData(sce2)$is.HVG)
X = as.matrix(t(assay(sce2[index, ], "logcounts")))
p3 = plot_eigens(sampleID, X)


sampleID = 151510
filename = paste0("application/DLPFCdata/", sampleID, "/data/", sampleID, "_counts.RData")
load(file = filename)
sce2 = spatialPreprocess(sce1, n.PCs = 15, n.HVGs = 2000)
index = which(rowData(sce2)$is.HVG)
X = as.matrix(t(assay(sce2[index, ], "logcounts")))
p4 = plot_eigens(sampleID, X)


sampleID = 151669
filename = paste0("application/DLPFCdata/", sampleID, "/data/", sampleID, "_counts.RData")
load(file = filename)
sce2 = spatialPreprocess(sce1, n.PCs = 15, n.HVGs = 2000)
index = which(rowData(sce2)$is.HVG)
X = as.matrix(t(assay(sce2[index, ], "logcounts")))
p5 = plot_eigens(sampleID, X)


sampleID = 151670
filename = paste0("application/DLPFCdata/", sampleID, "/data/", sampleID, "_counts.RData")
load(file = filename)
sce2 = spatialPreprocess(sce1, n.PCs = 15, n.HVGs = 2000)
index = which(rowData(sce2)$is.HVG)
X = as.matrix(t(assay(sce2[index, ], "logcounts")))
p6 = plot_eigens(sampleID, X)



sampleID = 151671
filename = paste0("application/DLPFCdata/", sampleID, "/data/", sampleID, "_counts.RData")
load(file = filename)
sce2 = spatialPreprocess(sce1, n.PCs = 15, n.HVGs = 2000)
index = which(rowData(sce2)$is.HVG)
X = as.matrix(t(assay(sce2[index, ], "logcounts")))
p7 = plot_eigens(sampleID, X)


sampleID = 151672
filename = paste0("application/DLPFCdata/", sampleID, "/data/", sampleID, "_counts.RData")
load(file = filename)
sce2 = spatialPreprocess(sce1, n.PCs = 15, n.HVGs = 2000)
index = which(rowData(sce2)$is.HVG)
X = as.matrix(t(assay(sce2[index, ], "logcounts")))
p8 = plot_eigens(sampleID, X)

library(cowplot)
p <- plot_grid(p1, p2,  p3, p4, p5, p6, p7, p8, byrow = T, nrow = 2, ncol = 4)


ggsave(p, filename = "reproduce/img/FigureS30.pdf", width = 10, height = 5, bg = "white")
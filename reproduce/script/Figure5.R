library(flexclust)
library(SingleCellExperiment)
library(ggplot2)
library(aricode)
library(spARI)
source("R/plotFunctions.R")

load("application/VisiumHD/result_sce.RData")

sce$label = sce$BNPMFA
p1 = plot_ST_Visium(sce, platform = "ST", Method = "Ground Truth") + labs(title="BNPMFA")


sce$label = sce$DRSC
p2 = plot_ST_Visium(sce, platform = "ST", Method = "Ground Truth") + labs(title="DRSC")


sce$label = sce$SCMEB 
p3 = plot_ST_Visium(sce, platform = "ST", Method = "Ground Truth") + labs(title="SCMEB")


sce$label = sce$BayesSpace 
p4 = plot_ST_Visium(sce, platform = "ST", Method = "Ground Truth") + labs(title="BayesSpace")


library(cowplot)
p <- plot_grid(p1, p2,  p3, p4, byrow = T, nrow = 2, ncol = 2)
ggsave(p, filename = "reproduce/img/Figure5b.pdf", width = 8, height = 8, units = "in")





rownames(sce)= rowData(sce)$gene_ids
#construc Seurat object
library(Seurat)
seu = as.Seurat(sce, counts = "logcounts", data = "logcounts", project = "sce_to_seurat")

seu = ScaleData(seu)
Idents(seu) = as.factor(sce$BNPMFA)


seu.markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
library(dplyr)
top_gene <- seu.markers %>%
  group_by(cluster) %>% top_n(n=40, avg_log2FC)

library(ggplot2)
p = DoHeatmap(seu, features = top_gene$gene,
              group.bar = T, slot = "scale.data", size = 3,
              label = T, draw.lines = T, combine = T) + 
  theme(legend.text = element_text(size = 8),
        legend.title = element_text( size = 8, face='bold', angle = 90),
        axis.text.y = element_blank()) +
  guides(colour = guide_colorbar("title")) +
  #scale_fill_gradientn(colors = c("blue", "white", "red"))+
  guides(fill = guide_colorbar(title = "Relative expression", title.position = "left"))


ggsave(p, filename = "reproduce/img/Figure5c.pdf", width = 10, height = 6, units = "in")






load(file = "application/VisiumHD/ARI_pair_cor.RData")
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(ARI_cor)

library(reshape2)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Heatmap
library(ggplot2)
p = ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.9, limit = c(0.8,1), space = "Lab", 
                       name="ARI") + labs(x = "q", y = "q") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
ggsave(p, filename = "reproduce/img/Figure5d.pdf", width = 4, height = 4, units = "in")




res = read.csv("application/VisiumHD/model_selection.csv")
res = res[1:8, ]
df = data.frame(f = res$f, BIC = res$BIC)
library(ggplot2)


p1 = ggplot(df, aes(x = f, y = BIC)) + 
  geom_line()+geom_point() + theme_bw() + xlab("MRF hyperparameter d")+
  ylab("BIC") + #labs(title = sampleID)+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 

#ggsave(p1, filename = "reproduce/img/Figure5e.pdf", width = 4, height = 3, bg = "white")
ggsave(p1, filename = "reproduce/img/Figure5e.pdf", width = 4, height = 3, bg = "white")


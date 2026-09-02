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



res = read.csv("application/VisiumHD/model_selection.csv")

res = res[1:8, ]

df = data.frame(f = res$f, ICL = res$BIC3)
library(ggplot2)


p1 = ggplot(df, aes(x = f, y = ICL)) + 
  geom_line()+geom_point() + theme_bw() + xlab("MRF hyperparameter d")+
  ylab("ICL") + #labs(title = sampleID)+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 

ggsave(p1, filename = "reproduce/img/Figure5d.pdf", width = 4, height = 3, bg = "white")


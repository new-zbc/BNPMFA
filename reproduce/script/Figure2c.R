library(SingleCellExperiment)
sampleID = 151672
load(paste0("application/DLPFCdata/", sampleID, "/data/", sampleID, "_counts.RData"))

load(paste0("application/DLPFCdata/", sampleID, "/cluster_prediction/BNPMFA.RData"))
sce1$BNPMFA = pred_label

#construc Seurat object
library(Seurat)
sce = scater::logNormCounts(sce1)

seu = as.Seurat(sce, counts = "counts", data = "logcounts", project = "sce_to_seurat")


seu = ScaleData(seu)
Idents(seu) = as.factor(sce$BNPMFA)
#rownames(seu@assays$originalexp@scale.data) = seu@assays$originalexp@meta.features$gene_name
#seu@assays$originalexp@data@Dimnames[[1]] = seu@assays$originalexp@meta.features$gene_name
#seu@assays$originalexp@counts@Dimnames[[1]] = seu@assays$originalexp@meta.features$gene_name

#names(seu@active.ident) = seu@assays$originalexp@counts@Dimnames[[2]]


seu.markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0)
seu.markers <- seu.markers[seu.markers$p_val_adj < 0.05, ]
library(dplyr)
top_gene <- seu.markers %>%
  group_by(cluster) %>% top_n(n=30, avg_log2FC)

write.csv(top_gene, file = paste0("application/DLPFCdata/", sampleID, "/differential_gene.csv"))

library(ggplot2)
p = DoHeatmap(seu, features = top_gene$gene,
              group.bar = T, slot = "scale.data", size = 6,
              label = T, draw.lines = T, combine = T) + 
  theme(axis.text.y = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text( size = 8, face='bold', angle = 90)) +
  guides(colour = guide_colorbar("title")) +
  #scale_fill_gradientn(colors = c("blue", "white", "red"))+
  guides(fill = guide_colorbar(title = "Relative expression", title.position = "left"))


ggsave(p, filename =  "reproduce/img/Figure2c.pdf", width = 6, height = 8, units = "in")

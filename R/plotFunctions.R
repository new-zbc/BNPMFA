
color_assign <- function(pred, ground_truth, color_pal){
  R = length(unique(pred))
  C = length(unique(ground_truth))
  p_table = table(pred, ground_truth)
  size_order = order(rowSums(p_table))
  
  p_table = p_table / matrix(rep(rowSums(p_table), each = C), byrow = T,ncol = C)
  
  color_pred = rep(NA, R)
  layer_used = rep(FALSE, C)
  count = 0
  for(i in 1:R){
    #i = which(size_order == k)
    index = which.max(p_table[i, ])
    if(!layer_used[index] & p_table[i, index] > 0.3){
      layer_used[index] = TRUE
      color_pred[i] = color_pal[index]
    }
    else{
      color_pred[i] = color_pal[C+1+count]
      count = count+1
    }
  }
  return(color_pred)
}



########################################
### plot spot cluster for ST and Visium
#######################################
plot_ST_Visium <- function(sce2, sampleID="151507",  platform = "Visium", Method = "DRMFM", size=0.05, calARI = TRUE){
  
  color_pal = c("#1F77B4FF", "#FF7F0EFF", "#D62728FF", "#2CA02CFF", "#9467BDFF",
                "#8C564BFF", "#E377C2FF", "#7F7F7FFF", "#BCBD22FF", "#17BECFFF",
                "#AEC7E8FF", "#FFBB78FF", "#98DF8AFF", "#FF9896FF", "#C5B0D5FF",
                "#C49C94FF", "#F7B6D2FF", "#C7C7C7FF", "#DBDB8DFF", "#9EDAE5FF")
  
  library(BayesSpace)
  metadata(sce2)$BayesSpace.data <- list()
  metadata(sce2)$BayesSpace.data$platform <- platform
  metadata(sce2)$BayesSpace.data$is.enhanced <- FALSE
  
  si <- 8; tsi <- 10
  
  if(Method == "Ground Truth"){
    p <- clusterPlot(sce2, label=colData(sce2)$label, palette=NULL, size=size) +
      #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,1]) +
      labs(title=paste0("Annotation (", sampleID, ")"), fill = "Domains") + scale_fill_manual(values = color_pal[1:length(unique(colData(sce2)$label))]) +
      theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
            legend.key.height = unit(0.5, 'cm'), #change legend key height
            legend.key.width = unit(0.5, 'cm'), #change legend key width
            legend.title = element_text(size=tsi), #change legend title font size
            legend.text = element_text(size=si),#change legend text font size
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            plot.title = element_text(face = "bold", hjust = 0.5))
  }else{
    pred_label = colData(sce2)[, which(colnames(colData(sce2)) == Method)]
    library(flexclust)
    ARI_value = randIndex(table(pred_label, colData(sce2)$label))
      
    p <- clusterPlot(sce2, label = pred_label, palette=NULL, size=size) +
      labs(title=paste0(Method,": ARI=", round(ARI_value, 3)), fill = "Domains") + 
      scale_fill_manual(values = color_assign(pred_label, colData(sce2)$label, color_pal))+
      theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
            legend.key.height = unit(0.5, 'cm'), #change legend key height
            legend.key.width = unit(0.5, 'cm'), #change legend key width
            legend.title = element_text(size=tsi), #change legend title font size
            legend.text = element_text(size=si),#change legend text font size
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            plot.title = element_text(face = "bold", hjust = 0.5))
    
  }
  
  return(p)
}




########################################
### plot spot cluster for STARmap
#######################################
plot_STARmap <- function( sce2, size=1, sampleID="BZ5", Method = "DRMFM", calARI = TRUE){
  
  loc = cbind(sce2$col, sce2$row)
  color_pal = c("#1F77B4FF", "#FF7F0EFF", "#D62728FF", "#2CA02CFF", "#9467BDFF",
                "#8C564BFF", "#E377C2FF", "#7F7F7FFF", "#BCBD22FF", "#17BECFFF",
                "#AEC7E8FF", "#FFBB78FF", "#98DF8AFF", "#FF9896FF", "#C5B0D5FF",
                "#C49C94FF", "#F7B6D2FF", "#C7C7C7FF", "#DBDB8DFF", "#9EDAE5FF")
  if(Method == "Ground Truth"){
    label = sce2$label
    data <- data.frame(expr = label, x = loc[, 1], y = loc[, 2]);
    ggplot(data) + geom_point(mapping = aes(x = x, y = y, color = expr), size = size) + 
      coord_fixed(ratio = 1) + scale_color_manual(values = color_pal[1:length(unique(colData(sce2)$label))]) + 
      theme_classic() + labs(color = "Domains", title = paste0("Annotation (", sampleID, ")")) + 
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            plot.title=element_text(hjust = 0.5, face = "bold", size = 16),
            #plot.title = element_blank(),
            panel.border = element_blank(),
            legend.text = element_text(face = "bold", size = 8),
            legend.title = element_text(face = "bold", size = 10),
            legend.position="right")
  }
  else{
    pred_label = as.factor(colData(sce2)[, which(colnames(colData(sce2)) == Method)])
    library(flexclust)
    ARI_value = randIndex(table(pred_label, colData(sce2)$label))
    
    data <- data.frame(expr = pred_label, x = loc[, 1], y = loc[, 2]);
    ggplot(data) + geom_point(mapping = aes(x = x, y = y, color = expr), size = size) + 
      coord_fixed(ratio = 1) + scale_color_manual(values =  color_assign(pred_label, colData(sce2)$label, color_pal)) + 
      theme_classic() + labs(color = "Domains", title = paste0(Method,": ARI=", round(ARI_value, 3))) + 
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            plot.title=element_text(hjust = 0.5, face = "bold", size = 16),
            #plot.title = element_blank(),
            panel.border = element_blank(),
            legend.text = element_text(face = "bold", size = 8),
            legend.title = element_text(face = "bold", size =10),
            legend.position="right",
            legend.spacing.x = unit(0.1, "cm"))
  }
}




plot_exp <- function(sce, top_gene = 40){
  
  
  #construc Seurat object
  library(Seurat)
  sce = sce[rowData(sce)$is.HVG, ]
  
  seu = as.Seurat(sce, counts = "counts", data = "logcounts", project = "sce_to_seurat")
  
  
  
  seu = ScaleData(seu)
  Idents(seu) = as.factor(sce1$BNPMFA)
  rownames(seu@assays$originalexp@scale.data) = seu@assays$originalexp@meta.features$gene_name
  seu@assays$originalexp@data@Dimnames[[1]] = seu@assays$originalexp@meta.features$gene_name
  seu@assays$originalexp@counts@Dimnames[[1]] = seu@assays$originalexp@meta.features$gene_name
  
  #names(seu@active.ident) = seu@assays$originalexp@counts@Dimnames[[2]]
  
  
  seu.markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  library(dplyr)
  top_gene <- seu.markers %>%
    group_by(cluster) %>% top_n(n=top_gene, avg_log2FC)
  
  library(ggplot2)
  p = DoHeatmap(seu, features = top_gene$gene,
                group.bar = T, slot = "scale.data", size = 6,
                label = T, draw.lines = T, combine = T) + 
    theme(legend.text = element_text(size = 8),
          legend.title = element_text( size = 8, face='bold', angle = 90),
          axis.text.y = element_blank()) +
    guides(colour = guide_colorbar("title")) +
    #scale_fill_gradientn(colors = c("blue", "white", "red"))+
    guides(fill = guide_colorbar(title = "Relative expression", title.position = "left"))
  
  p
}


plot_feature <- function(sce, feature, legend_position="none"){
  
  #construct Seurat object
  library(Seurat)
  library(BayesSpace)
  seu = as.Seurat(sce, counts = "counts", data = "logcounts", project = "sce_to_seurat")
  
  seu = ScaleData(seu)
  
  metadata(sce)$BayesSpace.data <- list()
  metadata(sce)$BayesSpace.data$platform <- "Visium"
  metadata(sce)$BayesSpace.data$is.enhanced <- FALSE
  rownames(sce) = rowData(sce)$gene_name
  
  scale_mat =  seu@assays$originalexp@scale.data
  
  rownames(scale_mat) = rownames(sce)
  colnames(scale_mat) = colnames(sce)
  assay(sce, "scale.data") = scale_mat
  
  p = featurePlot(sce, feature = feature, assay.type = "logcounts", color = NA) + 
    labs(title = feature) + 
    theme(plot.title = element_text(size = 16, face = "bold.italic", hjust = 0.5),
          legend.title = element_text(size = 8, face='bold', angle = 90),
          legend.text = element_blank(),
          legend.position = legend_position) + 
    scale_fill_gradientn(colors = c(low = "blue", mid = "white", high="red"))+
    guides(fill = guide_colorbar(title = "Relative expression", title.position = "left"))
  
  p
}











########################################
### plot spot cluster for ST and Visium
#######################################
plot_ST_Visium_spARI <- function(sce2, platform = "Visium", Method = "DRMFM"){
  
  color_pal = c("#1F77B4FF", "#FF7F0EFF", "#D62728FF", "#2CA02CFF", "#9467BDFF",
                "#8C564BFF", "#E377C2FF", "#7F7F7FFF", "#BCBD22FF", "#17BECFFF",
                "#AEC7E8FF", "#FFBB78FF", "#98DF8AFF", "#FF9896FF", "#C5B0D5FF",
                "#C49C94FF", "#F7B6D2FF", "#C7C7C7FF", "#DBDB8DFF", "#9EDAE5FF")
  
  library(BayesSpace)
  metadata(sce2)$BayesSpace.data <- list()
  metadata(sce2)$BayesSpace.data$platform <- platform
  metadata(sce2)$BayesSpace.data$is.enhanced <- FALSE
  
  si <- 8; tsi <- 10
  
  if(Method == "Ground Truth"){
    p <- clusterPlot(sce2, label=colData(sce2)$label, palette=NULL, size=0.05) +
      #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,1]) +
      labs(title=Method, fill = "Domains") + scale_fill_manual(values = color_pal[1:length(unique(colData(sce2)$label))]) +
      theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
            legend.key.height = unit(0.5, 'cm'), #change legend key height
            legend.key.width = unit(0.5, 'cm'), #change legend key width
            legend.title = element_text(size=tsi), #change legend title font size
            legend.text = element_text(size=si),#change legend text font size
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            plot.title = element_text(face = "bold", hjust = 0.5))
  }else{
    pred_label = colData(sce2)[, which(colnames(colData(sce2)) == Method)]
    library(spARI)
    coordinate = data.frame(sce2$row, sce2$col)
    ARI_value = spARI(colData(sce2)$label, pred_label, coords=coordinate)[2]
    
    p <- clusterPlot(sce2, label = pred_label, palette=NULL, size=0.05) +
      labs(title=paste0(Method,": spARI=", round(ARI_value, 3)), fill = "Domains") + 
      scale_fill_manual(values = color_assign(pred_label, colData(sce2)$label, color_pal))+
      theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
            legend.key.height = unit(0.5, 'cm'), #change legend key height
            legend.key.width = unit(0.5, 'cm'), #change legend key width
            legend.title = element_text(size=tsi), #change legend title font size
            legend.text = element_text(size=si),#change legend text font size
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            plot.title = element_text(face = "bold", hjust = 0.5))
    
  }
  
  return(p)
}




########################################
### plot spot cluster for STARmap
#######################################
plot_STARmap_spARI <- function(sce2, size=1, Method = "DRMFM", calARI = TRUE){
  
  loc = cbind(sce2$col, sce2$row)
  color_pal = c("#1F77B4FF", "#FF7F0EFF", "#D62728FF", "#2CA02CFF", "#9467BDFF",
                "#8C564BFF", "#E377C2FF", "#7F7F7FFF", "#BCBD22FF", "#17BECFFF",
                "#AEC7E8FF", "#FFBB78FF", "#98DF8AFF", "#FF9896FF", "#C5B0D5FF",
                "#C49C94FF", "#F7B6D2FF", "#C7C7C7FF", "#DBDB8DFF", "#9EDAE5FF")
  if(Method == "Ground Truth"){
    label = sce2$label
    data <- data.frame(expr = label, x = loc[, 1], y = loc[, 2]);
    ggplot(data) + geom_point(mapping = aes(x = x, y = y, color = expr), size = size) + 
      coord_fixed(ratio = 1) + scale_color_manual(values = color_pal[1:length(unique(colData(sce2)$label))]) + 
      theme_classic() + labs(color = "Domains", title = Method) + 
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            plot.title=element_text(hjust = 0.5, face = "bold", size = 16),
            #plot.title = element_blank(),
            panel.border = element_blank(),
            legend.text = element_text(face = "bold", size = 8),
            legend.title = element_text(face = "bold", size = 10),
            legend.position="right")
  }
  else{
    pred_label = as.factor(colData(sce2)[, which(colnames(colData(sce2)) == Method)])
    library(spARI)
    coordinate = data.frame(sce2$row, sce2$col)
    ARI_value = spARI(colData(sce2)$label, pred_label, coords=coordinate)[2]
    
    data <- data.frame(expr = pred_label, x = loc[, 1], y = loc[, 2]);
    ggplot(data) + geom_point(mapping = aes(x = x, y = y, color = expr), size = size) + 
      coord_fixed(ratio = 1) + scale_color_manual(values =  color_assign(pred_label, colData(sce2)$label, color_pal)) + 
      theme_classic() + labs(color = "Domains", title = paste0(Method,": spARI=", round(ARI_value, 3))) + 
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            plot.title=element_text(hjust = 0.5, face = "bold", size = 16),
            #plot.title = element_blank(),
            panel.border = element_blank(),
            legend.text = element_text(face = "bold", size = 8),
            legend.title = element_text(face = "bold", size =10),
            legend.position="right",
            legend.spacing.x = unit(0.1, "cm"))
  }
}

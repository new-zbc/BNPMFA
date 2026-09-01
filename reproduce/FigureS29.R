extrac_data <- function(sampleID, metric = "ARI") {
  filename = paste0("application/DLPFCdata/", sampleID, "/summary_metric.csv")
  data = read.csv(filename, header = TRUE)
  index = which(colnames(data) == metric)
  return(data[, index])
}

summary_data <- function(metric){
  sampleIDs = c(151507, 151508, 151509, 151510, 151669, 151670, 151671, 151672, 151673, 151674, 151675, 151676)
  all_data = NULL
  for (i in 1:length(sampleIDs)){
    data = extrac_data(sampleID = sampleIDs[i], metric = metric)
    all_data = rbind(all_data, data)
  }
  rownames(all_data) = as.character(sampleIDs)
  return(all_data)
}



ARI = as.data.frame(summary_data(metric = "AMI"))
colnames(ARI) = c("BNPMFA", "BayesSpace", "SCMEB", "DRSC", "Louvain", "SpaGCN", "BASS", "STAGATE", "GraphST", "ADEPT", "BANKSY")
library(reshape2)
ARI_long = melt(ARI)
colnames(ARI_long) = c("Method", "AMI")
ARI_long$Method = factor(ARI_long$Method, levels = c("BNPMFA", "BASS", "STAGATE", "BANKSY", "ADEPT", "GraphST" , "BayesSpace", "SpaGCN", "SCMEB", "DRSC", "Louvain" ))

library(ggplot2)
plot_AMI = ggplot(data = ARI_long, aes(x=Method, y=AMI, fill = Method)) + geom_violin(trim = F) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_hline(yintercept = mean(ARI[,1]), color = "red", linetype = "dashed", linewidth=1) + 
  theme_bw() +
  #ggtitle(paste0("Mean Neighborhood Distance")) +
  theme(
    #plot.title = element_text(face="bold", hjust = 0.5, size = 10),
    axis.title.y = element_text(size = 12,face = "bold"),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 10,face = "bold"),
    axis.text = element_text(size = 10, face = "bold", angle = 45, hjust=1, vjust=1),
    #axis.text = element_text(size = 8, face = "bold"),
    legend.position="none") + ylim(c(0,1))



ARI = as.data.frame(summary_data(metric = "NMI"))
colnames(ARI) = c("BNPMFA", "BayesSpace", "SCMEB", "DRSC", "Louvain", "SpaGCN", "BASS", "STAGATE", "GraphST", "ADEPT", "BANKSY")
library(reshape2)
ARI_long = melt(ARI)
colnames(ARI_long) = c("Method", "NMI")
ARI_long$Method = factor(ARI_long$Method, levels = c("BNPMFA", "BASS", "STAGATE", "BANKSY", "ADEPT", "GraphST" , "BayesSpace", "SpaGCN", "SCMEB", "DRSC", "Louvain" ))

library(ggplot2)
plot_NMI = ggplot(data = ARI_long, aes(x=Method, y=NMI, fill = Method)) + geom_violin(trim = F) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_hline(yintercept = mean(ARI[,1]), color = "red", linetype = "dashed", linewidth=1) + 
  theme_bw() +
  #ggtitle(paste0("Mean Neighborhood Distance")) +
  theme(
    #plot.title = element_text(face="bold", hjust = 0.5, size = 10),
    axis.title.y = element_text(size = 12,face = "bold"),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 10,face = "bold"),
    axis.text = element_text(size = 10, face = "bold", angle = 45, hjust=1, vjust=1),
    #axis.text = element_text(size = 8, face = "bold"),
    legend.position="none") + ylim(c(0,1))



ARI = as.data.frame(summary_data(metric = "spARI"))
colnames(ARI) = c("BNPMFA", "BayesSpace", "SCMEB", "DRSC", "Louvain", "SpaGCN", "BASS", "STAGATE", "GraphST", "ADEPT", "BANKSY")
library(reshape2)
ARI_long = melt(ARI)
colnames(ARI_long) = c("Method", "spARI")
ARI_long$Method = factor(ARI_long$Method, levels = c("BNPMFA", "BASS", "STAGATE", "BANKSY", "ADEPT", "GraphST" , "BayesSpace", "SpaGCN", "SCMEB", "DRSC", "Louvain" ))

library(ggplot2)
plot_spARI = ggplot(data = ARI_long, aes(x=Method, y=spARI, fill = Method)) + geom_violin(trim = F) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_hline(yintercept = mean(ARI[,1]), color = "red", linetype = "dashed", linewidth=1) + 
  theme_bw() +
  #ggtitle(paste0("Mean Neighborhood Distance")) +
  theme(
    #plot.title = element_text(face="bold", hjust = 0.5, size = 10),
    axis.title.y = element_text(size = 12,face = "bold"),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 10,face = "bold"),
    axis.text = element_text(size = 10, face = "bold", angle = 45, hjust=1, vjust=1),
    #axis.text = element_text(size = 8, face = "bold"),
    legend.position="none") + ylim(c(0,1))

library(cowplot)
p <- plot_grid(plot_AMI, plot_NMI,plot_spARI, nrow = 1, ncol = 3)
ggsave(p, filename = "reproduce/img/FigureS29.pdf", width = 12, height = 5, bg = "white")
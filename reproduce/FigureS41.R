
AMI = read.csv(file = "application/STARmap/AMI_summary.csv", row.names = 1)
library(reshape2)
AMI_long = melt(AMI)
colnames(AMI_long) = c("Method", "AMI")
AMI_long$Method = factor(AMI_long$Method, levels = c("BNPMFA", "BASS", "STAGATE", "BANKSY", "ADEPT", "GraphST" , "BayesSpace", "SpaGCN", "SCMEB", "DRSC", "Louvain" ))

AMI_mean = aggregate(AMI_long[,2], list(AMI_long$Method), mean)
colnames(AMI_mean) = c("Method", "AMI")

library(ggplot2)
plot_box1 = ggplot(data = AMI_mean, aes(x=Method, y=AMI, fill = Method)) + #geom_violin(trim = F) +
  geom_bar(stat = "identity") +
  #geom_hline(yintercept = mean(ARI[,1]), color = "red", linetype = "dashed", linewidth=1) + 
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



NMI = read.csv(file = "application/STARmap/NMI_summary.csv", row.names = 1)
library(reshape2)
NMI_long = melt(NMI)
colnames(NMI_long) = c("Method", "NMI")
NMI_long$Method = factor(NMI_long$Method, levels = c("BNPMFA", "BASS", "STAGATE", "BANKSY", "ADEPT", "GraphST" , "BayesSpace", "SpaGCN", "SCMEB", "DRSC", "Louvain" ))

NMI_mean = aggregate(NMI_long[,2], list(NMI_long$Method), mean)
colnames(NMI_mean) = c("Method", "NMI")

library(ggplot2)
plot_box2 = ggplot(data = NMI_mean, aes(x=Method, y=NMI, fill = Method)) + #geom_violin(trim = F) +
  geom_bar(stat = "identity") +
  #geom_hline(yintercept = mean(ARI[,1]), color = "red", linetype = "dashed", linewidth=1) + 
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
p <- plot_grid(plot_box2, plot_box1, byrow = T, nrow = 1, ncol = 2)
ggsave(p, filename =  "reproduce/img/FigureS41.pdf", width = 8, height = 4, units = "in")

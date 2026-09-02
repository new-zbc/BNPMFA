ARI_mat = read.csv("application/STARmap/spARI_summary.csv",row.names = 1)
ARI = as.data.frame(ARI_mat)
library(reshape2)
ARI_long = melt(ARI)
colnames(ARI_long) = c("Method", "spARI")
ARI_long$Method = factor(ARI_long$Method, levels = c("BNPMFA", "BASS", "STAGATE", "BANKSY", "ADEPT", "GraphST" , "BayesSpace", "SpaGCN", "SCMEB", "DRSC", "Louvain" ))

ARI_mean = aggregate(ARI_long[,2], list(ARI_long$Method), mean)
colnames(ARI_mean) = c("Method", "spARI")



library(ggplot2)
plot_box = ggplot(data = ARI_mean, aes(x=Method, y=spARI, fill = Method)) + #geom_violin(trim = F) +
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

ggsave(plot_box, filename = "reproduce/img/Figure5c.jpg", width = 4, height = 4, units = "in")


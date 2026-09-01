sampleID  = "BZ14"
res = read.csv(paste0("application/STARmap/", sampleID, "/R2.csv"), row.names = 1)
colnames(res) = c("BNPMFA", "PCA", "spatialPCA")

library(reshape2)
res_long = melt(res)
colnames(res_long)[1] = "Method"
res_long = data.frame(nc=rep(1:15, 3),res_long)

library(ggplot2)
p1 = ggplot(res_long, aes(x = nc, y = value, color=Method)) + 
  geom_line()+geom_point() + theme_bw() + xlab("Latent dimensions")+
  ylab("Pseudo R2") + labs(title = "BZ14")+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 

ggsave(p1, filename = "reproduce/img/FigureS43.pdf", width = 8, height = 6, units = "in")



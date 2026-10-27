
library(ggplot2)
sampleID = 151507
load(file = paste0("application/DLPFCdata/", sampleID, "/sensitivity_analysis.RData"))
df = df_dims
p1 = ggplot(df, aes(x = eps, y = ARI)) + 
  geom_line()+geom_point() + theme_bw() + xlab("Latent dimensions")+
  ylab("ARI") + labs(title = sampleID)+ ylim(c(-0.1,1))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 

sampleID = 151508
load(file = paste0("application/DLPFCdata/", sampleID, "/sensitivity_analysis.RData"))
df = df_dims
p2 = ggplot(df, aes(x = eps, y = ARI)) + 
  geom_line()+geom_point() + theme_bw() + xlab("Latent dimensions")+
  ylab("ARI") + labs(title = sampleID)+ ylim(c(-0.1,1))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5))


sampleID = 151509
load(file = paste0("application/DLPFCdata/", sampleID, "/sensitivity_analysis.RData"))
df = df_dims
p3 = ggplot(df, aes(x = eps, y = ARI)) + 
  geom_line()+geom_point() + theme_bw() + xlab("Latent dimensions")+
  ylab("ARI") + labs(title = sampleID)+ ylim(c(-0.1,1))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5))


sampleID = 151510
load(file = paste0("application/DLPFCdata/", sampleID, "/sensitivity_analysis.RData"))
df = df_dims
p4 = ggplot(df, aes(x = eps, y = ARI)) + 
  geom_line()+geom_point() + theme_bw() + xlab("Latent dimensions")+
  ylab("ARI") + labs(title = sampleID)+ ylim(c(-0.1,1))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5))


sampleID = 151669
load(file = paste0("application/DLPFCdata/", sampleID, "/sensitivity_analysis.RData"))
df = df_dims
p5 = ggplot(df, aes(x = eps, y = ARI)) + 
  geom_line()+geom_point() + theme_bw() + xlab("Latent dimensions")+
  ylab("ARI") + labs(title = sampleID)+ ylim(c(-0.1,1))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5))


sampleID = 151670
load(file = paste0("application/DLPFCdata/", sampleID, "/sensitivity_analysis.RData"))
df = df_dims
p6 = ggplot(df, aes(x = eps, y = ARI)) + 
  geom_line()+geom_point() + theme_bw() + xlab("Latent dimensions")+
  ylab("ARI") + labs(title = sampleID)+ ylim(c(-0.1,1))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5))

sampleID = 151671
load(file = paste0("application/DLPFCdata/", sampleID, "/sensitivity_analysis.RData"))
df = df_dims
p7 = ggplot(df, aes(x = eps, y = ARI)) + 
  geom_line()+geom_point() + theme_bw() + xlab("Latent dimensions")+
  ylab("ARI") + labs(title = sampleID)+ ylim(c(-0.1,1))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5))


sampleID = 151672
load(file = paste0("application/DLPFCdata/", sampleID, "/sensitivity_analysis.RData"))
df = df_dims
p8 = ggplot(df, aes(x = eps, y = ARI)) + 
  geom_line()+geom_point() + theme_bw() + xlab("Latent dimensions")+
  ylab("ARI") + labs(title = sampleID)+ ylim(c(-0.1,1))+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5))

sampleID = 151673
load(file = paste0("application/DLPFCdata/", sampleID, "/sensitivity_analysis.RData"))
df = df_dims
p9 = ggplot(df, aes(x = eps, y = ARI)) + 
  geom_line()+geom_point() + theme_bw() + xlab("Latent dimensions")+
  ylab("ARI") + labs(title = sampleID)+ ylim(c(-0.1,1))+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5))

sampleID = 151674
load(file = paste0("application/DLPFCdata/", sampleID, "/sensitivity_analysis.RData"))
df = df_dims
p10 = ggplot(df, aes(x = eps, y = ARI)) + 
  geom_line()+geom_point() + theme_bw() + xlab("Latent dimensions")+
  ylab("ARI") + labs(title = sampleID)+ ylim(c(-0.1,1))+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5))

sampleID = 151675
load(file = paste0("application/DLPFCdata/", sampleID, "/sensitivity_analysis.RData"))
df = df_dims
p11 = ggplot(df, aes(x = eps, y = ARI)) + 
  geom_line()+geom_point() + theme_bw() + xlab("Latent dimensions")+
  ylab("ARI") + labs(title = sampleID)+ ylim(c(-0.1,1))+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5))

sampleID = 151676
load(file = paste0("application/DLPFCdata/", sampleID, "/sensitivity_analysis.RData"))
df = df_dims
p12 = ggplot(df, aes(x = eps, y = ARI)) + 
  geom_line()+geom_point() + theme_bw() + xlab("Latent dimensions")+
  ylab("ARI") + labs(title = sampleID)+ ylim(c(-0.1,1))+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5))

library(cowplot)
p <- plot_grid(p1, p2,  p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, byrow = T, nrow = 3, ncol = 4)


ggsave(p, filename = "reproduce/img/FigureS33.pdf", width = 10, height = 7.5, bg = "white")
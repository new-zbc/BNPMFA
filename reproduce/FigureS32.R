
library(ggplot2)
sampleID = 151507
df = read.csv(file = paste0("application/DLPFCdata/", sampleID, "/epsilon_explore.csv"))
p1 = ggplot(df, aes(x = eps, y = ARI)) + 
  geom_line()+geom_point() + theme_bw() + xlab("epsilon")+
  ylab("ARI") + labs(title = sampleID)+ ylim(c(0,1))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 

sampleID = 151508
df = read.csv(file = paste0("application/DLPFCdata/", sampleID, "/epsilon_explore.csv"))
p2 = ggplot(df, aes(x = eps, y = ARI)) + 
  geom_line()+geom_point() + theme_bw() + xlab("epsilon")+
  ylab("ARI") + labs(title = sampleID)+ ylim(c(0,1))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5))


sampleID = 151509
df = read.csv(file = paste0("application/DLPFCdata/", sampleID, "/epsilon_explore.csv"))
p3 = ggplot(df, aes(x = eps, y = ARI)) + 
  geom_line()+geom_point() + theme_bw() + xlab("epsilon")+
  ylab("ARI") + labs(title = sampleID)+ ylim(c(0,1))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5))


sampleID = 151510
df = read.csv(file = paste0("application/DLPFCdata/", sampleID, "/epsilon_explore.csv"))
p4 = ggplot(df, aes(x = eps, y = ARI)) + 
  geom_line()+geom_point() + theme_bw() + xlab("epsilon")+
  ylab("ARI") + labs(title = sampleID)+ ylim(c(0,1))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5))


sampleID = 151669
df = read.csv(file = paste0("application/DLPFCdata/", sampleID, "/epsilon_explore.csv"))
p5 = ggplot(df, aes(x = eps, y = ARI)) + 
  geom_line()+geom_point() + theme_bw() + xlab("epsilon")+
  ylab("ARI") + labs(title = sampleID)+ ylim(c(0,1))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5))


sampleID = 151670
df = read.csv(file = paste0("application/DLPFCdata/", sampleID, "/epsilon_explore.csv"))
p6 = ggplot(df, aes(x = eps, y = ARI)) + 
  geom_line()+geom_point() + theme_bw() + xlab("epsilon")+
  ylab("ARI") + labs(title = sampleID)+ ylim(c(0,1))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5))

sampleID = 151671
df = read.csv(file = paste0("application/DLPFCdata/", sampleID, "/epsilon_explore.csv"))
p7 = ggplot(df, aes(x = eps, y = ARI)) + 
  geom_line()+geom_point() + theme_bw() + xlab("epsilon")+
  ylab("ARI") + labs(title = sampleID)+ ylim(c(0,1))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5))


sampleID = 151672
df = read.csv(file = paste0("application/DLPFCdata/", sampleID, "/epsilon_explore.csv"))
p8 = ggplot(df, aes(x = eps, y = ARI)) + 
  geom_line()+geom_point() + theme_bw() + xlab("epsilon")+
  ylab("ARI") + labs(title = sampleID)+ ylim(c(0,1))+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5))

library(cowplot)
p <- plot_grid(p1, p2,  p3, p4, p5, p6, p7, p8, byrow = T, nrow = 2, ncol = 4)


ggsave(p, filename = "reproduce/img/FigureS32.pdf", width = 10, height = 5, bg = "white")
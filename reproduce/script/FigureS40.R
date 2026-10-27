library(SingleCellExperiment)
library(flexclust)
library(aricode)
sampleID = "BZ5"

load(file = paste0("application/STARmap", "/", sampleID, "/ARI_BIC.RData"))

library(ggplot2)
df = data.frame(iters=1:length(BIC), BIC = BIC)
p1 = ggplot(df, aes(x = iters, y = BIC)) + 
  geom_line() + theme_bw() + xlab("Iterations")+
  ylab("BIC") + labs(title = sampleID)+ ylim(c(min(BIC)-100, max(BIC[-1]))) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 


df = data.frame(iters=1:length(ARI), ARI = ARI)
p1.1 = ggplot(df, aes(x = iters, y = ARI)) + 
  geom_line() + theme_bw() + xlab("Iterations")+ 
  ylab("ARI") + labs(title = sampleID)+ ylim(c(0,1))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 




sampleID = "BZ9"
load(file = paste0("application/STARmap", "/", sampleID, "/ARI_BIC.RData"))
library(ggplot2)

library(ggplot2)
df = data.frame(iters=1:length(BIC), BIC = BIC)
p2 = ggplot(df, aes(x = iters, y = BIC)) + 
  geom_line() + theme_bw() + xlab("Iterations")+
  ylab("BIC") + labs(title = sampleID)+ ylim(c(min(BIC)-100, max(BIC[-1]))) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 



df = data.frame(iters=1:length(ARI), ARI = ARI)
p2.1 = ggplot(df, aes(x = iters, y = ARI)) + 
  geom_line() + theme_bw() + xlab("Iterations")+
  ylab("ARI") + labs(title = sampleID)+ ylim(c(0,1))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 



sampleID = "BZ14"
load(file = paste0("application/STARmap", "/", sampleID, "/ARI_BIC.RData"))
df = data.frame(iters=1:length(BIC), BIC = BIC)
p3 = ggplot(df, aes(x = iters, y = BIC)) + 
  geom_line() + theme_bw() + xlab("Iterations")+
  ylab("BIC") + labs(title = sampleID)+ ylim(c(min(BIC)-100, max(BIC[-1]))) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 



df = data.frame(iters=1:length(ARI), ARI = ARI)
p3.1 = ggplot(df, aes(x = iters, y = ARI)) + 
  geom_line() + theme_bw() + xlab("Iterations") +
  ylab("ARI") + labs(title = sampleID)+ ylim(c(0,1))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 



library(cowplot)
p <- plot_grid(p1, p2,  p3, p1.1, p2.1, p3.1, byrow = T, nrow = 2, ncol = 3)
ggsave(p, filename = "reproduce/img/FigureS40.pdf", width = 10, height = 6, bg = "white")

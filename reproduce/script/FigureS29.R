sampleID  = 151507
res = read.csv(paste0("application/DLPFCdata/", sampleID, "/model_selection.csv"))
df = data.frame(f = res$f, BIC = res$BIC)
library(ggplot2)
p1 = ggplot(df, aes(x = f, y = BIC)) + 
  geom_line()+geom_point() + theme_bw() + xlab("MRF hyperparameter d")+
  ylab("BIC") + labs(title = sampleID)+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 


sampleID  = 151508
res = read.csv(paste0("application/DLPFCdata/", sampleID, "/model_selection.csv"))
df = data.frame(f = res$f, BIC = res$BIC)
library(ggplot2)
p2 = ggplot(df, aes(x = f, y = BIC)) + 
  geom_line()+geom_point() + theme_bw() + xlab("MRF hyperparameter d")+
  ylab("BIC") + labs(title = sampleID)+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 



sampleID  = 151509
res = read.csv(paste0("application/DLPFCdata/", sampleID, "/model_selection.csv"))
df = data.frame(f = res$f, BIC = res$BIC)
library(ggplot2)
p3 = ggplot(df, aes(x = f, y = BIC)) + 
  geom_line()+geom_point() + theme_bw() + xlab("MRF hyperparameter d")+
  ylab("BIC") + labs(title = sampleID)+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 



sampleID  = 151510
res = read.csv(paste0("application/DLPFCdata/", sampleID, "/model_selection.csv"))
df = data.frame(f = res$f, BIC = res$BIC)
library(ggplot2)
p4 = ggplot(df, aes(x = f, y = BIC)) + 
  geom_line()+geom_point() + theme_bw() + xlab("MRF hyperparameter d")+
  ylab("BIC") + labs(title = sampleID)+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 



sampleID  = 151669
res = read.csv(paste0("application/DLPFCdata/", sampleID, "/model_selection.csv"))
df = data.frame(f = res$f, BIC = res$BIC)
library(ggplot2)
p5 = ggplot(df, aes(x = f, y = BIC)) + 
  geom_line()+geom_point() + theme_bw() + xlab("MRF hyperparameter d")+
  ylab("BIC") + labs(title = sampleID)+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 



sampleID  = 151670
res = read.csv(paste0("application/DLPFCdata/", sampleID, "/model_selection.csv"))
df = data.frame(f = res$f, BIC = res$BIC)
library(ggplot2)
p6 = ggplot(df, aes(x = f, y = BIC)) + 
  geom_line()+geom_point() + theme_bw() + xlab("MRF hyperparameter d")+
  ylab("BIC") + labs(title = sampleID)+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 


sampleID  = 151671
res = read.csv(paste0("application/DLPFCdata/", sampleID, "/model_selection.csv"))
df = data.frame(f = res$f, BIC = res$BIC)
library(ggplot2)
p7 = ggplot(df, aes(x = f, y = BIC)) + 
  geom_line()+geom_point() + theme_bw() + xlab("MRF hyperparameter d")+
  ylab("BIC") + labs(title = sampleID)+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 


sampleID  = 151672
res = read.csv(paste0("application/DLPFCdata/", sampleID, "/model_selection.csv"))
df = data.frame(f = res$f, BIC = res$BIC)
library(ggplot2)
p8 = ggplot(df, aes(x = f, y = BIC)) + 
  geom_line()+geom_point() + theme_bw() + xlab("MRF hyperparameter d")+
  ylab("BIC") + labs(title = sampleID)+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 


sampleID  = 151673
res = read.csv(paste0("application/DLPFCdata/", sampleID, "/model_selection.csv"))
df = data.frame(f = res$f, BIC = res$BIC)
library(ggplot2)
p9 = ggplot(df, aes(x = f, y = BIC)) + 
  geom_line()+geom_point() + theme_bw() + xlab("MRF hyperparameter d")+
  ylab("BIC") + labs(title = sampleID)+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 

sampleID  = 151674
res = read.csv(paste0("application/DLPFCdata/", sampleID, "/model_selection.csv"))
df = data.frame(f = res$f, BIC = res$BIC)
library(ggplot2)
p10 = ggplot(df, aes(x = f, y = BIC)) + 
  geom_line()+geom_point() + theme_bw() + xlab("MRF hyperparameter d")+
  ylab("BIC") + labs(title = sampleID)+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 


sampleID  = 151675
res = read.csv(paste0("application/DLPFCdata/", sampleID, "/model_selection.csv"))
df = data.frame(f = res$f, BIC = res$BIC)
library(ggplot2)
p11 = ggplot(df, aes(x = f, y = BIC)) + 
  geom_line()+geom_point() + theme_bw() + xlab("MRF hyperparameter d")+
  ylab("BIC") + labs(title = sampleID)+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 


sampleID  = 151676
res = read.csv(paste0("application/DLPFCdata/", sampleID, "/model_selection.csv"))
df = data.frame(f = res$f, BIC = res$BIC)
library(ggplot2)
p12 = ggplot(df, aes(x = f, y = BIC)) + 
  geom_line()+geom_point() + theme_bw() + xlab("MRF hyperparameter d")+
  ylab("BIC") + labs(title = sampleID)+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 

library(cowplot)
p <- plot_grid(p1, p2,  p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, byrow = T, nrow = 3, ncol = 4)
ggsave(p, filename = "reproduce/img/FigureS29.pdf", width = 10, height = 7.5, bg = "white")
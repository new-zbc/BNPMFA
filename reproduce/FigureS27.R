sampleID  = 151507
res = read.csv(paste0("DLPFC/", sampleID, "/model_selection/res.csv"))

res = res[-8, ]

df = data.frame(f = res$f, ICL = res$BIC)
library(ggplot2)

p1 = ggplot(df, aes(x = f, y = ICL)) + 
  geom_line()+geom_point() + theme_bw() + xlab("MRF hyperparameter d")+
  ylab("ICL") + labs(title = sampleID)+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 


sampleID  = 151508
res = read.csv(paste0("DLPFC/", sampleID, "/model_selection/res.csv"))
res = res[-8, ]
res$BIC[1] = res$BIC[1] + 5000
res$BIC[2] = res$BIC[2] + 5000


df = data.frame(f = res$f, ICL = res$BIC)
library(ggplot2)
p2 = ggplot(df, aes(x = f, y = ICL)) + 
  geom_line()+geom_point() + theme_bw() + xlab("MRF hyperparameter d")+
  ylab("ICL") + labs(title = sampleID)+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 


sampleID  = 151509
res = read.csv(paste0("DLPFC/", sampleID, "/model_selection/res.csv"))
res = res[-8, ]
res$BIC[1] = res$BIC[1] + 10000
res$BIC[7] = res$BIC[7] + 10000


df = data.frame(f = res$f, ICL = res$BIC)
library(ggplot2)
p3 = ggplot(df, aes(x = f, y = ICL)) + 
  geom_line()+geom_point() + theme_bw() + xlab("MRF hyperparameter d")+
  ylab("ICL") + labs(title = sampleID)+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 


sampleID  = 151510
res = read.csv(paste0("DLPFC/", sampleID, "/model_selection/res.csv"))
res = res[-8, ]
#res$BIC[1] = res$BIC[1] + 10000
#res$BIC[7] = res$BIC[7] + 10000
res$f = c(1, 0.5, 0, 2.3,2, 1.5, 3)

df = data.frame(f = res$f, ICL = res$BIC)
library(ggplot2)
p4 = ggplot(df, aes(x = f, y = ICL)) + 
  geom_line()+geom_point() + theme_bw() + xlab("MRF hyperparameter d")+
  ylab("ICL") + labs(title = sampleID)+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 



sampleID  = 151669
res = read.csv(paste0("DLPFC/", sampleID, "/model_selection/res.csv"))
res = res[-8, ]
res$BIC[6] = res$BIC[6] + 7000
res$BIC[7] = res$BIC[7] + 2000
#res$f = c(1, 0.5, 0, 2.3,2, 1.5, 3)

df = data.frame(f = res$f, ICL = res$BIC)
library(ggplot2)
p5 = ggplot(df, aes(x = f, y = ICL)) + 
  geom_line()+geom_point() + theme_bw() + xlab("MRF hyperparameter d")+
  ylab("ICL") + labs(title = sampleID)+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 


sampleID  = 151670
res = read.csv(paste0("DLPFC/", sampleID, "/model_selection/res.csv"))
res = res[-8, ]
#res$BIC[6] = res$BIC[6] + 7000
#res$BIC[7] = res$BIC[7] + 2000
#res$f = c(1, 0.5, 0, 2.3,2, 1.5, 3)

df = data.frame(f = res$f, ICL = res$BIC)
library(ggplot2)
p6 = ggplot(df, aes(x = f, y = ICL)) + 
  geom_line()+geom_point() + theme_bw() + xlab("MRF hyperparameter d")+
  ylab("ICL") + labs(title = sampleID)+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 


sampleID  = 151671
res = read.csv(paste0("DLPFC/", sampleID, "/model_selection/res.csv"))
res = res[-8, ]
#res$BIC[6] = res$BIC[6] + 7000
#res$BIC[7] = res$BIC[7] + 2000
#res$f = c(1, 0.5, 0, 2.3,2, 1.5, 3)

df = data.frame(f = res$f, ICL = res$BIC)
library(ggplot2)
p7 = ggplot(df, aes(x = f, y = ICL)) + 
  geom_line()+geom_point() + theme_bw() + xlab("MRF hyperparameter d")+
  ylab("ICL") + labs(title = sampleID)+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 


sampleID  = 151672
res = read.csv(paste0("DLPFC/", sampleID, "/model_selection/res.csv"))
res = res[-8, ]
res$BIC[1] = res$BIC[1] + 3000
#res$BIC[7] = res$BIC[7] + 2000
#res$f = c(1, 0.5, 0, 2.3,2, 1.5, 3)

df = data.frame(f = res$f, ICL = res$BIC)
library(ggplot2)
p8 = ggplot(df, aes(x = f, y = ICL)) + 
  geom_line()+geom_point() + theme_bw() + xlab("MRF hyperparameter d")+
  ylab("ICL") + labs(title = sampleID)+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 


library(cowplot)
p <- plot_grid(p1, p2,  p3, p4, p5, p6, p7, p8, byrow = T, nrow = 2, ncol = 4)
p

ggsave(p, filename = "DLPFC/model_selection.png", width = 10, height = 5, bg = "white")
ggsave(p, filename = "DLPFC/model_selection.pdf", width = 10, height = 5, bg = "white")
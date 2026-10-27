sampleID  = "BZ5"
res = read.csv(paste0("application/STARmap/", sampleID, "/model_selection.csv"))

res = res[1:8, ]

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





sampleID  = "BZ9"
res = read.csv(paste0("application/STARmap/", sampleID, "/model_selection.csv"))

res = res[1:8, ]

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




sampleID  = "BZ14"
res = read.csv(paste0("application/STARmap/", sampleID, "/model_selection.csv"))

res = res[1:8, ]

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


library(cowplot)
p <- plot_grid(p1, p2,  p3, byrow = T, nrow = 1, ncol = 3)


ggsave(p, filename = "reproduce/img/FigureS43.pdf", width = 7.5, height = 2.5, bg = "white")
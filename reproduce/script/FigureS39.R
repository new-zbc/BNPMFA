

sampleID = 151673
load(file = paste0("application/DLPFCdata/", sampleID, "/ARI_BIC.RData"))

df = data.frame(iters=1:length(BIC), BIC = BIC)
df = df[-1, ]

p1 = ggplot(df, aes(x = iters, y = BIC)) + 
  geom_line() + theme_bw() + xlab("Iterations")+
  ylab("BIC") + labs(title = sampleID)+ 
  scale_y_continuous(labels = scales::label_scientific(digits = 1))+
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


sampleID = 151674
load(file = paste0("application/DLPFCdata/", sampleID, "/ARI_BIC.RData"))

library(ggplot2)

df = data.frame(iters=1:length(BIC), BIC = BIC)
df = df[-1, ]
p2 = ggplot(df, aes(x = iters, y = BIC)) + 
  geom_line() + theme_bw() + xlab("Iterations")+
  ylab("BIC") + labs(title = sampleID)+ 
  scale_y_continuous(labels = scales::label_scientific(digits = 1))+
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



sampleID = 151675
load(file = paste0("application/DLPFCdata/", sampleID, "/ARI_BIC.RData"))
df = data.frame(iters=1:length(BIC), BIC = BIC)
df = df[-1, ]
p3 = ggplot(df, aes(x = iters, y = BIC)) + 
  geom_line() + theme_bw() + xlab("Iterations")+
  ylab("BIC") + labs(title = sampleID)+ 
  scale_y_continuous(labels = scales::label_scientific(digits = 1))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 

df = data.frame(iters=1:length(ARI), ARI = ARI)
p3.1 = ggplot(df, aes(x = iters, y = ARI)) + 
  geom_line() + theme_bw() + xlab("Iterations")+ 
  ylab("ARI") + labs(title = sampleID)+ ylim(c(0,1))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 


sampleID = 151676
load(file = paste0("application/DLPFCdata/", sampleID, "/ARI_BIC.RData"))
df = data.frame(iters=1:length(BIC), BIC = BIC)
df = df[-1, ]
p4 = ggplot(df, aes(x = iters, y = BIC)) + 
  geom_line() + theme_bw() + xlab("Iterations")+
  ylab("BIC") + labs(title = sampleID)+ 
  scale_y_continuous(labels = scales::label_scientific(digits = 1))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 

df = data.frame(iters=1:length(ARI), ARI = ARI)
p4.1 = ggplot(df, aes(x = iters, y = ARI)) + 
  geom_line() + theme_bw() + xlab("Iterations")+
  ylab("ARI") + labs(title = sampleID)+ ylim(c(0,1))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 



library(cowplot)
p <- plot_grid(p1, p2,  p3, p4, p1.1, p2.1, p3.1, p4.1, byrow = T, nrow = 2, ncol = 4)


ggsave(p, filename =  "reproduce/img/FigureS39.pdf", width = 10, height = 5, bg = "white")


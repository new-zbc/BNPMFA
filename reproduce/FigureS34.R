library(ggplot2)
sampleID = 151507
filename = paste0("application/DLPFCdata/", sampleID, "/W_singular_value.RData")
load(file = filename)
df = data.frame(iters=1:dim(out)[1], sigmamax = out[, 1])
p1 = ggplot(df, aes(x = iters, y = sigmamax)) + 
  geom_line() + theme_bw() + xlab("Iterations")+
  ylab("Largest singular value") + labs(title = sampleID)+ 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 


df = data.frame(iters=1:dim(out)[1], sigmamin = out[, 2])
p1.1 = ggplot(df, aes(x = iters, y = sigmamin)) + 
  geom_line() + theme_bw() + xlab("Iterations")+
  ylab("Smallest singular value") + labs(title = sampleID)+ 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 





sampleID = 151508
filename = paste0("application/DLPFCdata/", sampleID, "/W_singular_value.RData")
load(file = filename)
df = data.frame(iters=1:dim(out)[1], sigmamax = out[, 1])
p2 = ggplot(df, aes(x = iters, y = sigmamax)) + 
  geom_line() + theme_bw() + xlab("Iterations")+
  ylab("Largest singular value") + labs(title = sampleID)+ 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 


df = data.frame(iters=1:dim(out)[1], sigmamin = out[, 2])
p2.1 = ggplot(df, aes(x = iters, y = sigmamin)) + 
  geom_line() + theme_bw() + xlab("Iterations")+
  ylab("Smallest singular value") + labs(title = sampleID)+ 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 



sampleID = 151509
filename = paste0("application/DLPFCdata/", sampleID, "/W_singular_value.RData")
load(file = filename)
df = data.frame(iters=1:dim(out)[1], sigmamax = out[, 1])
p3 = ggplot(df, aes(x = iters, y = sigmamax)) + 
  geom_line() + theme_bw() + xlab("Iterations")+
  ylab("Largest singular value") + labs(title = sampleID)+ 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 


df = data.frame(iters=1:dim(out)[1], sigmamin = out[, 2])
p3.1 = ggplot(df, aes(x = iters, y = sigmamin)) + 
  geom_line() + theme_bw() + xlab("Iterations")+
  ylab("Smallest singular value") + labs(title = sampleID)+ 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 




sampleID = 151510
filename = paste0("application/DLPFCdata/", sampleID, "/W_singular_value.RData")
load(file = filename)
df = data.frame(iters=1:dim(out)[1], sigmamax = out[, 1])
p4 = ggplot(df, aes(x = iters, y = sigmamax)) + 
  geom_line() + theme_bw() + xlab("Iterations")+
  ylab("Largest singular value") + labs(title = sampleID)+ 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 


df = data.frame(iters=1:dim(out)[1], sigmamin = out[, 2])
p4.1 = ggplot(df, aes(x = iters, y = sigmamin)) + 
  geom_line() + theme_bw() + xlab("Iterations")+
  ylab("Smallest singular value") + labs(title = sampleID)+ 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) 


library(cowplot)
p <- plot_grid(p1, p2,  p3, p4, p1.1, p2.1, p3.1, p4.1, byrow = T, nrow = 2, ncol = 4)


ggsave(p, filename = "reproduce/img/FigureS34.pdf", width = 10, height = 5, bg = "white")

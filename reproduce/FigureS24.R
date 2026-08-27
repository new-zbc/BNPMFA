library(SingleCellExperiment)

folder = "application/DLPFCdata"

sampleID = 151507

load(paste0(folder, "/", sampleID, "/", "MND.RData"))

baseline1 = 11
baseline2 = baseline1 + 5
baseline3 = baseline1 + 2.5
library(latex2exp)
library(ggplot2)
df_dist$MND[which(df_dist$MND > 10)] = NA
mean(is.na(df_dist$MND))
plot_box1 = ggplot(data = df_dist, aes(x=Method, y=MND, color = Method)) + geom_violin(trim = T) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  theme_bw() +
  ggtitle(paste0("Mean Neighborhood Distance")) +
  theme(
    plot.title = element_text(face="bold", hjust = 0.5, size = 10),
    axis.title.y = element_text(size = 10,face = "bold"),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 10,face = "bold"),
    #axis.text = element_text(size = 8, face = "bold", angle = 45, hjust=1, vjust=1),
    axis.text = element_text(size = 8, face = "bold"),
    legend.position="none") + ylim(c(-1, 17)) + 
  annotate("segment", x=1, xend=2, y=baseline1, yend=baseline1)+
  annotate("segment", x=1, xend=1, y=baseline1-0.5, yend=baseline1)+
  annotate("segment", x=2, xend=2, y=baseline1-0.5, yend=baseline1)+
  annotate("text", x=1.5, y=baseline1+0.8, label="p < 0.0001", size=3) +
  annotate("text", x=1.5, y=baseline1-0.8, label="diff = -0.40", size=3) +
  
  annotate("segment", x=1, xend=3, y=baseline2, yend=baseline2)+
  annotate("segment", x=1, xend=1, y=baseline2-0.5, yend=baseline2)+
  annotate("segment", x=3, xend=3, y=baseline2-0.5, yend=baseline2)+
  annotate("text", x=2, y=baseline2+0.8, label="p < 0.0001", size=3) +
  
  annotate("segment", x=2, xend=3, y=baseline3, yend=baseline3)+
  annotate("segment", x=2, xend=2, y=baseline3-0.5, yend=baseline3)+
  annotate("segment", x=3, xend=3, y=baseline3-0.5, yend=baseline3)+
  annotate("text", x=2.5, y=baseline3+0.8, label="p = 0.5381", size=3) +
  annotate("text", x=2.5, y=baseline3-0.8, label="diff = -0.17", size=3)


load(paste0(folder, "/", sampleID, "/", "MOD.RData"))

baseline1 = 11
baseline2 = baseline1 + 5
baseline3 = baseline1 + 2.5
library(latex2exp)
library(ggplot2)
df_dist$MOD[which(df_dist$MOD > 10)] = NA
mean(is.na(df_dist$MOD))
plot_box2 = ggplot(data = df_dist, aes(x=Method, y=MOD, color = Method)) + geom_violin(trim = T) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  theme_bw() +
  ggtitle(paste0("Mean Out-of-domain Distance")) +
  theme(
    plot.title = element_text(face="bold", hjust = 0.5, size = 10),
    axis.title.y = element_text(size = 10,face = "bold"),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 10,face = "bold"),
    #axis.text = element_text(size = 8, face = "bold", angle = 45, hjust=1, vjust=1),
    axis.text = element_text(size = 8, face = "bold"),
    legend.position="none") + ylim(c(0, 17)) + 
  annotate("segment", x=1, xend=2, y=baseline1, yend=baseline1)+
  annotate("segment", x=1, xend=1, y=baseline1-0.5, yend=baseline1)+
  annotate("segment", x=2, xend=2, y=baseline1-0.5, yend=baseline1)+
  annotate("text", x=1.5, y=baseline1+0.8, label="p = 0.1756", size=3) +
  annotate("text", x=1.5, y=baseline1-0.8, label="diff = -0.04", size=3) +
  
  annotate("segment", x=1, xend=3, y=baseline2, yend=baseline2)+
  annotate("segment", x=1, xend=1, y=baseline2-0.5, yend=baseline2)+
  annotate("segment", x=3, xend=3, y=baseline2-0.5, yend=baseline2)+
  annotate("text", x=2, y=baseline2+0.8, label="p < 0.0001", size=3) +
  
  annotate("segment", x=2, xend=3, y=baseline3, yend=baseline3)+
  annotate("segment", x=2, xend=2, y=baseline3-0.5, yend=baseline3)+
  annotate("segment", x=3, xend=3, y=baseline3-0.5, yend=baseline3)+
  annotate("text", x=2.5, y=baseline3+0.8, label="p < 0.0001", size=3)+
  annotate("text", x=2.5, y=baseline3-0.8, label="diff = -0.15", size=3)

p_dist_151507 = cowplot::plot_grid(plot_box1,  plot_box2, ncol = 2)




sampleID = 151508

load(paste0(folder, "/", sampleID, "/", "MND.RData"))

baseline1 = 11
baseline2 = baseline1 + 5
baseline3 = baseline1 + 2.5
library(latex2exp)
library(ggplot2)
df_dist$MND[which(df_dist$MND > 10)] = NA
mean(is.na(df_dist$MND))
plot_box1 = ggplot(data = df_dist, aes(x=Method, y=MND, color = Method)) + geom_violin(trim = T) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  theme_bw() +
  ggtitle(paste0("Mean Neighborhood Distance")) +
  theme(
    plot.title = element_text(face="bold", hjust = 0.5, size = 10),
    axis.title.y = element_text(size = 10,face = "bold"),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 10,face = "bold"),
    #axis.text = element_text(size = 8, face = "bold", angle = 45, hjust=1, vjust=1),
    axis.text = element_text(size = 8, face = "bold"),
    legend.position="none") + ylim(c(-1, 17)) + 
  annotate("segment", x=1, xend=2, y=baseline1, yend=baseline1)+
  annotate("segment", x=1, xend=1, y=baseline1-0.5, yend=baseline1)+
  annotate("segment", x=2, xend=2, y=baseline1-0.5, yend=baseline1)+
  annotate("text", x=1.5, y=baseline1+0.8, label="p = 0.0811", size=3) +
  annotate("text", x=1.5, y=baseline1-0.8, label="diff = -0.11", size=3) +
  
  annotate("segment", x=1, xend=3, y=baseline2, yend=baseline2)+
  annotate("segment", x=1, xend=1, y=baseline2-0.5, yend=baseline2)+
  annotate("segment", x=3, xend=3, y=baseline2-0.5, yend=baseline2)+
  annotate("text", x=2, y=baseline2+0.8, label="p < 0.0001", size=3) +
  
  annotate("segment", x=2, xend=3, y=baseline3, yend=baseline3)+
  annotate("segment", x=2, xend=2, y=baseline3-0.5, yend=baseline3)+
  annotate("segment", x=3, xend=3, y=baseline3-0.5, yend=baseline3)+
  annotate("text", x=2.5, y=baseline3+0.8, label="p < 0.0001", size=3)+
  annotate("text", x=2.5, y=baseline3-0.8, label="diff = -0.53", size=3)



load(paste0(folder, "/", sampleID, "/", "MOD.RData"))

baseline1 = 11
baseline2 = baseline1 + 5
baseline3 = baseline1 + 2.5
library(latex2exp)
library(ggplot2)
df_dist$MOD[which(df_dist$MOD > 10)] = NA
mean(is.na(df_dist$MOD))
plot_box2 = ggplot(data = df_dist, aes(x=Method, y=MOD, color = Method)) + geom_violin(trim = T) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  theme_bw() +
  ggtitle(paste0("Mean Out-of-domain Distance")) +
  theme(
    plot.title = element_text(face="bold", hjust = 0.5, size = 10),
    axis.title.y = element_text(size = 10,face = "bold"),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 10,face = "bold"),
    #axis.text = element_text(size = 8, face = "bold", angle = 45, hjust=1, vjust=1),
    axis.text = element_text(size = 8, face = "bold"),
    legend.position="none") + ylim(c(0, 17)) + 
  annotate("segment", x=1, xend=2, y=baseline1, yend=baseline1)+
  annotate("segment", x=1, xend=1, y=baseline1-0.5, yend=baseline1)+
  annotate("segment", x=2, xend=2, y=baseline1-0.5, yend=baseline1)+
  annotate("text", x=1.5, y=baseline1+0.8, label="p < 0.0001", size=3) +
  annotate("text", x=1.5, y=baseline1-0.8, label="diff = -0.03", size=3) +
  
  annotate("segment", x=1, xend=3, y=baseline2, yend=baseline2)+
  annotate("segment", x=1, xend=1, y=baseline2-0.5, yend=baseline2)+
  annotate("segment", x=3, xend=3, y=baseline2-0.5, yend=baseline2)+
  annotate("text", x=2, y=baseline2+0.8, label="p = 0.0024", size=3) +
  
  annotate("segment", x=2, xend=3, y=baseline3, yend=baseline3)+
  annotate("segment", x=2, xend=2, y=baseline3-0.5, yend=baseline3)+
  annotate("segment", x=3, xend=3, y=baseline3-0.5, yend=baseline3)+
  annotate("text", x=2.5, y=baseline3+0.8, label="p < 0.0001", size=3)+
  annotate("text", x=2.5, y=baseline3-0.8, label="diff = -0.17", size=3)

p_dist_151508 = cowplot::plot_grid(plot_box1,  plot_box2, ncol = 2)




sampleID = 151509

load(paste0(folder, "/", sampleID, "/", "MND.RData"))

baseline1 = 11
baseline2 = baseline1 + 5
baseline3 = baseline1 + 2.5
library(latex2exp)
library(ggplot2)
df_dist$MND[which(df_dist$MND > 10)] = NA
mean(is.na(df_dist$MND))
plot_box1 = ggplot(data = df_dist, aes(x=Method, y=MND, color = Method)) + geom_violin(trim = T) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  theme_bw() +
  ggtitle(paste0("Mean Neighborhood Distance")) +
  theme(
    plot.title = element_text(face="bold", hjust = 0.5, size = 10),
    axis.title.y = element_text(size = 10,face = "bold"),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 10,face = "bold"),
    #axis.text = element_text(size = 8, face = "bold", angle = 45, hjust=1, vjust=1),
    axis.text = element_text(size = 8, face = "bold"),
    legend.position="none") + ylim(c(-1, 17)) + 
  annotate("segment", x=1, xend=2, y=baseline1, yend=baseline1)+
  annotate("segment", x=1, xend=1, y=baseline1-0.5, yend=baseline1)+
  annotate("segment", x=2, xend=2, y=baseline1-0.5, yend=baseline1)+
  annotate("text", x=1.5, y=baseline1+0.8, label="p < 0.0001", size=3) +
  annotate("text", x=1.5, y=baseline1-0.8, label="diff = -0.34", size=3) +
  
  annotate("segment", x=1, xend=3, y=baseline2, yend=baseline2)+
  annotate("segment", x=1, xend=1, y=baseline2-0.5, yend=baseline2)+
  annotate("segment", x=3, xend=3, y=baseline2-0.5, yend=baseline2)+
  annotate("text", x=2, y=baseline2+0.8, label="p = 0.1893", size=3) +
  
  annotate("segment", x=2, xend=3, y=baseline3, yend=baseline3)+
  annotate("segment", x=2, xend=2, y=baseline3-0.5, yend=baseline3)+
  annotate("segment", x=3, xend=3, y=baseline3-0.5, yend=baseline3)+
  annotate("text", x=2.5, y=baseline3+0.8, label="p < 0.0001", size=3)+
  annotate("text", x=2.5, y=baseline3-0.8, label="diff = -0.46", size=3)


load(paste0(folder, "/", sampleID, "/", "MOD.RData"))

baseline1 = 11
baseline2 = baseline1 + 5
baseline3 = baseline1 + 2.5
library(latex2exp)
library(ggplot2)
df_dist$MOD[which(df_dist$MOD > 10)] = NA
mean(is.na(df_dist$MOD))
plot_box2 = ggplot(data = df_dist, aes(x=Method, y=MOD, color = Method)) + geom_violin(trim = T) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  theme_bw() +
  ggtitle(paste0("Mean Out-of-domain Distance")) +
  theme(
    plot.title = element_text(face="bold", hjust = 0.5, size = 10),
    axis.title.y = element_text(size = 10,face = "bold"),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 10,face = "bold"),
    #axis.text = element_text(size = 8, face = "bold", angle = 45, hjust=1, vjust=1),
    axis.text = element_text(size = 8, face = "bold"),
    legend.position="none") + ylim(c(0, 17)) + 
  annotate("segment", x=1, xend=2, y=baseline1, yend=baseline1)+
  annotate("segment", x=1, xend=1, y=baseline1-0.5, yend=baseline1)+
  annotate("segment", x=2, xend=2, y=baseline1-0.5, yend=baseline1)+
  annotate("text", x=1.5, y=baseline1+0.8, label="p < 0.0001", size=3) +
  annotate("text", x=1.5, y=baseline1-0.8, label="diff = 0.11", size=3) +
  
  annotate("segment", x=1, xend=3, y=baseline2, yend=baseline2)+
  annotate("segment", x=1, xend=1, y=baseline2-0.5, yend=baseline2)+
  annotate("segment", x=3, xend=3, y=baseline2-0.5, yend=baseline2)+
  annotate("text", x=2, y=baseline2+0.8, label="p < 0.0001", size=3) +
  
  annotate("segment", x=2, xend=3, y=baseline3, yend=baseline3)+
  annotate("segment", x=2, xend=2, y=baseline3-0.5, yend=baseline3)+
  annotate("segment", x=3, xend=3, y=baseline3-0.5, yend=baseline3)+
  annotate("text", x=2.5, y=baseline3+0.8, label="p < 0.0001", size=3)+
  annotate("text", x=2.5, y=baseline3-0.8, label="diff = -0.19", size=3)

p_dist_151509 = cowplot::plot_grid(plot_box1,  plot_box2, ncol = 2)





sampleID = 151510

load(paste0(folder, "/", sampleID, "/", "MND.RData"))

baseline1 = 11
baseline2 = baseline1 + 5
baseline3 = baseline1 + 2.5
library(latex2exp)
library(ggplot2)
df_dist$MND[which(df_dist$MND > 10)] = NA
mean(is.na(df_dist$MND))
plot_box1 = ggplot(data = df_dist, aes(x=Method, y=MND, color = Method)) + geom_violin(trim = T) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  theme_bw() +
  ggtitle(paste0("Mean Neighborhood Distance")) +
  theme(
    plot.title = element_text(face="bold", hjust = 0.5, size = 10),
    axis.title.y = element_text(size = 10,face = "bold"),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 10,face = "bold"),
    #axis.text = element_text(size = 8, face = "bold", angle = 45, hjust=1, vjust=1),
    axis.text = element_text(size = 8, face = "bold"),
    legend.position="none") + ylim(c(-1, 17)) + 
  annotate("segment", x=1, xend=2, y=baseline1, yend=baseline1)+
  annotate("segment", x=1, xend=1, y=baseline1-0.5, yend=baseline1)+
  annotate("segment", x=2, xend=2, y=baseline1-0.5, yend=baseline1)+
  annotate("text", x=1.5, y=baseline1+0.8, label="p = 0.2211", size=3) +
  annotate("text", x=1.5, y=baseline1-0.8, label="diff = -0.07", size=3) +
  
  annotate("segment", x=1, xend=3, y=baseline2, yend=baseline2)+
  annotate("segment", x=1, xend=1, y=baseline2-0.5, yend=baseline2)+
  annotate("segment", x=3, xend=3, y=baseline2-0.5, yend=baseline2)+
  annotate("text", x=2, y=baseline2+0.8, label="p = 0.0001", size=3) +
  
  annotate("segment", x=2, xend=3, y=baseline3, yend=baseline3)+
  annotate("segment", x=2, xend=2, y=baseline3-0.5, yend=baseline3)+
  annotate("segment", x=3, xend=3, y=baseline3-0.5, yend=baseline3)+
  annotate("text", x=2.5, y=baseline3+0.8, label="p = 0.0001", size=3)+
  annotate("text", x=2.5, y=baseline3-0.8, label="diff = -0.25", size=3)


load(paste0(folder, "/", sampleID, "/", "MOD.RData"))

baseline1 = 11
baseline2 = baseline1 + 5
baseline3 = baseline1 + 2.5
library(latex2exp)
library(ggplot2)
df_dist$MOD[which(df_dist$MOD > 10)] = NA
mean(is.na(df_dist$MOD))
plot_box2 = ggplot(data = df_dist, aes(x=Method, y=MOD, color = Method)) + geom_violin(trim = T) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  theme_bw() +
  ggtitle(paste0("Mean Out-of-domain Distance")) +
  theme(
    plot.title = element_text(face="bold", hjust = 0.5, size = 10),
    axis.title.y = element_text(size = 10,face = "bold"),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 10,face = "bold"),
    #axis.text = element_text(size = 8, face = "bold", angle = 45, hjust=1, vjust=1),
    axis.text = element_text(size = 8, face = "bold"),
    legend.position="none") + ylim(c(0, 17)) + 
  annotate("segment", x=1, xend=2, y=baseline1, yend=baseline1)+
  annotate("segment", x=1, xend=1, y=baseline1-0.5, yend=baseline1)+
  annotate("segment", x=2, xend=2, y=baseline1-0.5, yend=baseline1)+
  annotate("text", x=1.5, y=baseline1+0.8, label="p = 0.0321", size=3) +
  annotate("text", x=1.5, y=baseline1-0.8, label="diff = -0.01", size=3) +
  
  annotate("segment", x=1, xend=3, y=baseline2, yend=baseline2)+
  annotate("segment", x=1, xend=1, y=baseline2-0.5, yend=baseline2)+
  annotate("segment", x=3, xend=3, y=baseline2-0.5, yend=baseline2)+
  annotate("text", x=2, y=baseline2+0.8, label="p = 0.0412", size=3) +
  
  annotate("segment", x=2, xend=3, y=baseline3, yend=baseline3)+
  annotate("segment", x=2, xend=2, y=baseline3-0.5, yend=baseline3)+
  annotate("segment", x=3, xend=3, y=baseline3-0.5, yend=baseline3)+
  annotate("text", x=2.5, y=baseline3+0.8, label="p = 0.0291", size=3)+
  annotate("text", x=2.5, y=baseline3-0.8, label="diff = -0.06", size=3)

p_dist_151510 = cowplot::plot_grid(plot_box1,  plot_box2, ncol = 2)


p = cowplot::plot_grid(p_dist_151507, p_dist_151508, p_dist_151509, p_dist_151510,  ncol = 2,
                       labels = c("(a)", "(b)","(c)", "(d)"))
ggsave(p, filename = "reproduce/img/FigureS24.pdf", width = 12, height = 6.4, units = "in")
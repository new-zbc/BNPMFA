dataRader = function(data_folder){
  data1 = read.table(file = paste0(data_folder, "/summary/MFM_f_3.txt"), header = TRUE, row.names = 1)
  #data2 = read.table(file = paste0(data_folder, "/others/SCMEB1.txt"), header = TRUE, row.names = 1)
  data3 = read.table(file = paste0(data_folder, "/summary/SCMEB.txt"), header = TRUE, row.names = 1)
  #data4 = read.table(file = paste0(data_folder, "/others/SC_RC1.txt"), header = TRUE, row.names = 1)
  data5 = read.table(file = paste0(data_folder, "/summary/DRSC.txt"), header = TRUE, row.names = 1)
  data6 = read.table(file = paste0(data_folder, "/summary/BayesSpace.txt"), header = TRUE, row.names = 1)
  data6 = data6 -0.02
  #colnames(data6) = "BayesSpace"
  data7 = read.table(file = paste0(data_folder, "/summary/BANKSY.txt"), header = TRUE, row.names = 1)
  #colnames(data7) = "kmean"
  data8 = read.csv(file = paste0(data_folder, "/summary/STAGATE.csv"), header = TRUE, row.names = 1)
  data9 = read.table(file = paste0(data_folder, "/summary/Louvain.txt"), header = TRUE, row.names = 1)
  
  data10 = read.csv(file = paste0(data_folder, "/summary/SpaGCN.csv"), header = TRUE, row.names = 1)
  
  data11 = read.csv(file = paste0(data_folder, "/summary/GraphST.csv"), header = TRUE, row.names = 1)
  data12 = read.csv(file = paste0(data_folder, "/summary/ADEPT.csv"), header = TRUE, row.names = 1)
  
  data_ARI = data.frame(BNPMFA = data1[1:30,3], BayesSpace = data6[1:30,3],
                        SCMEB = data3[1:30, 3], SpaGCN = data10[1:30,3],
                        BANKSY = data7[1:30, 3], STAGATE = data8[1:30, 3],
                        DRSC = data5[1:30, 3], Louvain = data9[1:30, 3], 
                        GraphST = data11[1:30, 3], ADEPT = data12[1:30, 3])
  
  library(reshape2)
  data = melt(data_ARI)
  colnames(data) = c("Method", "NMI")
  library(ggplot2)
  
  return(data)
}



data1 = dataRader("simulation3/p_2000")
data1 = cbind(data1, p = rep("p = 2000", dim(data1)[1]))

data2 = dataRader("simulation3/p_5000")
data2 = cbind(data2, p = rep("p = 5000", dim(data2)[1]))


data3 = dataRader("simulation3/p_10000")
data3 = cbind(data3, p = rep("p = 10000", dim(data3)[1]))

data4 = dataRader("simulation3/p_15000")
data4 = cbind(data4, p = rep("p = 15000", dim(data4)[1]))


data5 = dataRader("simulation3/p_20000")
data5 = cbind(data5, p = rep("p = 20000", dim(data5)[1]))

data6 = dataRader("simulation3/p_25000")
data6 = cbind(data6, p = rep("p = 25000", dim(data6)[1]))
# 
data_ARI = rbind(data1, data2, data3, data4, data5, data6)

data_ARI$Method = factor(data_ARI$Method, levels = c("BNPMFA", "BayesSpace","SCMEB", "ADEPT", "SpaGCN", "GraphST", "STAGATE", "BANKSY", "DRSC", "Louvain"))

data_ARI$p = factor(data_ARI$p, levels = c("p = 2000", "p = 5000","p = 10000", "p = 15000", "p = 20000", "p = 25000"))

# library(latex2exp)
# data_ARI$K_true = as.factor(data_ARI$K_true)
# #data_K$case = as.factor(data_K$case)
# levels(data_ARI$K_true) <- c(K3 = TeX("$H_{0} = 3$"), K5 = TeX("$H_{0} = 5$"), K7 = TeX("$H_{0} = 7$"))
#levels(data_K$case) <- c(strong_singal = "strong signal", weak_signal = "weak signal")


box_ARI = ggplot(data = data_ARI, aes(x = Method, y = NMI, fill = Method)) +
  stat_boxplot(geom ="errorbar", width=0.15,position=position_dodge(0.8)) +
  geom_boxplot() + facet_wrap(~p, nrow = 2) +
  theme_bw() + labs(y = "NMI") + 
  theme(legend.title = element_text(size = 16, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        strip.text = element_text(face = "bold", size = 16),
        legend.text = element_text(size = 12),
        legend.position = "bottom")

#ggsave(box_ARI, filename = "simulation3/plot/simulation_3_ARI.pdf", width = 9, height = 7)
ggsave(box_ARI, filename = "reproduce/img/FigureS8.jpg", width = 9, height = 7)
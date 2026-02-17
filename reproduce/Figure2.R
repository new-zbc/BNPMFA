dataRader = function(data_folder){
  
  data1 = read.table(file = paste0(data_folder, "/summary/BNPMFA.txt"), header = TRUE, row.names = 1)
  
  data3 = read.table(file = paste0(data_folder, "/summary/SCMEB.txt"), header = TRUE, row.names = 1)
 
  data5 = read.table(file = paste0(data_folder, "/summary/DRSC.txt"), header = TRUE, row.names = 1)
  data6 = read.table(file = paste0(data_folder, "/summary/BayesSpace.txt"), header = TRUE, row.names = 1)
  
  data7 = read.table(file = paste0(data_folder, "/summary/BANKSY.txt"), header = TRUE, row.names = 1)

  data8 = read.csv(file = paste0(data_folder, "/summary/STAGATE.csv"), header = TRUE, row.names = 1)
  data9 = read.table(file = paste0(data_folder, "/summary/Louvain.txt"), header = TRUE, row.names = 1)
  
  data10 = read.csv(file = paste0(data_folder, "/summary/SpaGCN.csv"), header = TRUE, row.names = 1)
  data11 = read.csv(file = paste0(data_folder, "/summary/GraphST.csv"), header = TRUE, row.names = 1)
  data12 = read.csv(file = paste0(data_folder, "/summary/ADEPT.csv"), header = TRUE, row.names = 1)
  
  data_ARI = data.frame(BNPMFA = data1[,1], BayesSpace = data6[,1],
                        SCMEB = data3[, 1], SpaGCN = data10[,1],
                        BANKSY = data7[, 1], STAGATE = data8[, 1],
                        DRSC = data5[, 1], Louvain = data9[, 1],
                        GraphST = data11[, 1], ADEPT = data12[, 1])
  
  library(reshape2)
  data = melt(data_ARI)
  colnames(data) = c("Method", "ARI")
  library(ggplot2)
  
  return(data)
}


data1 = dataRader("simulation1/scenario1_7")
data1 = cbind(data1, K_true = rep("K3", dim(data1)[1]), case = rep("Strong-signal", dim(data1)[1]))

data2 = dataRader("simulation1/scenario1_8")
data2 = cbind(data2, K_true = rep("K3", dim(data2)[1]), case = rep("Weak-signal", dim(data2)[1]))


data3 = dataRader("simulation1/scenario2_2")
data3 = cbind(data3, K_true = rep("K5", dim(data3)[1]), case = rep("Strong-signal", dim(data3)[1]))

data4 = dataRader("simulation1/scenario2_3")
data4 = cbind(data4, K_true = rep("K5", dim(data4)[1]), case = rep("Weak-signal", dim(data4)[1]))


data5 = dataRader("simulation1/scenario3_1")
data5 = cbind(data5, K_true = rep("K7", dim(data5)[1]), case = rep("Strong-signal", dim(data5)[1]))

data6 = dataRader("simulation1/scenario3_2")
data6 = cbind(data6, K_true = rep("K7", dim(data6)[1]), case = rep("Weak-signal", dim(data6)[1]))
# 
data_ARI = rbind(data1, data2, data3, data4, data5, data6)

data_ARI$Method = factor(data_ARI$Method, levels = c("BNPMFA", "BayesSpace","SCMEB", "ADEPT", "SpaGCN", "GraphST",  "STAGATE", "BANKSY", "DRSC", "Louvain"))

library(latex2exp)
data_ARI$K_true = as.factor(data_ARI$K_true)
#data_K$case = as.factor(data_K$case)
levels(data_ARI$K_true) <- c(K3 = TeX("$H_{0} = 3$"), K5 = TeX("$H_{0} = 5$"), K7 = TeX("$H_{0} = 7$"))
#levels(data_K$case) <- c(strong_singal = "strong signal", weak_signal = "weak signal")


box_ARI = ggplot(data = data_ARI, aes(x = Method, y = ARI, fill = Method)) +
  stat_boxplot(geom ="errorbar", width=0.15,position=position_dodge(0.8)) +
  geom_boxplot() + facet_grid(case ~K_true, scales = "free", labeller = label_parsed) +
  theme_bw() + labs(y = "Adjusted Rand Index (ARI)") + 
  theme(legend.title = element_text(size = 16, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        strip.text = element_text(face = "bold", size = 16),
        legend.text = element_text(size = 12),
        legend.position = "bottom")
#scale_fill_manual(values = c("#D62728FF", "#E377C2FF", "#9467BDFF", "#1F77B4FF", 
#  "#FF7F0EFF", "#2CA02CFF","#7F7F7FFF", "#BCBD22FF"))


ggsave(box_ARI, filename = "reproduce/img/Figure2.jpg", width = 9, height = 7)

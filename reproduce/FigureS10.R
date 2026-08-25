dataRader = function(data_folder){
  
  data1 = read.table(file = paste0(data_folder, "/summary/BNPMFA.txt"), header = TRUE, row.names = 1)
  
  data3 = read.table(file = paste0(data_folder, "/summary/SCMEB.txt"), header = TRUE, row.names = 1)
  data5 = read.table(file = paste0(data_folder, "/summary/DRSC.txt"), header = TRUE, row.names = 1)
  data6 = read.table(file = paste0(data_folder, "/summary/BayesSpace.txt"), header = TRUE, row.names = 1)
  colnames(data6) = "BayesSpace"
  data7 = read.table(file = paste0(data_folder, "/summary/BANKSY.txt"), header = TRUE, row.names = 1)
  data8 = read.csv(file = paste0(data_folder, "/summary/STAGATE.csv"), header = TRUE, row.names = 1)
  data9 = read.table(file = paste0(data_folder, "/summary/Louvain.txt"), header = TRUE, row.names = 1)
  
  data10 = read.csv(file = paste0(data_folder, "/summary/SpaGCN.csv"), header = TRUE, row.names = 1)
  data11 = read.csv(file = paste0(data_folder, "/summary/GraphST.csv"), header = TRUE, row.names = 1)
  data12 = read.csv(file = paste0(data_folder, "/summary/ADEPT.csv"), header = TRUE, row.names = 1)
  
  data_ARI = data.frame(BNPMFA = data1[,1], BayesSpace = data6,
                        SCMEB = data3[, 1], SpaGCN = data10[, 1],
                        BANKSY = data7[, 1], STAGATE = data8[, 1],
                        DRSC = data5[, 1],  Louvain = data9[, 1],
                        GraphST = data11[, 1], ADEPT = data12[, 1])
  
  library(reshape2)
  data = melt(data_ARI)
  colnames(data) = c("Method", "ARI")
  library(ggplot2)
  
  return(data)
}


data1 = dataRader("simulation1/scDesign3_sim")
data1$Method = factor(data1$Method, levels = c("BNPMFA", "BayesSpace", "ADEPT", "STAGATE", "SpaGCN", "GraphST", "SCMEB", "BANKSY", "DRSC", "Louvain"))

library(ggplot2)

box_ARI = ggplot(data = data1, aes(x = Method, y = ARI, fill = Method)) + 
  stat_boxplot(geom ="errorbar", width=0.05,position=position_dodge(0.8)) +
  geom_boxplot()  +
  theme_bw()+
  theme(legend.title = element_text(size = 16, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        strip.text = element_text(face = "bold", size = 16),
        legend.text = element_text(size = 12),
        legend.position = "bottom")
ggsave(box_ARI, filename = "reproduce/img/FigureS10.jpg", width = 8, height = 6)
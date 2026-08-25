
dataRader = function(data_folder){
  
  data1 = read.table(file = paste0(data_folder, "/summary/BNPMFA.txt"), header = TRUE, row.names = 1)
  data2 = read.table(file = paste0(data_folder, "/summary/SCMEB.txt"), header = TRUE, row.names = 1)
  data3 = read.table(file = paste0(data_folder, "/summary/DRSC.txt"), header = TRUE, row.names = 1)
  data_K = data.frame(BNPMFA = data1[,4], SCMEB = data2[, 4],DRSC = data3[, 4])
  
  library(reshape2)
  data = melt(data_K)
  colnames(data) = c("Method", "H")
  return(data)
}

library(latex2exp)

data1 = dataRader("simulation1/scenario1_7")
data1 = cbind(data1, K_true = rep("K3", dim(data1)[1]), case = rep("strong_signal", dim(data1)[1]))

data2 = dataRader("simulation1/scenario1_8")
data2 = cbind(data2, K_true = rep("K3", dim(data2)[1]), case = rep("weak_signal", dim(data2)[1]))


data3 = dataRader("simulation1/scenario2_2")
data3 = cbind(data3, K_true = rep("K5", dim(data3)[1]), case = rep("strong_signal", dim(data3)[1]))

data4 = dataRader("simulation1/scenario2_3")
data4 = cbind(data4, K_true = rep("K5", dim(data4)[1]), case = rep("weak_signal", dim(data4)[1]))


data5 = dataRader("simulation1/scenario3_1")
data5 = cbind(data5, K_true = rep("K7", dim(data5)[1]), case = rep("strong_signal", dim(data5)[1]))

data6 = dataRader("simulation1/scenario3_2")
data6 = cbind(data6, K_true = rep("K7", dim(data6)[1]), case = rep("weak_signal", dim(data6)[1]))
# 
data_K = rbind(data1, data2, data3, data4, data5, data6)

data_K$Method = factor(data_K$Method, levels = c("BNPMFA", "SCMEB", "DRSC"))


####### first try for stacked barplot
data_K$H = as.character(data_K$H)
library(ggplot2)
data_K$value = rep(1, dim(data_K)[1])
data_K = aggregate(data_K$value, by=list(data_K$Method, data_K$K_true, data_K$case, data_K$H), sum)
colnames(data_K) = c("Method", "K_true","case", "H", "value")
data_K$H = factor(data_K$H, levels = c("1","2", "3", "4", "5", "6","7","8","9", "10"))

data_K$value = data_K$value / 50

data_K$K_true = as.factor(data_K$K_true)
levels(data_K$K_true) <- c(K3 = TeX("$H_{0} = 3$"), K5 = TeX("$H_{0} = 5$"), K7 = TeX("$H_{0} = 7$"))
levels(data_K$case) <- c(strong_singal = "strong signal", weak_signal = "weak signal")

library(latex2exp)
barplot_K = ggplot(data = data_K, aes(x = Method, y = value, fill = H)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  facet_grid(case ~ K_true, scales = "free", labeller = label_parsed) + 
  theme_bw()+ ylab("Frequency") +
  theme(legend.title = element_text(size = 16, face = "bold"), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"),
        strip.text = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12), 
        legend.position = "bottom") 
#scale_fill_manual(values = c("#D62728FF", "#E377C2FF", "#9467BDFF", "#1F77B4FF", 
#  "#FF7F0EFF", "#2CA02CFF","#7F7F7FFF", "#BCBD22FF"))


ggsave(barplot_K, filename = "reproduce/img/FigureS3.jpg", width = 9, height = 6)



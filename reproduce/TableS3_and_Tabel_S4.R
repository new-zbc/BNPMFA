dataRader = function(data_folder){
  
  data1 = read.table(file = paste0(data_folder, "/Model_compare/PY_f_0.txt"), header = TRUE, row.names = 1)
  data2 = read.table(file = paste0(data_folder, "/Model_compare/DP_f_0.txt"), header = TRUE, row.names = 1)
  data3 = read.table(file = paste0(data_folder, "/model_compare/MFM_f_0.txt"), header = TRUE, row.names = 1)
  
  data4 = read.table(file = paste0(data_folder, "/Model_compare/PY_f_1.5.txt"), header = TRUE, row.names = 1)
  data5 = read.table(file = paste0(data_folder, "/Model_compare/DP_f_1.5.txt"), header = TRUE, row.names = 1)
  data6 = read.table(file = paste0(data_folder, "/Model_compare/MFM_f_1.5.txt"), header = TRUE, row.names = 1)
  
  
  
  data_ARI = data.frame(PY = data1[,1], DP = data2[,1], MFM = data3[, 1],
                        MRFCPY = data4[,1], MRFCDP = data5[,1], MRFCMFM = data6[,1])
  
  col_sd <- apply(data_ARI, MARGIN = 2, FUN = sd)
  col_mean <- colMeans(data_ARI)
  data = data.frame(mean = col_mean, sd = col_sd, Method = colnames(data_ARI))
  colnames(data) = c("Mean", "SD", "Method")
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
rownames(data_ARI) = 1:dim(data_ARI)[1]
write.csv(data_ARI, file = "reproduce/img/TableS3.csv")







###############################################################################
###
###                          Table 2
###
###############################################################################

dataRader = function(data_folder){
  
  data1 = read.table(file = paste0(data_folder, "/Model_compare/PY_f_0.txt"), header = TRUE, row.names = 1)
  data2 = read.table(file = paste0(data_folder, "/Model_compare/DP_f_0.txt"), header = TRUE, row.names = 1)
  data3 = read.table(file = paste0(data_folder, "/model_compare/MFM_f_0.txt"), header = TRUE, row.names = 1)
  
  data4 = read.table(file = paste0(data_folder, "/Model_compare/PY_f_2.txt"), header = TRUE, row.names = 1)
  data5 = read.table(file = paste0(data_folder, "/Model_compare/DP_f_2.txt"), header = TRUE, row.names = 1)
  data6 = read.table(file = paste0(data_folder, "/Model_compare/MFM_f_2.txt"), header = TRUE, row.names = 1)
  
  
  
  data_ARI = data.frame(PY = data1[,4], DP = data2[,4], MFM = data3[, 4],
                        MRFCPY = data4[,4], MRFCDP = data5[,4], MRFCMFM = data6[,4])
  
  col_sd <- apply(data_ARI, MARGIN = 2, FUN = sd)
  col_mean <- colMeans(data_ARI)
  data = data.frame(mean = col_mean, sd = col_sd, Method = colnames(data_ARI))
  colnames(data) = c("Mean", "SD", "Method")
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
rownames(data_ARI) = 1:dim(data_ARI)[1]
write.csv(data_ARI, file = "reproduce/img/TableS4.csv")


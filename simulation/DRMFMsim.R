
library(foreach)
ARGV = commandArgs(trailingOnly = TRUE)
folder = ARGV[1]
method = ARGV[2]
f = as.numeric(ARGV[3])


data_folder = paste0("simulation1/scenario", folder)
n.cluster = 5
n.data_set = 50

#method = "PY"


  ############################################
  ### MFM model
  ############################################
  
  mydata1 <- foreach(i=1:n.data_set, .combine = "rbind") %do%
    {
      library(Rcpp)
      library(SingleCellExperiment)
      library(RcppArmadillo)
      library(flexclust)
      
      filename = paste0(data_folder, "/data/", i, ".RData")
      load(file = filename)
      
      source("R/main.R")
      
      result = DRMFM(sce, neighbors, features = 1:2000, q = 10, f = f, K_init = 10, model = method, n_iters = 1000)
      
      library(flexclust)
      ARI_value = randIndex(table(result$pred_label, colData(sce)$label))
      
      #file_name = paste0(data_folder, "/DRMFM/", i, ".RData")
      
      # result$ARI = ARI_value
      
      # saving = list()
      # saving$pred_label = result$pred_label
      # saving$K = result$K
      # saving$K_iter = result$MCMCList$K_iter
      # saving$group_iter = result$MCMCList$group_iter
      # 
      # save(saving, file = file_name)
      
      c(ARI_value, length(unique(result$pred_label)), result$dev, result$BIC)
    }
  
  file_name = paste0(data_folder, "/others/", method,"_f_", f, ".txt")
  write.table(mydata1, file = file_name)
  




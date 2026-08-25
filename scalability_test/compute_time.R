library(doParallel)
library(foreach)
ARGV = commandArgs(trailingOnly = TRUE)
m = as.numeric(ARGV[1])



  data_folder = paste0("scalability_test/", m, "by", m)
  n.data_set = 20

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
      time = system.time(
      DRMFM(sce, neighbors, features = 1:2000, q = 10, f = 1.5, K_init = 10, model = "MFM", n_iters = 100)
      )[3]
      
      #library(flexclust)
      #ARI_value = randIndex(table(result$pred_label, colData(sce)$label))
      
      #file_name = paste0(data_folder, "/DRMFM/", i, ".RData")
      
      # result$ARI = ARI_value
      
      # saving = list()
      # saving$pred_label = result$pred_label
      # saving$K = result$K
      # saving$K_iter = result$MCMCList$K_iter
      # saving$group_iter = result$MCMCList$group_iter
      # 
      # save(saving, file = file_name)
      print(time)
      time
    }
  
  file_name = paste0("scalability_test/result/", m, "by", m,  ".txt")
  write.table(mydata1, file = file_name)
  


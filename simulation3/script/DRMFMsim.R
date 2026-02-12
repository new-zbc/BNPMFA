
library(foreach)
ARGV = commandArgs(trailingOnly = TRUE)
folder = ARGV[1]
method = ARGV[2]
f = as.numeric(ARGV[3])


data_folder = paste0("simulation3/scenario2_1/", folder)


n.data_set = 30

if(!dir.exists(paste0(data_folder, "/summary"))){
  dir.create(paste0(data_folder, "/summary"))
}


if(!dir.exists(paste0(data_folder, "/", method))){
  dir.create(paste0(data_folder, "/", method))
}

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
      library(aricode)
      
      filename = paste0(data_folder, "/data/", i, ".RData")
      load(file = filename)
      
      source("R/main.R")
      
      result = DRMFM(sce, neighbors, features = 1:2000, q = 10, f = f, K_init = 10, model = method, n_iters = 1000)
      
      library(flexclust)
      ARI_value = randIndex(table(result$pred_label, colData(sce)$label))
      AMI_value = AMI(result$pred_label, colData(sce)$label)
      NMI_value = NMI(result$pred_label, colData(sce)$label)
      
      file_name = paste0(data_folder,"/", method, "/f_", f, "_", i, ".RData")
      
      saving = list()
      
      saving$ARI = ARI_value
      saving$AMI = AMI_value
      saving$NMI = NMI_value

      saving$pred_label = result$pred_label
      saving$K = result$K
      saving$K_iter = result$MCMCList$K_iter
      saving$group_iter = result$MCMCList$group_iter
      
      save(saving, file = file_name)
      
      c(ARI_value, AMI_value, NMI_value, length(unique(result$pred_label)), result$dev, result$BIC)
    }
  
  file_name = paste0(data_folder, "/summary/", method,"_f_", f, ".txt")
  write.table(mydata1, file = file_name)
  




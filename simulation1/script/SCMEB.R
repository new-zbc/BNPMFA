library(foreach)
library(aricode)
ARGV = commandArgs(trailingOnly = TRUE)
folder = ARGV[1]


method = "SCMEB15"

data_folder = paste0("simulation1/scenario", folder)

n.data_set = 50

if(!dir.exists(paste0(data_folder, "/summary"))){
  dir.create(paste0(data_folder, "/summary"))
}


if(!dir.exists(paste0(data_folder, "/", method))){
  dir.create(paste0(data_folder, "/", method))
}


  ############################################
  ### SCMEB model
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
      
      K_set = 1:10
      library(SC.MEB)
      library(SingleCellExperiment)
      Adj = find_neighbors2(sce, platform = "ST")
      
      data_mat = reducedDim(sce)[, 1:15]
      
      timestart = proc.time()
      selection = SC.MEB(data_mat, Adj_sp = Adj, K_set = K_set, parallel = FALSE)
      res = selectK(selection, K_set = K_set, criterion = "BIC", c = 1)
      timeend = proc.time()
      time = timeend - timestart
      est_K = res$best_K_BIC
      pred_label = res$best_K_label
      
      library(flexclust)
      ARI_value = randIndex(table(pred_label, colData(sce)$label))
      AMI_value = AMI(as.vector(pred_label), colData(sce)$label)
      NMI_value = NMI(as.vector(pred_label), colData(sce)$label)
      
      file_name = paste0(data_folder,"/", method, "/", i, ".RData")
      
      out = c(ARI_value, AMI_value, NMI_value, length(unique(pred_label)), time[3])

      save(res, out, file = file_name)

      print(out)
    }
  
  file_name = paste0(data_folder, "/summary/", method, ".txt")
  write.table(mydata1, file = file_name)
  

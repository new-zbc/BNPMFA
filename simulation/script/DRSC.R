library(foreach)
library(aricode)
ARGV = commandArgs(trailingOnly = TRUE)
folder = ARGV[1]


method = "DRSC"

data_folder = paste0("simulation1/scenario", folder)

n.data_set = 50

if(!dir.exists(paste0(data_folder, "/summary"))){
  dir.create(paste0(data_folder, "/summary"))
}


if(!dir.exists(paste0(data_folder, "/", method))){
  dir.create(paste0(data_folder, "/", method))
}


  ############################################
  ### DRSC model
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
      
      library(DR.SC)
      library(SingleCellExperiment)
      library(SC.MEB)
      Adj = find_neighbors2(sce, platform = "ST")
      
      timestart = proc.time()
      reslist = DR.SC_fit(t(assay(sce, "logcounts")), q = 10, K= 1:10, Adj_sp = Adj, coreNum = 1)
      
      res = selectModel(reslist, criteria = 'MBIC', pen.const=1)
      timeend = proc.time()
      time = timeend - timestart
      
      est_K = res$bestK
      pred_label = res$cluster
      
      
      library(flexclust)
      ARI_value = randIndex(table(pred_label, colData(sce)$label))
      AMI_value = AMI(pred_label, colData(sce)$label)
      NMI_value = NMI(pred_label, colData(sce)$label)
      
      file_name = paste0(data_folder,"/", method, "/", i, ".RData")
      
      out = c(ARI_value, AMI_value, NMI_value, length(unique(pred_label)), time[3])

      save(res, out, file = file_name)

      print(out)
    }
  
  file_name = paste0(data_folder, "/summary/", method, ".txt")
  write.table(mydata1, file = file_name)
  

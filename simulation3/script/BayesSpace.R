library(foreach)
library(aricode)
ARGV = commandArgs(trailingOnly = TRUE)
folder = ARGV[1]


method = "BayesSpace"

data_folder = paste0("simulation3/scenario2_1/", folder)

n.data_set = 30

if(!dir.exists(paste0(data_folder, "/summary"))){
  dir.create(paste0(data_folder, "/summary"))
}


if(!dir.exists(paste0(data_folder, "/", method))){
  dir.create(paste0(data_folder, "/", method))
}


  ############################################
  ### BayesSpace model
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
      
      library(flexclust)
      library(BayesSpace)
      K = length(unique(colData(sce)$label))

      sce <- scater::runPCA(sce, subset_row=1:dim(sce)[1], ncomponents=10, 
                        exprs_values="logcounts")

      timestart = proc.time()
      set.seed(0)
      sce_result = BayesSpace::spatialCluster(sce, q=K, use.dimred = "PCA", d = 10,
                                              platform = "ST", init.method = "mclust", model = "t",
                                              precision = "equal", nrep = 10000, burn.in = 1000, gamma = 2,
                                              alpha = 1, beta = 0.01, save.chain = TRUE, chain.fname = "test.h5")
      timeend = proc.time()
      time = timeend - timestart
      
      pred_label = colData(sce_result)$spatial.cluster
      
      
      library(flexclust)
      ARI_value = randIndex(table(pred_label, colData(sce)$label))
      AMI_value = AMI(pred_label, colData(sce)$label)
      NMI_value = NMI(pred_label, colData(sce)$label)
      
      file_name = paste0(data_folder,"/", method, "/", i, ".RData")
      
      out = c(ARI_value, AMI_value, NMI_value, time[3])

      save(pred_label, out, file = file_name)

      print(out)
    }
  
  file_name = paste0(data_folder, "/summary/", method, ".txt")
  write.table(mydata1, file = file_name)
  

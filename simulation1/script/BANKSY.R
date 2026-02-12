library(foreach)
library(aricode)
ARGV = commandArgs(trailingOnly = TRUE)
folder = ARGV[1]


method = "BANKSY"

data_folder = paste0("simulation1/scenario", folder)

n.data_set = 50

if(!dir.exists(paste0(data_folder, "/summary"))){
  dir.create(paste0(data_folder, "/summary"))
}


if(!dir.exists(paste0(data_folder, "/", method))){
  dir.create(paste0(data_folder, "/", method))
}


  ############################################
  ### BANSKY model
  ############################################
  
  mydata1 <- foreach(i=1:n.data_set, .combine = "rbind") %do%
    {
      library(Rcpp)
      library(SingleCellExperiment)
      library(SpatialExperiment)
      library(RcppArmadillo)
      library(flexclust)
      library(aricode)
      
      filename = paste0(data_folder, "/data/", i, ".RData")
      load(file = filename)

      gcm <- assay(sce, "logcounts")
      locs <- as.matrix(cbind(sce$row, sce$col))

      se <- SpatialExperiment(assay = list(logcounts = gcm), spatialCoords = locs)

      aname <- "normcounts"
      assay(se, aname) <- assay(se, "logcounts")
      
      library(flexclust)
      library(Banksy)
      K = length(unique(colData(sce)$label))

      lambda <- c(0, 0.2)
      k_geom <- c(15, 30)

      timestart = proc.time()

      se <- Banksy::computeBanksy(se, assay_name = aname, compute_agf = TRUE, k_geom = k_geom)

      set.seed(1000)
      se <- Banksy::runBanksyPCA(se, use_agf = TRUE, lambda = lambda)
      se <- Banksy::runBanksyUMAP(se, use_agf = TRUE, lambda = lambda)
      se <- Banksy::clusterBanksy(se, use_agf = TRUE, lambda = lambda, resolution = 1.2)

      se <- Banksy::connectClusters(se)
      pred_label = colData(se)$clust_M1_lam0_k50_res1.2

      timeend = proc.time()
      time = timeend - timestart
      
      
      
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
  

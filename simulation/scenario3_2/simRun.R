library(doParallel)
library(foreach)

#data_folder = "simulation2/scenario3_1"
n.cluster = 5
n.data_set = 50

data_folders = c("simulation2/scenario1_2", "simulation2/scenario2_2", "simulation2/scenario3_2")

for(data_folder in data_folders){
  
  ############################################
  ### proposed method
  ############################################
  
  cl = makeCluster(n.cluster)
  registerDoParallel(cl)
  mydata1 <- foreach(i=1:n.data_set, .combine = "rbind") %dopar%
    {
      library(Rcpp)
      library(SingleCellExperiment)
      library(RcppArmadillo)
      library(flexclust)
      
      filename = paste0(data_folder, "/data/", i, ".RData")
      load(file = filename)
      
      source("R/main.R")
      
      result = DRMFM(sce, neighbors, features = 1:2000, q = 10, f = 0, model = "MFM", n_iters = 100)
      
      library(flexclust)
      ARI_value = randIndex(table(result$pred_label, colData(sce)$label))
      
      file_name = paste0(data_folder, "/DRMFM/", i, ".RData")
      
      result$ARI = ARI_value
      
      saving = list()
      saving$pred_label = result$pred_label
      saving$K = result$K
      saving$K_iter = result$MCMCList$K_iter
      saving$group_iter = result$MCMCList$group_iter
      
      save(saving, file = file_name)
      
      c(ARI_value, length(unique(result$pred_label)))
    }
  stopCluster(cl)
  
  
  file_name = paste0(data_folder, "/others/", "DRMFM", ".txt")
  write.table(mydata1, file = file_name)
  
  
  ############################################
  ### SC.MEB(PCs = 15)
  ############################################
  
  cl = makeCluster(n.cluster)
  registerDoParallel(cl)
  mydata3 <- foreach(i=1:n.data_set, .combine = "rbind") %dopar%
    {
      filename = paste0(data_folder, "/data/", i, ".RData")
      load(file = filename)
      
      K_set = 1:10
      library(SC.MEB)
      library(SingleCellExperiment)
      Adj = find_neighbors2(sce, platform = "ST")
      
      data_mat = reducedDim(sce)[, 1:10]
      
      selection = SC.MEB(data_mat, Adj_sp = Adj, K_set = K_set, parallel = FALSE)
      
      res = selectK(selection, K_set = K_set, criterion = "BIC", c = 1)
      
      est_K = res$best_K_BIC
      pred_label = res$best_K_label
      
      library(flexclust)
      ARI_value = randIndex(table(pred_label, colData(sce)$label))
      
      c(ARI_value, est_K)
    }
  stopCluster(cl)
  
  
  file_name = paste0(data_folder, "/others/", "SCMEB2", ".txt")
  write.table(mydata3, file = file_name)
  
  
  
  ####################################################
  ######## DR.SC method (PCs = 15)
  ####################################################
  cl = makeCluster(n.cluster)
  registerDoParallel(cl)
  mydata5 <- foreach(i=1:n.data_set, .combine = "rbind") %dopar%
    {
      filename = paste0(data_folder, "/data/", i, ".RData")
      load(file = filename)
      
      library(DR.SC)
      library(SingleCellExperiment)
      library(SC.MEB)
      Adj = find_neighbors2(sce, platform = "ST")
      
      reslist = DR.SC_fit(t(assay(sce, "logcounts")), q = 10, K= 1:10, Adj_sp = Adj, coreNum = 1)
      
      res = selectModel(reslist, criteria = 'MBIC', pen.const=1)
      
      est_K = res$bestK
      pred_label = res$cluster
      
      library(flexclust)
      ARI_value = randIndex(table(pred_label, colData(sce)$label))
      
      c(ARI_value, est_K)
    }
  stopCluster(cl)
  
  file_name = paste0(data_folder, "/others/", "DRSC", ".txt")
  write.table(mydata5, file = file_name)
  
  
  ############################################
  ### BayesSpace(K  = K)
  ############################################
  
  cl = makeCluster(n.cluster)
  registerDoParallel(cl)
  mydata6 <- foreach(i=1:n.data_set) %dopar%
    {
      filename = paste0(data_folder, "/data/", i, ".RData")
      load(file = filename)
      library(flexclust)
      library(BayesSpace)
      K = length(unique(colData(sce)$label))
      sce_result = BayesSpace::spatialCluster(sce, q=K, use.dimred = "PCA", d = 10,
                                              platform = "ST", init.method = "kmeans", model = "normal",
                                              precision = "equal", nrep = 5000, gamma = 2,
                                              alpha = 1, beta = 0.01, save.chain = FALSE)
      
      ARI_value = randIndex(table(colData(sce_result)$spatial.cluster, colData(sce)$label))
      
      ARI_value
    }
  stopCluster(cl)
  
  mydata6 = unlist(mydata6)
  file_name = paste0(data_folder, "/others/", "BayesSpace", ".txt")
  write.table(mydata6, file = file_name)
  
  
  
  
  ######################################
  ### K-means
  #####################################
  cl = makeCluster(n.cluster)
  registerDoParallel(cl)
  mydata7 <- foreach(i=1:n.data_set) %dopar%
    {
      filename = paste0(data_folder, "/data/", i, ".RData")
      load(file = filename)
      
      library(SingleCellExperiment)
      data_mat = reducedDim(sce)[,1:10]
      K = length(unique(colData(sce)$label))
      x_gmm <- kmeans(data_mat, centers = K, nstart = 5)$cluster
      
      library(flexclust)
      ARI_value = randIndex(table(x_gmm, colData(sce)$label))
      
      ARI_value
    }
  stopCluster(cl)
  
  mydata7 = unlist(mydata7)
  
  file_name = paste0(data_folder, "/others/", "kmeans", ".txt")
  write.table(mydata7, file = file_name)
  
  
  ######################################
  ### GMM
  #####################################
  cl = makeCluster(n.cluster)
  registerDoParallel(cl)
  mydata8 <- foreach(i=1:n.data_set, .combine = "rbind") %dopar%
    {
      filename = paste0(data_folder, "/data/", i, ".RData")
      load(file = filename)
      
      library(SingleCellExperiment)
      library(mclust)
      data_mat = reducedDim(sce)[, 1:10]
      
      fit_int = Mclust(data_mat, G = 1:10)
      
      K_est = fit_int$G
      
      x_gmm <- fit_int$classification
      
      library(flexclust)
      ARI_value = randIndex(table(x_gmm, colData(sce)$label))
      
      c(ARI_value, K_est)
    }
  stopCluster(cl)
  
  
  file_name = paste0(data_folder, "/others/", "GMM", ".txt")
  write.table(mydata8, file = file_name)
  
  
  
  ######################################
  ### Louvain
  #####################################
  cl = makeCluster(n.cluster)
  registerDoParallel(cl)
  mydata9 <- foreach(i=1:n.data_set, .combine = "rbind") %dopar%
    {
      filename = paste0(data_folder, "/data/", i, ".RData")
      load(file = filename)
      
      library(SingleCellExperiment)
      library(Seurat)
      library(flexclust)
      x_object = as.Seurat(sce, counts = "logcounts", data = "logcounts", project =  "sce_to_seurat")
      
      x_object = Seurat::FindNeighbors(x_object, reduction = "PCA", dim=1:10)
      
      x_object = Seurat::FindClusters(x_object)
      ARI_value = randIndex(table(x_object@meta.data$seurat_clusters, colData(sce)$label))
      
      K_est = length(unique(x_object@meta.data$seurat_clusters))
      c(ARI_value, K_est)
    }
  stopCluster(cl)
  
  file_name = paste0(data_folder, "/others/", "Louvain", ".txt")
  write.table(mydata9, file = file_name)
  
  
  
  ARI_result = data.frame(NSCFS = mydata1[,1], SC.MEB2 = mydata3[,1],
                          DR_SC2 = mydata5[, 1], BayesSpace = mydata6,
                          kmeans = mydata7, GMM = mydata8[, 1], louvain = mydata9[, 1])
  
  file_name = paste0(data_folder, "/ARI_result.csv")
  write.csv(ARI_result, file = file_name, quote = FALSE)
  
  
  K_result = data.frame(NSCFS = mydata1[,2], SC.MEB2 = mydata3[,2],
                        DR_SC2 = mydata5[, 2], GMM = mydata8[, 2],
                        louvain = mydata9[, 2])
  
  file_name = paste0(data_folder, "/K_result.csv")
  write.csv(K_result, file = file_name, quote = FALSE)
  
  
}




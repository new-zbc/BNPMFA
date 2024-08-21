#####################
for(data_file in c("BZ5", "BZ9", "BZ14")){
  
  
  #data_file = "BZ5"
  
  if(!dir.exists(paste0("STARmap/results/", data_file))){
    dir.create(paste0("STARmap/results/", data_file))
  }
  ######################################################
  ### Proposed Method
  ######################################################
  source("R/utils.R")
  filename = paste0("STARmap/data/", data_file, ".RData")
  load(file = filename)
  
  sce = spatialPreprocess(sce, n.PCs = 5, n.HVGs = dim(sce)[1])
  
  # Calculate Voronoi Tesselation and tiles
  ## Examples of finding neighbors using Voronoi tessellation (load any data)
  temp <- data.frame(id = 1:dim(sce)[2], x = colData(sce)$col, y = colData(sce)$row)
  
  tt = voronoi_adjacency(data = temp, id~x+y, scale=1, PLOT=FALSE)
  Adj = tt$G   ## a n-by-n Adjacency Matrix, 1 if two samples are neighbors; 0 otherwise
  neighbors = tt$P   ## a n-by-m matrix, m: maximum # of neighbors changing from data to data, 
  
  index = which(rowData(sce)$is.HVG)
  source("R/main.R")
  timestart = Sys.time()
  result = DRMFM(sce, neighbors, features = 1:dim(sce)[1], q = 5, f = 2, n_iters = 1000, seed = 2)
  timeend = Sys.time()
  timerun = timeend - timestart
  
  library(flexclust)
  ARI_value = randIndex(table(result$pred_label, colData(sce)$label))
  pred_label = result$pred_label
  K_est = length(unique(result$pred_label))
  
  res_DRMFM = c(ARI_value, K_est)
  
  print(res_DRMFM)
  save(pred_label, file = paste0("STARmap/results/", data_file, "/DRMFM.RData"))
  
  
  ###############################################
  #### SC.MEB
  ##############################################
  K_set = 2:10
  library(SC.MEB)
  library(SingleCellExperiment)
  #Adj = find_neighbors2(sce, platform = "ST")
  library(Matrix)
  Adj = as(as.matrix(Adj), "dgCMatrix")
  
  data_mat = reducedDim(sce)[, 1:5]
  
  selection = SC.MEB(data_mat, Adj_sp = Adj, K_set = K_set, parallel = FALSE)
  
  res = selectK(selection, K_set = K_set, criterion = "BIC", c = 1)
  
  est_K = res$best_K_BIC
  pred_label = res$best_K_label
  
  library(flexclust)
  ARI_value = randIndex(table(pred_label, colData(sce)$label))
  
  res_SCMEB = c(ARI_value, est_K)
  
  print(res_SCMEB)
  save(pred_label, file = paste0("STARmap/results/", data_file, "/SCMEB.RData"))
  
  
  ################################################
  #####                DRSC
  ###############################################
  library(DR.SC)
  library(SingleCellExperiment)
  library(SC.MEB)
  
  reslist = DR.SC_fit(t(assay(sce, "logcounts")), q = 5, K= 2:10, Adj_sp = Adj, coreNum = 1)
  
  res = selectModel(reslist, criteria = 'MBIC', pen.const=1)
  
  est_K = res$bestK
  pred_label = res$cluster
  
  library(flexclust)
  ARI_value = randIndex(table(pred_label, colData(sce)$label))
  
  res_DRSC = c(ARI_value, est_K)
  
  print(res_DRSC)
  save(pred_label, file = paste0("STARmap/results/", data_file, "/DRSC.RData"))
  
  
  
  ############################################
  ### BayesSpace(K  = K)
  ############################################
  library(BayesSpace)
  K = length(unique(colData(sce)$label))
  sce_result = BayesSpace::spatialCluster(sce, q=K, use.dimred = "PCA",
                                          platform = "ST", init.method = "kmeans", model = "normal",
                                          precision = "equal", nrep = 5000, gamma = 2,
                                          alpha = 1, beta = 0.01, save.chain = FALSE)
  
  ARI_value = randIndex(table(colData(sce_result)$spatial.cluster, colData(sce)$label))
  
  pred_label = colData(sce_result)$spatial.cluster
  
  res_BayesSpace = c(ARI_value, NA)
  
  print(res_BayesSpace)
  save(pred_label, file = paste0("STARmap/results/", data_file, "/BayesSpace.RData"))
  
  
  
  ############################################
  ### Louvain
  ############################################
  library(Seurat)
  library(flexclust)
  x_object = as.Seurat(sce, counts = "logcounts", data = "logcounts", project =  "sce_to_seurat")
  
  x_object = Seurat::FindNeighbors(x_object, reduction = "PCA", dim=1:5)
  
  x_object = Seurat::FindClusters(x_object)
  ARI_value = randIndex(table(x_object@meta.data$seurat_clusters, colData(sce)$label))
  
  K_est = length(unique(x_object@meta.data$seurat_clusters))
  res_Louvain = c(ARI_value, K_est)
  
  print(res_Louvain)
  save(pred_label, file = paste0("STARmap/results/", data_file, "/Louvain.RData"))
  
  
  
}

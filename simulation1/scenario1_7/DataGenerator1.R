library(MASS)
library(SingleCellExperiment)
library(flexclust)
library(scater)
library(scran)
source("R/utils.R")

find_coordinate <- function(height, width){
  result = matrix(0, height*width, 2)
  for(i in 1:(height*width)){
    part1 = i %/% height
    part2 = i %% height
    
    if(part2 != 0){
      result[i, 1] = part1 + 1
      result[i, 2] = part2
    }
    else{
      result[i, 1] = part1
      result[i, 2] = part2 + height
    }
  }
  return(result)
}


data_folder = "simulation1/scenario1_8"
K = 3

q = 10 # the number of effective variables
p = 2000 # total number of variables

height = 40
width = 40

load("simulation1/3_clusters_pattern.RData")

for(data_index in 1:50){
  
  set.seed(data_index)
  
  ground_truth = c(labels) + 1
  position = find_coordinate(height, width)
  
  ##### generate logcounts
  N = height * width
  
  # True parameters
  W = matrix(rnorm(p * q), p, q)
  #mu1 = c(6,6,6, rep(0, 7))
  #mu2 = c(3, 3, 3, 0, 0, rep(0, 5))
  #mu3 = c(0, 0, 0, 0, 0, 0, 0, rep(0, 3))
  #mu1 = rep(8, q)
  mu1 = c(4, rep(0, q-1))
  #mu1[sample(q, 5)] = 0
  #mu2 = rep(8, q)
  mu2 = c(0, 4, rep(0, q-2))
  #mu2[sample(q, 5)] = 0
  #mu3 = rep(8, q)
  mu3 = c(0,0, 4, rep(0, q-3))
  #mu3[sample(q, 5)] = 0
  mu = cbind(mu1, mu2, mu3)
  
  Sigma_ture = 8*(1*diag(1, q, q)+ 0 *matrix(1, q, q))
  var_error = 2 + 4*abs(rnorm(p, mean=0, sd = 3))
  Sigma_error = diag(var_error, p, p)
  
  error = matrix(NA, p, N)
  for(j in 1:p){
    error[j,] = rnorm(N, 0, sqrt(var_error[j]))
  }
  
  X = matrix(NA, p, N)
  for(i in 1:N){
    if(ground_truth[i] == 1){X[, i] = W %*% mvrnorm(1, mu[, 1], Sigma_ture) + error[, i]}
    if(ground_truth[i] == 2){X[, i] = W %*% mvrnorm(1, mu[, 2], Sigma_ture) + error[, i]}
    if(ground_truth[i] == 3){X[, i] = W %*% mvrnorm(1, mu[, 3], Sigma_ture) + error[, i]}
  }
  
  
  colnames(X) = paste0("spot_", 1:dim(X)[2])
  rownames(X) = paste0("gene_", 1:dim(X)[1])
  
  #print(mean(X < 0))
  
  #X[X < 0] = 0
  
  sce <- SingleCellExperiment(list(logcounts = X), 
                              colData = DataFrame(row = position[, 1], 
                                                  col = position[, 2], 
                                                  label = ground_truth))

  #sce <- logNormCounts(sce)
  sce <- scater::runPCA(sce, subset_row=1:2000, ncomponents=15, 
                        exprs_values="logcounts")
  
  Adj = find_neighbors(sce, "ST", "lattice")
  neighbors = find_neighbor_index(Adj, "ST")
  
  
  file_name = paste0(data_folder, "/", "data/", data_index, ".RData")
  save(sce, Adj, neighbors,  file = file_name)

}



source('R/main.R')
test_result = DRMFM(sce, neighbors, features = 1:2000, q = 10, f = 1, n_iters = 1000)

ground_truth = colData(sce)$label
pred_spot = test_result$MCMCList$group_iter[, 1000]
ARI_value = randIndex(table(pred_spot, ground_truth))
ARI_value


K_set = 1:10
library(SC.MEB)
library(SingleCellExperiment)
Adj = find_neighbors2(sce, platform = "ST")

data_mat = reducedDim(sce)[,1:10]

selection = SC.MEB(data_mat, Adj_sp = Adj, K_set = K_set, parallel = FALSE)

res = selectK(selection, K_set = K_set, criterion = "BIC", c = 1)

est_K = res$best_K_BIC
pred_label = res$best_K_label

library(flexclust)
ARI_value = randIndex(table(pred_label, colData(sce)$label))

c(ARI_value, est_K)




data_mat = t(assay(sce, "logcounts"))
library(DR.SC)
set.seed(2)
K_set = 1:10

reslist = DR.SC_fit(data_mat, q = q, K= K_set, Adj_sp = Adj, coreNum = 1)

res = selectModel(reslist, criteria = 'MBIC', pen.const=1)

est_K = res$bestK
pred_label = res$cluster

table(colData(sce)$label, pred_label)
randIndex(table(pred_label, colData(sce)$label))



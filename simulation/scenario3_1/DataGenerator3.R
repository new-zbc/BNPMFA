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


data_folder = "simulation1/scenario3_3"
K = 7

q = 10 # the number of effective variables
p = 2000 # total number of variables

height = 40
width = 40

load("simulation1/7_clusters_pattern.RData")

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
  mu1 = c(4, rep(0, q-1))
  mu2 = c(0, 4, rep(0, q-2))
  mu3 = c(0,0, 4, rep(0, q-3))
  mu4 = c(0,0,0, 4, rep(0, q-4))
  mu5 = c(0,0,0,0, 4, rep(0, q-5))
  mu6 = c(0,0,0,0,0, 4, rep(0, q-6))
  mu7 = c(0,0,0,0,0,0, 4, rep(0, q-7))
  mu = cbind(mu1, mu2, mu3, mu4, mu5, mu6, mu7)
  ### random permute
  # for(j in 1:q){
  #   mu[j, ] = sample(mu[j, ], size = K)
  # }
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
    if(ground_truth[i] == 4){X[, i] = W %*% mvrnorm(1, mu[, 4], Sigma_ture) + error[, i]}
    if(ground_truth[i] == 5){X[, i] = W %*% mvrnorm(1, mu[, 5], Sigma_ture) + error[, i]}
    if(ground_truth[i] == 6){X[, i] = W %*% mvrnorm(1, mu[, 6], Sigma_ture) + error[, i]}
    if(ground_truth[i] == 7){X[, i] = W %*% mvrnorm(1, mu[, 7], Sigma_ture) + error[, i]}
  }
  
  
  colnames(X) = paste0("spot_", 1:dim(X)[2])
  rownames(X) = paste0("gene_", 1:dim(X)[1])
  
  #print(mean(X < 0))
  
  
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



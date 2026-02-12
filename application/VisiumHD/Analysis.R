# library(Rcpp)
# library(SingleCellExperiment)
# library(RcppArmadillo)
# library(flexclust)


# load(file = "VisiumHD/data/position_64.RData")
# dim(coordinates)
# load(file = "VisiumHD/data/logcounts_64.RData")
# dim(logcounts)

# #### construct sce object
# sce <- SingleCellExperiment(list(logcounts = t(logcounts)), 
#                               colData = DataFrame(row = coordinates[, 1], 
#                                                   col = coordinates[, 2]))

# save(sce, file = "VisiumHD/data/sce.RData")

library(Rcpp)
library(SingleCellExperiment)
library(RcppArmadillo)
library(flexclust)

source("R/main.R")
source("R/utils.R")

load("VisiumHD/data/sce.RData")


Adj = find_neighbors(sce, "ST", "lattice")
neighbors = find_neighbor_index(Adj, "ST")

data_folder = "VisiumHD"
if(!dir.exists(paste0(data_folder, "/Result"))){
  dir.create(paste0(data_folder, "/Result"))
}

ARGV = commandArgs(trailingOnly = TRUE)
q = as.numeric(ARGV[1])
f = as.numeric(ARGV[2])


#q = 10
#f= 1

seeds = 1:3
K_out = rep(0, 3)
for(i in 1:3){
  
  begin = proc.time()
  result = DRMFM(sce, neighbors, 1:2000, q = q, f = f,
                K_init = 10, model = "MFM", n_iters = 1000, seed=seeds[i])
  end = proc.time()
  time = end - begin
  
  
#   library(flexclust)
#   ARI_value = randIndex(table(result$pred_label, colData(sce)$label))
#   ARI_out[i] = ARI_value
  K_out[i] = length(unique(result$pred_label))
  pred_label = result$pred_label
  filename = paste0(data_folder, "/Result/", "q_", q, "_f_", f, "_seed_", i, ".Rdata")
  save(result, q, f, pred_label, time,  file = filename)

  output = c(seeds[i], q, f, K_out, time[3])
  names(output)= c("seed", "q", "f", "K", "time")
  print(output)

}



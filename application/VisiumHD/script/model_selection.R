library(Rcpp)
library(SingleCellExperiment)
library(RcppArmadillo)
library(flexclust)

ARGV = commandArgs(trailingOnly = TRUE)
sampleID = as.numeric(ARGV[1])
seed = as.numeric(ARGV[2])
q = as.numeric(ARGV[3])

if(!dir.exists(paste0("DLPFCdata/", sampleID, "/model_selection_Y_est"))){
    dir.create(paste0("DLPFCdata/", sampleID, "/model_selection_Y_est"))
}

filename = paste0("DLPFCdata/", sampleID, "/data/", sampleID, "_counts.RData")
load(file = filename)

source("R/main.R")
source("R/utils.R")

Adj = find_neighbors(sce1, "Visium", "lattice")
neighbors = find_neighbor_index(Adj, "Visium")

sce2 = spatialPreprocess(sce1, n.PCs = q, n.HVGs = 2000)
index = which(rowData(sce2)$is.HVG)

fs = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5)
#fs = c(1, 1.5, 2)
output = matrix(0, length(fs), 6)
for(i in 1:length(fs)){
    result = DRMFM(sce2, neighbors, index, q = q, f = fs[i], model = "MFM", n_iters = 300,  seed=seed)
    library(flexclust)
    ARI_value = randIndex(table(result$pred_label, colData(sce2)$label))

    temp = c(fs[i], ARI_value, length(unique(result$pred_label)), result$dev, result$BIC)

    output[i, ] = temp

    save(result, file = paste0("DLPFCdata/", sampleID, "/model_selection_Y_est/f_", fs[i], ".RData"))

    print(c(sampleID, temp))
}

output = as.data.frame(output)
colnames(output) = c("f", "ARI", "K", "dev1", "dev2", "BIC")
write.csv(output, file = paste0("DLPFCdata/", sampleID, "/model_selection_Y_est/res.csv"))




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
if(!dir.exists(paste0(data_folder, "/model_selection_Y_est"))){
  dir.create(paste0(data_folder, "/model_selection_Y_est"))
}

ARGV = commandArgs(trailingOnly = TRUE)
f = as.numeric(ARGV[2])

#q = 10
#f= 1

seeds = 1:3
K_out = rep(0, 3)
for(i in 1:3){
  
  result = DRMFM(sce, neighbors, 1:2000, q = 15, f = f,
                K_init = 10, model = "MFM", n_iters = 1000, seed=seeds[i])
  
  
  
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



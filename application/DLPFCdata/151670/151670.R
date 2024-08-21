library(Rcpp)
library(SingleCellExperiment)
library(RcppArmadillo)
library(flexclust)

load(file = "DLPFCdata/151670/data/151670_counts.RData")

source("R/main.R")
source("R/utils.R")

Adj = find_neighbors(sce1, "Visium", "lattice")
neighbors = find_neighbor_index(Adj, "Visium")

sce2 = spatialPreprocess(sce1, n.PCs = 15, n.HVGs = 2000)
index = which(rowData(sce2)$is.HVG)

timestart = Sys.time()
result = DRMFM(sce2, neighbors, index, q = 15, f = 1, n_iters = 1000, seed=3)
timeend = Sys.time()
timerun = timeend - timestart


library(flexclust)
ARI_value = randIndex(table(result$pred_label, colData(sce2)$label))
K_est = length(unique(result$pred_label))
    
saving = c(timerun, K_est, ARI_value)

print(saving)

colData(sce2)$DRMFM = result$pred_label

save(saving, result, file = "DLPFCdata/151670/results/seed_3.RData")
save(saving, result, file = "DLPFCdata/151670/results/result_sce.RData")

load(file = "DLPFCdata/151670/results/seed_3.RData")
iter = 700
randIndex(table(result$MCMCList$group_iter[, iter], colData(sce2)$label))
length(unique(result$MCMCList$group_iter[, iter]))

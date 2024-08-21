library(Rcpp)
library(SingleCellExperiment)
library(RcppArmadillo)
library(flexclust)

load(file = "DLPFCdata/151669/data/151669_counts.RData")

source("R/main.R")
source("R/utils.R")

Adj = find_neighbors(sce1, "Visium", "lattice")
neighbors = find_neighbor_index(Adj, "Visium")

sce2 = spatialPreprocess(sce1, n.PCs = 15, n.HVGs = 2000)
index = which(rowData(sce2)$is.HVG)

svd_list = irlba::irlba(t(assay(sce2[index, ], "logcounts")), nv = 100)
sum(svd_list$d[1] * svd_list$d[1]) / sum(svd_list$d * svd_list$d)


timestart = Sys.time()
result = DRMFM(sce2, neighbors, index, q = 15, f = 1, n_iters = 1000, seed=2)
timeend = Sys.time()
timerun = timeend - timestart

timestart = Sys.time() 
result_mis = DRMFM(sce2, neighbors, index, q = 15, f = 1, miss_data = T, n_iters = 800, seed=2)
timeend = Sys.time()
timerun = timeend - timestart

result = result_mis
library(flexclust)
ARI_value = randIndex(table(result$pred_label, colData(sce2)$label))
K_est = length(unique(result$pred_label))
    
saving = c(timerun, K_est, ARI_value)

print(saving)

colData(sce2)$DRMFM = result$pred_label

save(saving, result, file = "DLPFCdata/151669/results/seed_2.RData")
save(saving, result, file = "DLPFCdata/151669/results/result_sce.RData")

load(file = "DLPFCdata/151669/results/seed_2.RData")
iter = 90
randIndex(table(result$MCMCList$group_iter[, iter], colData(sce2)$label))
length(unique(result$MCMCList$group_iter[, iter]))

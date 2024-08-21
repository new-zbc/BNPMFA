library(Rcpp)
library(SingleCellExperiment)
library(RcppArmadillo)
library(flexclust)

load(file = "DLPFCdata/151510/data/151510_counts.RData")

source("R/main.R")
source("R/utils.R")

Adj = find_neighbors(sce1, "Visium", "lattice")
neighbors = find_neighbor_index(Adj, "Visium")

sce2 = spatialPreprocess(sce1, n.PCs = 15, n.HVGs = 2000)
index = which(rowData(sce2)$is.HVG)

timestart = Sys.time()
result = DRMFM(sce2, neighbors, index, q = 15, f = 1, n_iters = 1000, seed=4)
timeend = Sys.time()
timerun = timeend - timestart


library(flexclust)
ARI_value = randIndex(table(result$pred_label, colData(sce2)$label))
K_est = length(unique(result$pred_label))
    
saving = c(timerun, K_est, ARI_value)

print(saving)

colData(sce2)$DRMFM = result$pred_label

save(saving, result, file = "DLPFCdata/151510/results/seed_4.RData")
pred_label = result$pred_label
save(pred_label, file = "DLPFCdata/151510/others/DRMFM.RData")
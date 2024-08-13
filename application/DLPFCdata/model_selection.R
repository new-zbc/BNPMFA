library(Rcpp)
library(SingleCellExperiment)
library(RcppArmadillo)
library(flexclust)

ARGV = commandArgs(trailingOnly = TRUE)
sampleID = as.numeric(ARGV[1])
seed = as.numeric(ARGV[2])
q = as.numeric(ARGV[3])

if(!dir.exists(paste0("DLPFCdata/", sampleID, "/model_selection"))){
    dir.create(paste0("DLPFCdata/", sampleID, "/model_selection"))
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
output = matrix(0, length(fs), 6)
for(i in 1:length(fs)){
    result = DRMFM(sce2, neighbors, index, q = q, f = fs[i], model = "MFM", n_iters = 300,  seed=seed)
    library(flexclust)
    ARI_value = randIndex(table(result$pred_label, colData(sce2)$label))

    temp = c(fs[i], ARI_value, length(unique(result$pred_label)), result$dev, result$BIC)

    output[i, ] = temp

    save(result, file = paste0("DLPFCdata/", sampleID, "/model_selection/f_", fs[i], ".RData"))

    print(c(sampleID, temp))
}

output = as.data.frame(output)
colnames(output) = c("f", "ARI", "K", "dev1", "dev2", "BIC")
write.csv(output, file = paste0("DLPFCdata/", sampleID, "/model_selection/res.csv"))
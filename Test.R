sampleID = 151671
f = 1.5

load(paste0("revision/", sampleID, "/", "f_",f, ".RData"))
mcmclist = result$MCMCList
pred_label = result$pred_label

table(pred_label)

index = which(pred_label == 3)

N = length(pred_label)
K = result$K

prob_pred = rep(0, N)
for(i in 1:N){
  prob_pred[i] = mean((mcmclist$group_iter[i, 1000:2000] + 1) == pred_label[i])
}

hist(prob_pred)

prob_pred[index]


load(paste0("revision/", sampleID, "/", sampleID, "_counts.RData"))

sce1$sizeFactor

library(SingleCellExperiment)
sort(colSums(assay(sce1, "counts")))[index]

sort(sce1$sizeFactor)[1:20]

source("revision/R/utils.R")
sce2 = spatialPreprocess(sce1, n.PCs = 15, n.HVGs = 2000)
sort(colSums(assay(sce2[rowData(sce2)$is.HVG, ], "counts")))

X = as.matrix(assay(sce2[rowData(sce2)$is.HVG, ], "logcounts"))
sort(colMeans(X == 0), decreasing = T)[index]

outlier = X[, index]


library(Rcpp)
library(SingleCellExperiment)
library(RcppArmadillo)
library(flexclust)

sampleID = 151507
data_folder = paste0("DLPFCdata/", sampleID)

if(!dir.exists(paste0(data_folder, "/others"))){
  dir.create(paste0(data_folder, "/others"))
}

load(file = paste0(data_folder, "/data/", sampleID, "_counts.RData"))

#######################################################################
#############   SC.MEB    #############################################
#######################################################################
set.seed(0)
K_set = 2:10
library(SC.MEB)
library(SingleCellExperiment)
Adj = find_neighbors2(sce1, platform = "Visium")

data_mat = reducedDim(sce1)[, 1:15]

selection = SC.MEB(data_mat, Adj_sp = Adj, K_set = K_set, parallel = FALSE)

res = selectK(selection, K_set = K_set, criterion = "BIC", c = 1)

est_K = res$best_K_BIC
pred_label = res$best_K_label

library(flexclust)
ARI_value = randIndex(table(pred_label, colData(sce1)$label))
message("SC.MEB\n")
print(c(ARI_value, est_K))

save(pred_label, file = paste0(data_folder, "/others/", "SCMEB.RData"))
SCMEB_label = pred_label


#######################################################################
#############   DRSC    #############################################
#######################################################################
set.seed(0)
library(DR.SC)
library(SingleCellExperiment)
library(SC.MEB)
Adj = find_neighbors2(sce1, platform = "Visium")

data_mat = t(assay(sce1[rowData(sce1)$is.HVG, ], "logcounts"))
reslist = DR.SC_fit(data_mat, q = 15, K= 2:10, Adj_sp = Adj, coreNum = 1)

res = selectModel(reslist, criteria = 'MBIC', pen.const=1)

est_K = res$bestK
pred_label = res$cluster

library(flexclust)
ARI_value = randIndex(table(pred_label, colData(sce1)$label))

message("DRSC\n")
print(c(ARI_value, est_K))

save(pred_label, file = paste0(data_folder, "/others/", "DRSC.RData"))
DRSC_label = pred_label


#######################################################################
#############   BayesSpace    #########################################
#######################################################################
set.seed(0)
library(BayesSpace)
K = length(unique(colData(sce1)$label))
sce_result = BayesSpace::spatialCluster(sce1, q=K, use.dimred = "PCA",
                                            platform = "Visium", init.method = "kmeans", model = "normal",
                                            precision = "equal", nrep = 5000, gamma = 3,
                                            alpha = 1, beta = 0.01, save.chain = FALSE)

ARI_value = randIndex(table(colData(sce_result)$spatial.cluster, colData(sce1)$label))

pred_label = colData(sce_result)$spatial.cluster
message("BayesSpace\n")
print(c(ARI_value))

save(pred_label, file = paste0(data_folder, "/others/", "BayesSpace.RData"))
BayesSpace_label = pred_label


#######################################################################
#############   Louvain    ############################################
#######################################################################
set.seed(0)
library(Seurat)
library(flexclust)
x_object = as.Seurat(sce1, counts = "logcounts", data = "logcounts", project =  "sce_to_seurat")

x_object = Seurat::FindNeighbors(x_object, reduction = "PCA", dim=1:15)

x_object = Seurat::FindClusters(x_object)
ARI_value = randIndex(table(x_object@meta.data$seurat_clusters, colData(sce1)$label))

pred_label = x_object@meta.data$seurat_clusters
est_K = length(unique(x_object@meta.data$seurat_clusters))
message("Louvain\n")
print(c(ARI_value, est_K))

save(pred_label, file = paste0(data_folder, "/others/", "Louvain.RData"))
Louvain_label = as.numeric(pred_label) + 1

############ summary
result_data_frame = data.frame(SCMEB = SCMEB_label, 
                               DRSC = DRSC_label,
                               BayesSpace = BayesSpace_label,
                               Louvain = Louvain_label)
write.table(result_data_frame, file = paste0(data_folder, "/others/", "summary.txt"))
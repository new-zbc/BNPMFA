library(Rcpp)
library(SingleCellExperiment)
library(RcppArmadillo)
library(flexclust)

ARGV = commandArgs(trailingOnly = TRUE)
sampleID = as.numeric(ARGV[1])
#sampleID = 151507
data_folder = paste0("DLPFCdata/", sampleID)

if(!dir.exists(paste0(data_folder, "/others"))){
  dir.create(paste0(data_folder, "/others"))
}

load(file = paste0(data_folder, "/data/", sampleID, "_counts.RData"))

#######################################################################
#############   Louvain    ############################################
#######################################################################
set.seed(1)
library(Seurat)
library(flexclust)


K_true = length(unique(sce1$label))

rs = c(0.2, 0.25, 0.3, 0.35,  0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3)


delta_old = 100
output = NA
for(i in 1:length(rs)){
  print(rs[i])
  x_object = as.Seurat(sce1, counts = "logcounts", data = "logcounts", project =  "sce_to_seurat")
  x_object = Seurat::FindNeighbors(x_object, reduction = "PCA", dim=1:15)
  x_object = Seurat::FindClusters(x_object, resolution = rs[i])
   ARI_value = randIndex(table(x_object@meta.data$seurat_clusters, colData(sce1)$label))
   pred_label = x_object@meta.data$seurat_clusters
   est_K = length(unique(x_object@meta.data$seurat_clusters))
   delta_new = abs(est_K - K_true)
   if(delta_new == 0){
    output = pred_label
    break}
   else{
    if(delta_new < delta_old){
        output = pred_label
        delta_old = delta_new
    }
   }
}

pred_label = as.numeric(output)
est_K = length(unique(pred_label))
message("Louvain\n")
print(c(ARI_value, est_K))

save(pred_label, file = paste0(data_folder, "/others/", "Louvain.RData"))
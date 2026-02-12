library(foreach)
library(aricode)
library(SingleCellExperiment)
ARGV = commandArgs(trailingOnly = TRUE)
folder = ARGV[1]

method = "Louvain"

data_folder = paste0("simulation3/scenario2_1/", folder)

n.data_set = 50

if(!dir.exists(paste0(data_folder, "/summary"))){
  dir.create(paste0(data_folder, "/summary"))
}


if(!dir.exists(paste0(data_folder, "/", method))){
  dir.create(paste0(data_folder, "/", method))
}



mydata1 <- foreach(i=1:n.data_set, .combine = "rbind") %do%
{
  filename = paste0(data_folder, "/data/", i, ".RData")
  load(file = filename)
#   sce = sce1
#   sce$row = round(sce$row)
#   sce$col = round(sce$col)
#   library(scran)
#   sce <- scater::logNormCounts(sce)
#   dec <- scran::modelGeneVar(sce, assay.type = "logcounts")
#   top <- scran::getTopHVGs(dec, n = 2000)
#   sce <- scater::runPCA(sce, subset_row = top, ncomponents=15, 
#                         exprs_values="logcounts")
  
  K = length(unique(colData(sce)$label))
  
  
  set.seed(0)
    
  library(Seurat)
  library(flexclust)
  x_object = as.Seurat(sce, counts = "logcounts", data = "logcounts", project =  "sce_to_seurat")

  x_object = Seurat::FindNeighbors(x_object, reduction = "PCA", dim=1:15)
    
    timestart = proc.time()
    for(resolution in c(0.05, 0.1, 0.15, 0.2, 0.3, 0.4 ,0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5)){
      x_object = Seurat::FindClusters(x_object, resolution = resolution, algorithm = 1, random.seed=0)
      
      pred_label = x_object@meta.data$seurat_clusters
      
      
      K_est = length(unique(pred_label))
      
      if(K_est == K){break}
    }
    timeend = proc.time()
    time = timeend - timestart
    
    library(flexclust)
    ARI_value = randIndex(table(pred_label, colData(sce)$label))
    AMI_value = AMI(pred_label, colData(sce)$label)
    NMI_value = NMI(pred_label, colData(sce)$label)
      
    file_name = paste0(data_folder,"/", method, "/", i, ".RData")
      
    out = c(ARI_value, AMI_value, NMI_value, time[3])

    save(pred_label, out, file = file_name)

    print(out)
}

file_name = paste0(data_folder, "/summary/", method, ".txt")
write.table(mydata1, file = file_name)
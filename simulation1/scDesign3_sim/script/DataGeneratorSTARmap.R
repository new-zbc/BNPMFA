library(SingleCellExperiment)
source("R/utils.R")
filename = paste0("STARmap/data/BZ5.RData")
load(file = filename)
  
sce = spatialPreprocess(sce, n.PCs = 5, n.HVGs = dim(sce)[1])
  
  # Calculate Voronoi Tesselation and tiles
  ## Examples of finding neighbors using Voronoi tessellation (load any data)
temp <- data.frame(id = 1:dim(sce)[2], x = colData(sce)$col, y = colData(sce)$row)
tt = voronoi_adjacency(data = temp, id~x+y, scale=1, PLOT=FALSE)
Adj = tt$G   ## a n-by-n Adjacency Matrix, 1 if two samples are neighbors; 0 otherwise
neighbors = tt$P

sce_original <- scater::logNormCounts(sce)
#dec <- scran::modelGeneVar(sce_original, assay.type = "logcounts")
#top <- scran::getTopHVGs(dec, n = 1000)
# # sce <- scater::runPCA(sce, subset_row = top, ncomponents=15, 
# #                       exprs_values="logcounts")

sce_original = sce_original[, !is.na(sce_original$label)]

library(scDesign3)

family = "STARmap_nb"

if(!dir.exists(paste0("scDesign3/", family))){
  dir.create(paste0("scDesign3/", family))
}

data_folder = paste0("scDesign3/", family)
dir.create(paste0(data_folder, "/data"))

for(data_index in 31:50){
  set.seed(data_index)
  example <- scdesign3(sce = sce_original, 
                       assay_use = "counts",
                       celltype = "label",
                       pseudotime = NULL,
                       other_covariates = NULL,
                       spatial = c("row", "col"),
                       mu_formula = "s(row, col, bs = 'gp', k= 100)",
                       n_cores = 5,
                       sigma_formula = "1",
                       family_use = "nb",
                       corr_formula = "1",
                       copula = "gaussian")
  
  sce <- SingleCellExperiment(list(counts = example$new_count), 
                              colData = DataFrame(row = example$new_covariate$row, 
                                                  col = example$new_covariate$col, 
                                                  label = example$new_covariate$label))
  
  sce <- scater::logNormCounts(sce)
  sce <- scater::runPCA(sce, subset_row=1:dim(sce)[1], ncomponents=10, 
                        exprs_values="logcounts")
  
  Adj = find_neighbors(sce, "ST", "lattice")
  neighbors = find_neighbor_index(Adj, "ST")
  
  
  file_name = paste0(data_folder, "/", "data/", data_index, ".RData")
  save(sce, Adj, neighbors,  file = file_name)
}



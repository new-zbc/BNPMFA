#data_file = "BZ5"
ARGV = commandArgs(trailingOnly = TRUE)
data_file = ARGV[1]
  
  if(!dir.exists(paste0("STARmap/model_selection/", data_file))){
    dir.create(paste0("STARmap/model_selection/", data_file))
  }
  ######################################################
  ### Proposed Method
  ######################################################
  source("R/utils.R")
  filename = paste0("STARmap/data/", data_file, ".RData")
  load(file = filename)
  
  sce = spatialPreprocess(sce, n.PCs = 5, n.HVGs = dim(sce)[1])
  
  # Calculate Voronoi Tesselation and tiles
  ## Examples of finding neighbors using Voronoi tessellation (load any data)
  temp <- data.frame(id = 1:dim(sce)[2], x = colData(sce)$col, y = colData(sce)$row)
  
  tt = voronoi_adjacency(data = temp, id~x+y, scale=1, PLOT=FALSE)
  Adj = tt$G   ## a n-by-n Adjacency Matrix, 1 if two samples are neighbors; 0 otherwise
  neighbors = tt$P   ## a n-by-m matrix, m: maximum # of neighbors changing from data to data, 
  
  index = which(rowData(sce)$is.HVG)
  source("R/main.R")
  fs = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)
  output = matrix(0, length(fs), 6)
  for(i in 1:length(fs)){
    timestart = Sys.time()
  result = DRMFM(sce, neighbors, features = 1:dim(sce)[1], q = 5, f = fs[i], model = "MFM", n_iters = 200, seed = 2)
  timeend = Sys.time()
  timerun = timeend - timestart
  
  library(flexclust)
  ARI_value = randIndex(table(result$pred_label, colData(sce)$label))
  pred_label = result$pred_label
  K_est = length(unique(result$pred_label))

  temp = c(fs[i], ARI_value, K_est, result$dev, result$BIC)

  output[i, ] = temp
  
 save(result, file = paste0("STARmap/model_selection/", data_file, "/f_", fs[i], ".RData"))
  }

output = as.data.frame(output)
colnames(output) = c("f", "ARI", "K", "dev1", "dev2", "BIC")
write.csv(output, file = paste0("STARmap/model_selection/", data_file, "/res.csv"))
print(paste0(data_file " Done"))


  
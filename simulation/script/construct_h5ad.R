library(SingleCellExperiment)
library(Seurat)
#library(SeuratData)
library(SeuratDisk)


data_name = "simulation1/scenario3_2"

if(!dir.exists(paste0(data_name,"/SpaGCN"))){
  dir.create(paste0(data_name,"/SpaGCN"))
}

if(!dir.exists(paste0(data_name,"/SpaGCN/data"))){
  dir.create(paste0(data_name,"/SpaGCN/data"))
}


for(i in 1:50){
  load(paste0( data_name,"/data/", i, ".RData"))
  
  
  
  seu <- Seurat::as.Seurat(sce, counts = "logcounts", assay = NULL, project = "SingleCellExperiment")
  
  
  SaveH5Seurat(seu, filename = paste0(data_name, "/SpaGCN/data/", i, ".h5Seurat"), overwrite = T)
  Convert(paste0(data_name, "/SpaGCN/data/", i, ".h5Seurat"), dest = "h5ad")
  unlink(paste0(data_name, "/SpaGCN/data/", i, ".h5Seurat"))
}



# snp_mart = useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")

# snp_ids = c("rs16828074", "rs17232800")
# snp_attributes = c("refsnp_id", "chr_name", "chrom_start")

# snp_locations = getBM(attributes=snp_attributes, filters="snp_filter", 
#                       values=snp_ids, mart=snp_mart)

# snp_locations


# ensembl <- useMart("ensembl")
# ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
# getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'strand'),
#       filters=c('hgnc_symbol'),
#       values=list('A1BG-AS1'),
#       mart=ensembl)
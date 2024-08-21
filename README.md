# BNPMFA

### Overview

 BNPMFA is a generalized Bayesian nonparametric clustering framework for spatial transcriptomic data. We propose a novel Bayesian nonparametric mixture of factor analysis (BNPMFA) model, which incorporates a Markov random field-constrained Gibbs-type prior for partitioning high-dimensional spatial omics data. This new prior effectively integrates
the spatial constraints inherent in SRT data while simultaneously inferring cluster membership and determining the optimal number of spatial domains.

### User manual

The required R packages:

* Rcpp
* SingleCellExperiment
* RcppArmadillo
* flexclust
* scater

The data structure used for BNPMFA is a `SingleCellExperiment` object containing assay "counts" and col data "row" and "col". "row" and "col" are the coordinates of each spot. Here we use DLPFC 151672 as an example. 

```R
# Load data
library(SingleCellExperiment)
source("R/utils.R")
load("application/DLPFCdata/151672/data/151672_counts.RData")
sce1
# class: SingleCellExperiment 
# dim: 33538 3888 
# metadata(0):
# assays(2): counts logcounts
# rownames(33538): ENSG00000243485 ENSG00000237613 ... ENSG00000277475
# ENSG00000268674
# rowData names(2): gene_name is.HVG
# colnames(3888): AAACAAGTATCTCCCA AAACACCAATAACTGC ... TTGTTTCCATACAACT
#   TTGTTTGTGTAAATTC
# colData names(7): row col ... sizeFactor spaGCN
# reducedDimNames(1): PCA
# mainExpName: NULL
# altExpNames(0):
```

Data preprocess

```R
# log-normalization for counts data
sce1 = scater::logNormCounts(sce1)
# find highly variable genes
sce2 = spatialPreprocess(sce1, n.HVGs = 2000)
index = which(rowData(sce2)$is.HVG)
```

Construct neighborhood adjacent matrix

```R
Adj = find_neighbors(sce1, "Visium", "lattice")
# Adj is the adjacent matrix
neighbors = find_neighbor_index(Adj, "Visium")
# neighbors is a n-by-num_neighbors matrix where each row contains the index of neighbors for each spots. It is transformed from adjacent matrix
```

Run BNPMFA

```R
source("R/main.R")
result = DRMFM(sce=sce2, neighbors=neighbors, features=index, 
               q = 15, f = 1, K_init = 10, model = "MFM" n_iters = 1000, seed=1)
```

These are the definition of parameters:

`sce` is the `SingleCellExperiment` object with "logcounts" assay.

`neighbors` is the neighborhood information for each spot.

`features` is the index for highly variable genes.

`q` is the latent dimension of factor analysis model.

`f` is the Markov random field hyperparameters.

`K_init` is the initial values for the number of clusters.

`model` is the choice for different partition models. We only have three choices where "DP" represents Dirichlet process, "PY" represents Pitman-Yor process, and "MFM" represents mixture-of-finite-mixtures model. 

`n_iters` is the total number of MCMC iterations. 

`seed` is the seed to reproduce result. 

```R
# cluster result
pred_label = result$pred_label
library(flexclust)
# compute ARI
ARI_value = randIndex(table(pred_label, colData(sce2)$label))
```


library(Rcpp)
library(RcppArmadillo)
library(SingleCellExperiment)
library(irlba)

sourceCpp("R/MFMmcmc.cpp")
sourceCpp("R/DPmcmc.cpp")
sourceCpp("R/PYmcmc.cpp")

DRMFM <- function(sce, neighbors, features, q=5, f=1, K_init = 10, n_iters = 100, model = "MFM", seed = 2){
  
  tau = 1
  tau_w = 1
  set.seed(seed)
  data_set = as.matrix(assay(sce, "logcounts"))
  data_set = data_set[features, ]
  N = dim(data_set)[2]
  p = dim(data_set)[1]
  
  K_start = K_init
  
  svd_list = irlba::irlba(t(data_set), nv = q)
  
  W_t = svd_list$v 
  
  Y_t = svd_list$u %*% diag(svd_list$d[1:q])
  #Y_t = Y
  
  group_t = kmeans(Y_t, centers = K_start)$cluster
  
  var_inv_t = rep(1, p)
  
  mu_t = matrix(0, q, K_start)
  
  cov_t = diag(1, q, q)
  cov_inv_t = diag(1, q, q)
  
  if(model == "MFM"){
    result = MFMmcmc(group_t-1, W_t, var_inv_t, t(Y_t), mu_t, cov_t, cov_inv_t, 
                     data_set, G = neighbors, f = f, tau=tau, tau_w=tau_w, max_iters = n_iters)
  }else if(model == "DP"){
    result = DPmcmc(group_t-1, W_t, var_inv_t, t(Y_t), mu_t, cov_t, cov_inv_t, 
                     data_set, G = neighbors, f = f, tau=tau, tau_w=tau_w, max_iters = n_iters)
  }else if(model == "PY"){
    result = PYmcmc(group_t-1, W_t, var_inv_t, t(Y_t), mu_t, cov_t, cov_inv_t, 
                    data_set, G = neighbors, f = f, tau=tau, tau_w=tau_w, max_iters = n_iters)
  }else {stop("model == FMF or DP or PY")}
  
  
  max_class = max(result$K_iter)
  prob_pred = matrix(0, N, max_class)
  for(i in 1:N){
    for(j in 1:max_class){
      prob_pred[i, j] = mean((result$group_iter[i, ] + 1) == j)
    }
  }
  
  pred_label = as.vector(apply(prob_pred, 1,  function(x) which(x == max(x))[1]))
  K = length(unique(pred_label))
  output = list()
  output$pred_label = pred_label
  output$K = K
  output$MCMCList = result
  
  dev1 = result$dev_iter[1, n_iters]
  dev2 = result$dev_iter[2, n_iters]
 
  group_est = result$group_iter[, n_iters] + 1
  n.par = length(unique(group_est)) * q + q*(q-1)/2 + p + p *q
  
  BIC = -2*(dev1 + dev2) + n.par * log(N)
  output$BIC = BIC
  output$dev = c(dev1, dev2)
  
  return(output)
  
}





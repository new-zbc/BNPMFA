#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <Rmath.h>
#include <iostream>

using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends("RcppArmadillo")]]


// calculate the density of Poisson distribution
double dpois(int x, double lambda, bool log_value = true){
  if(lambda <= 0){
    throw std::logic_error("lambda should be positive");
  }
  double r;
  if(x == 0){
    r = - lambda;
  }else{
    vec s1 = arma::linspace(1, x, x);
    r = -sum(log(s1)) + x * log(lambda) - lambda;
  }
  if(log_value){return r;}
  else{return exp(r);}
}



double dmvnorm(const vec x, const vec mu, const mat sigma, bool log_value = true){
  double result = 0;
  result += - as_scalar((x - mu).t() * inv(sigma) * (x - mu))/2;
  result += - (log(2*datum::pi) * x.n_elem + log(det(sigma))) / 2;
  if(log_value){return result;}
  else{return exp(result);}
}



void update_y(mat &Y_t, const mat &W_t, const vec &group_t, const vec &var_inv_t, 
              const mat &mu_t, const mat &cov_inv_t, const mat X){
  
  mat Lambda = diagmat(var_inv_t);
  mat M1 = W_t.t() * Lambda;
  mat sigma_mat = inv(M1 * W_t + cov_inv_t);
  int N = group_t.n_elem;
  for(int i=0; i < N; i++){
    int label = group_t(i);
    mat sigma_i = sigma_mat;
    vec mu_i = M1 * X.col(i) + cov_inv_t * mu_t.col(label);
    vec y_i = mvnrnd(sigma_i * mu_i, sigma_i, 1);
    Y_t.col(i) = y_i;
  }
}


void update_W(mat &W_t, const mat &Y_t, const vec &var_inv_t, const mat &X, const double tau){
  int q = Y_t.n_rows;
  int p = W_t.n_rows;
  mat sigma = inv(Y_t * Y_t.t() + tau * arma::eye(q, q));
  for(int j = 0; j < p; j++){
    colvec x_j = X.row(j).t();
    vec mu_j = Y_t * x_j;
    vec w_j = mvnrnd(sigma * mu_j, sigma, 1);
    W_t.row(j) = w_j.t();
  }
}



void update_sigma(vec &var_inv_t, const mat &Y_t, const mat &W_t, const mat &X,  const double &tau,
                  const double &delta1, const double &delta2){
  int p = var_inv_t.n_elem;
  int q = Y_t.n_rows;
  int N = X.n_cols;
  for(int j = 0; j < p; j++){
    rowvec s_j = X.row(j) - W_t.row(j) * Y_t;
    double residual1 = as_scalar(s_j * s_j.t() / 2.0 );
    double residual2 = tau * sum(W_t.row(j) % W_t.row(j))/2.0;
    double sample = randg(1, distr_param(delta1 + N/2.0 + q/2.0, 1.0 /(delta2 + residual1 + residual2)))(0);
    var_inv_t(j) = sample;
  }
}


// function to update mean of each group 
void update_mu(mat &mu_t, const vec &group_t, const mat &cov_inv_t, const mat &cov_t,
               const mat &Y_t, const double &tau){
  int K = mu_t.n_cols;
  for(int k = 0; k < K; k++){
    uvec index = find(group_t == k);
    int n_k = index.n_elem;
    vec Y_k_sum = sum(Y_t.cols(index), 1);
    
    vec mu= Y_k_sum / (tau + n_k);
    
    vec mu_k = mvnrnd(mu, cov_t / (n_k +tau), 1);
    mu_t.col(k) = mu_k;
  }
}


// function to update covariance of each group
void update_cov(mat &cov_inv_t, mat &cov_t, const vec &group_t, const mat &Y_t, const mat &mu_t,
                const double &tau){
  int N = group_t.n_elem;
  int K = mu_t.n_cols;
  int q = cov_t.n_cols;
  mat S = zeros<mat>(q, q);
  for(int k = 0; k < K; k++){
    uvec index = find(group_t == k);
    mat Y_k = Y_t.cols(index);
    vec mu_k = mu_t.col(k);
    Y_k = (Y_k.each_col() - mu_k);
    S += Y_k * Y_k.t();
    S += tau * mu_k * mu_k.t();
  }
  
  int eta = N + K;
  cov_t = arma::iwishrnd(S, eta);
  //cov_t = diagmat(S);
  cov_inv_t = inv(cov_t);
}




//function to calculate MRF energy function
double MRF_energy(const int &k, const int &l, const vec &group_t, const umat &G,
                  const double f, const bool log_value = true){
  
  uvec G_l = G.row(l).t();
  uvec neighbor = G_l.elem(find(G_l > 0)) - 1;
  
  uvec neighbor_id = find(group_t(neighbor) == k);
  double result = f * neighbor_id.n_elem;
  
  if(log_value){return result;}
  else{ return exp(result); }
}


//function to calculate probability of existing clusters
double prob_existing(const int k, const int l, const vec group_t, const mat mu_t, const mat cov_t,
                     const vec y, const umat G, const double f, const double GAMMA){
  double result = 0;
  uvec index = find(group_t == k);
  int n_k;
  if(group_t(l) == k){
    n_k = index.n_elem - 1;
  }else{ n_k = index.n_elem;}
  
  result += log(n_k + 0.0001) + MRF_energy(k, l, group_t, G, f);
  
  result += dmvnorm(y, mu_t.col(k), cov_t);
  
  return result;
}

// function to calculate probability of a new cluster
double prob_new(const int K, const vec y, const mat cov_t, const double GAMMA, 
                const double tau){
  int q = y.n_elem;
  vec mu = zeros<vec>(q);
  
  double h = dmvnorm(y, mu, cov_t*(1 + 1.0/tau));
  
  double result = max((h + log(GAMMA) ), -999999.0);
  return result;
}



void update_group(vec &group_t, mat &mu_t, mat &cov_t, const mat &Y_t, 
                  const int &Kmax, const umat &G, const double &f, const double &GAMMA, 
                  const double &tau, const double &temperature){
  
  int N = Y_t.n_cols;
  int q = Y_t.n_rows;
  vec clusters_name = unique(group_t);
  int K = clusters_name.n_elem;
  
  for(int i = 0; i < N; i++){
    int label_old = group_t(i);
    int label_new;
    uvec group_i = find(group_t == label_old);
    int c_size = group_i.n_elem;
    
    
    if(c_size > 1){
      vec prob_i = zeros<vec>(K+1);
      for(int k = 0; k < K; k++){
        prob_i(k) = prob_existing(k, i, group_t, mu_t, cov_t, Y_t.col(i), G, f, GAMMA);
      }
      
      prob_i(K) = prob_new(K, Y_t.col(i), cov_t, GAMMA, tau);
      
      
      prob_i = (prob_i - max(prob_i))/temperature;
      prob_i = exp(prob_i) / sum(exp(prob_i));
      
      label_new = as_scalar(Rcpp::RcppArmadillo::sample(regspace(0, K), 1, false, prob_i));
      
      if(label_new >= K && K < Kmax){
        
        //initial parameters for new cluster
        mat mu_t_new = zeros<mat>(q, (K+1));
        mu_t_new.head_cols(K) = mu_t;
        mu_t_new.col(K) = mvnrnd(Y_t.col(i)/(1 + tau), cov_t/(1 + tau), 1);
        mu_t.swap(mu_t_new);
        
        
        group_t(i) = label_new;
        clusters_name = unique(group_t);
        K = clusters_name.n_elem;
      }
      else if(label_new >= K && K >= Kmax){
        
        group_t(i) = label_old;
        clusters_name = unique(group_t);
        K = clusters_name.n_elem;
      }
      else{
        group_t(i) = label_new;
        clusters_name = unique(group_t);
        K = clusters_name.n_elem;
      }
      
    }
    
    else{
      
      //The cluster is a singleton, only K choices
      vec prob_i = vec(K+1, fill::zeros);
      for(int k = 0; k < K; k++){
        prob_i(k) = prob_existing(k, i, group_t, mu_t, cov_t, Y_t.col(i), G, f, GAMMA);
      }
      
      prob_i(K) = prob_new(K, Y_t.col(i), cov_t, GAMMA, tau);
      
      prob_i = (prob_i - max(prob_i))/temperature;
      prob_i = exp(prob_i) / sum(exp(prob_i));
      
      label_new = as_scalar(Rcpp::RcppArmadillo::sample(regspace(0, K), 1, false, prob_i));
      
      if(label_new == label_old || label_new == K){
        group_t(i) = label_old;
        clusters_name = unique(group_t);
        K = clusters_name.n_elem;
      }
      else{
        uvec index = find(clusters_name != label_old);
        group_t(i) = label_new;
        for( int id = 0; id < N; id++){
          if(group_t(id) > label_old){
            group_t(id) += -1;
          }
        }
        clusters_name = unique(group_t);
        K = clusters_name.n_elem;
        
        if( K >= 1){
          
          mu_t = mu_t.cols(index);
          
        }
      }
    }
  }
}


vec Dev(const vec &group_t, const mat &W_t, const vec &var_inv_t, const mat &Y_t, 
        const mat &mu_t, const mat &cov_t, const mat &cov_inv_t, const mat &X){
  vec output = zeros<vec>(2);
  int N = X.n_cols;
  double dev1 = 0;
  double dev2 = 0;
  for(int i=0; i < N; i++){
    int label = group_t(i);
    vec mu = W_t * Y_t.col(i);
    vec loglike = arma::log_normpdf(X.col(i), mu, 1/sqrt(var_inv_t));
    dev1 += sum(loglike);
    dev2 += dmvnorm(Y_t.col(i), mu_t.col(label), cov_t);
  }
  output(0) = dev1;
  output(1) = dev2;
  return(output);
}



// [[Rcpp::export]]
Rcpp::List DPmcmc(vec &group_t, mat &W_t, vec &var_inv_t, mat &Y_t, mat &mu_t, mat &cov_t, mat &cov_inv_t, 
                   const mat &X, const umat &G, const double &f, const double &tau, const double &tau_w, const double GAMMA=1,
                   const int max_iters = 100, const int seed = 12569){
  int N = X.n_cols;
  int p = X.n_rows;
  int q = Y_t.n_rows;
  int Kmax = 30;
  
  arma_rng::set_seed(seed);
  
  //store the parameters
  vec K_iter = zeros<vec>(max_iters);
  mat group_iter = zeros<mat>(N, max_iters);
  cube W_iter = zeros<cube>(p, q, max_iters);
  mat var_iter = zeros<mat>(p, max_iters);
  cube mu_iter = zeros<cube>(q, Kmax, max_iters);
  cube cov_iter= zeros<cube>(q, q, max_iters);
  //store deviance
  mat dev_iter = zeros<mat>(2, max_iters);
  
  
  //store the initialized pars
  vec clusters = unique(group_t);
  int K = clusters.n_elem;
  K_iter(0) = K;
  group_iter.col(0) = group_t;
  W_iter.slice(0) = W_t;
  var_iter.col(0) = 1.0 / var_inv_t;
  mu_iter.slice(0).head_cols(K) = mu_t;
  cov_iter.slice(0) = cov_t;
  dev_iter.col(0) = Dev(group_t, W_t, var_inv_t, Y_t, mu_t, cov_t, cov_inv_t, X);
  
  
  for(int t = 1; t < max_iters; t++){
    double temperature = pow(0.999, t);
    
    update_W(W_t, Y_t, var_inv_t, X, tau_w);
    
    update_sigma(var_inv_t, Y_t, W_t, X, tau_w, 0.01, 0.01);
    
    update_mu(mu_t, group_t, cov_inv_t, cov_t, Y_t, tau);
    
    update_cov(cov_inv_t, cov_t, group_t, Y_t, mu_t, tau);
    
    
    update_group(group_t, mu_t, cov_t, Y_t, Kmax, G, f, GAMMA, tau, temperature);
    
    
    update_y(Y_t, W_t, group_t, var_inv_t, mu_t, cov_inv_t, X);
    
    
    vec clusters = unique(group_t);
    int K = clusters.n_elem;
    K_iter(t) = K;
    group_iter.col(t) = group_t;
    W_iter.slice(t) = W_t;
    var_iter.col(t) = 1.0 / var_inv_t;
    mu_iter.slice(t).head_cols(K) = mu_t;
    cov_iter.slice(t) = cov_t;
    dev_iter.col(t) = Dev(group_t, W_t, var_inv_t, Y_t, mu_t, cov_t, cov_inv_t, X);
  }
  
  return List::create(Named("K_iter") = K_iter,
                      Named("group_iter") = group_iter, 
                      Named("W_iter") = W_iter, 
                      Named("var_iter") = var_iter,
                      Named("mu_iter") = mu_iter,
                      Named("cov_iter") = cov_iter,
                      Named("dev_iter") = dev_iter);
}








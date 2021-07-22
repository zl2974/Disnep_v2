#### perfect transition matrix

perfect = function(signals,gene_int,n_signal,n_noise,corr){
  
  source("R/diffus_vec.R")
  
  a = matrix(corr,nrow = n_signal,ncol = n_signal)
  
  b = matrix(1e-3,nrow=n_noise,ncol = n_noise)
  
  adjacency = as.matrix(Matrix::bdiag(a,b))
  
  diag(adjacency) = 0
  
  result = as_tibble(diffus_vec(signals,adjacency,type = "statistics",beta = 0.75,iter = 100))
  
  return(result)
}
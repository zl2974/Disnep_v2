#### perfect transition matrix

perfect = function(signals,gene_int,n_signal,n_noise,corr,beta){
  
  source("R/diffus_vec.R")
  source("R/diffus_matrix.R")
  source("R/post_process.R")
  
  adjacency = matrix(1e-3,nrow=(n_signal+n_noise),ncol = (n_signal+n_noise))
  
  adjacency[1:n_signal,] = corr/20
  
  adjacency[,1:n_signal] = corr/20
  
  adjacency[1:n_signal,1:n_signal] = corr/5
  
  adjacency[(n_signal/2):n_signal,(n_signal/2):n_signal] = corr/10
  
  adjacency[1:(n_signal/2),1:(n_signal/2)] = corr/2
  
  diag(adjacency) = 0
  
  rownames(adjacency) = colnames(gene_int)
  colnames(adjacency) = colnames(gene_int)
  
  adjacency = diffus_matrix(gene_int,adjacency,alpha = beta)
  
  adjacency = post_process(adjacency)
  
  result = as_tibble(diffus_vec(signals,adjacency,type = "statistics",beta = beta,iter = 100))
  
  return(result)
}
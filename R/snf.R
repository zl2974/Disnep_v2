### snf

snf = function(signals,
               gene_int,
               gene_expr,
               k = 200,
               iter = 10,
               beta,
               n_signal,
               n_noise,
               corr,
               ...) {
  
  # normalize gene interaction
  sortedColumn = sort(colnames(gene_int))
  
  gene_int = gene_int[sortedColumn,sortedColumn]
  
  gene_int = gene_int / 999
  
  
  # generate adjacency matrix
  
  gene_expr = gene_expr[,sortedColumn]
  
  #adjacency = abs(cor(gene_expr))
  
  adjacency = matrix(1e-3,nrow=(n_signal+n_noise),ncol = (n_signal+n_noise))
  
  adjacency[1:n_signal,] = corr/20
  
  adjacency[,1:n_signal] = corr/20
  
  adjacency[1:n_signal,1:n_signal] = corr/5
  
  adjacency[(n_signal/2):n_signal,(n_signal/2):n_signal] = corr/10
  
  adjacency[1:(n_signal/2),1:(n_signal/2)] = corr/2
  
  diag(adjacency) = 0
  
  rownames(adjacency) = colnames(gene_int)
  colnames(adjacency) = colnames(gene_int)
  
  # SNF fuse
  snet = SNFtool::SNF(list(gene_int, adjacency), K = k, t = iter)
  
  source("R/post_process.R")
  snet_post = post_process(snet, percent = 0.9)
  
  snet_post = snet_post[sortedColumn,sortedColumn]
  
  
  # Random Walk With Restart
  source("R/diffus_vec.R")
  
  signals = arrange(signals,gene)
  
  score = diffus_vec(signals = signals,snet = snet_post,type = "statistics",iter = 100,beta=beta)
  
  return(score)
  
}
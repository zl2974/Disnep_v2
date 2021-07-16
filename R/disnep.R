##### DISNEP

disnep = function(signals,
                  gene_int,
                  gene_expr,
                  type,
                  alpha = 0.75,
                  beta = 0.75,
                  iter = 100,
                  difference = 1e-6) {
  source("R/diffus_matrix.R")
  source("R/diffus_vec.R")
  source("R/post_process.R")
  
  ### make sure that the gene row col doesn't change
  
  sorted_colname = sort(colnames(gene_expr))
  
  adjacency = abs(cor(gene_expr, method = "pearson"))
  
  adjacency = adjacency[sorted_colname,sorted_colname]
  
  diag(adjacency) = 0
  
  ### sort and build enhanced network
  
  gene_int = gene_int[sorted_colname,sorted_colname]
  
  se = diffus_matrix(s0 = gene_int, adjacency = adjacency, alpha, iter)
  
  se = se[sorted_colname,sorted_colname]
  
  snet_post = post_process(se, percent = 0.9)
  
  ### prioritizing
  
  snet_post = snet_post[sorted_colname,sorted_colname]
  
  signals = dplyr::arrange(signals,gene)
  
  score = as_tibble(diffus_vec(signals, snet_post, type, iter))
  
  return(score)
  
}

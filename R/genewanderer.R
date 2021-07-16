#### genewanderer

genewanderer = function(signals,
                        snet,
                        type,
                        beta = 0.75,
                        iter = 10,
                        difference = 1e-6) {
  source("R/diffus_vec.R")
  
  score = as_tibble(diffus_vec(signals, snet, type, beta, iter, difference))
  
  return(score)
  
}
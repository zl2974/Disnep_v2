#### genewanderer

genewanderer = function(signals,
                        snet,
                        type,
                        beta = 0.75,
                        iter = 10,
                        difference = 1e-6) {
  source("R/diffus_vec.R")
  source("R/post_process.R")
  
  snet = post_process(se = snet,percent = .9)
  
  score = as_tibble(
    diffus_vec(
      signals = signals,
      snet = snet,
      type = type,
      beta = beta,
      iter = iter,
      difference = difference
    )
  )
  
  return(score)
  
}
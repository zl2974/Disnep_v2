### Random Walk with restart

RWR = function(signals,
               trans,
               beta,
               normalize_type = 1,
               .max_iter = 100,
               .tol = 1e-6
){
  ### column normalize transition matrix
  
  diag(trans) = 0
  
  source("R/normalize.R")
  trans = normalize(trans,type = normalize_type)
  
  ### check signals tibble
  
  signals = arrange(signals,gene)
  
  p0 = signals$score
  
  current = signals$score
  
  previous = current+Inf
  
  ### random walk with restart
  
  .iter = 1
  
  while(abs(norm(current-previous))>.tol & .iter < .max_iter){
    
    previous = current
    
    current = beta * trans%*%previous + (1-beta)*p0
    
  }
  
  ### collect result
  
  result = tibble(gene = signals$gene,
                  score = current) %>% 
    arrange(score)
  
  return(result)
  
  
}
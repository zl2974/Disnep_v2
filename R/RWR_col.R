RWR_col = function(.data,
                   .transit,
                   .beta = 0.75,
                   .tol = 1e-6,
                   .max_iter = 100,
                   ...) {
  
  source("R/normalize.R")
  diag(.transit) = 0
  .transit = normalize(.transit,1)
  
  .data = as.matrix(.data)
  
  current = as.matrix(.data)
  
  previous = current - (-Inf)
  
  t0 = Sys.time()
  
  .iter = 0
  
  while (abs(norm(previous - current)) > .tol & .iter<.max_iter) {
    previous = current
    
    .iter = .iter + 1
    
    current = .beta * .transit %*% previous + (1-.beta) * .data
    
  }
  
  return(current)
}
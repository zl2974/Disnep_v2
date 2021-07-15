### simulate case

simulate_case = function(n_signal,
                         n_noise,
                         strong_signal_mu,
                         median_signal_mu,
                         noise_mu,
                         strong_corr,
                         median_corr,
                         noise_corr,
                         sigma,
                         n_case,
                         ...){
    
    mu_strong_signal = rnorm(n_signal/2,strong_signal_mu,1)
    
    mu_median_signal = rnorm(n_signal/2,median_signal_mu,1)
    
    mu_noise = rnorm(n_noise,noise_mu,1)
    
    Mu = c(mu_strong_signal,mu_median_signal,mu_noise)
    
    sigma_strong = matrix(runif((n_signal/2)^2,min = 0,max = strong_corr),n_signal/2)
    
    sigma_median = matrix(runif((n_signal/2)^2,min = 0,max = median_corr),n_signal/2)
    
    sigma_noise = matrix(runif(n_noise^2,min = 0,max = noise_corr),n_noise)
    
    Sigma = as.matrix(Matrix::bdiag(sigma_strong,sigma_median,sigma_noise))
    
    diag(Sigma) = 1
    
    case_data = MASS::mvrnorm(n_case,mu = Mu,Sigma = sigma^2*Sigma)
    
    colnames(case_data) = stringr::str_c("V", 1001:(1000 + n_signal + n_noise))
    
    return(case_data)
}
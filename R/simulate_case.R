### simulate case

simulate_case = function(n_signal,
                         n_noise,
                         strong_signal_mu,
                         median_signal_mu,
                         noise_mu,
                         Sigma,
                         n_case,
                         ...){
    
    mu_strong_signal = rnorm(n_signal/2,strong_signal_mu,1)
    
    mu_median_signal = rnorm(n_signal/2,median_signal_mu,1)
    
    mu_noise = rnorm(n_noise,noise_mu,1)
    
    Mu = c(mu_strong_signal,mu_median_signal,mu_noise)
    
    case_data = MASS::mvrnorm(n_case,mu = Mu,Sigma = Sigma)
    
    colnames(case_data) = stringr::str_c("V", 1001:(1000 + n_signal + n_noise))
    
    return(case_data)
}
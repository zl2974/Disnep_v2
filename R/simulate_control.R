#### simulate control data
simulate_control = function(n_signal,
                            n_noise,
                            noise_mu,
                            noise_corr,
                            sigma,
                            n_control,
                            ...){
    
    mu_noise = rnorm(n_signal+n_noise,noise_mu,1)
    
    Mu = c(mu_noise)
    
    sigma_noise = matrix(runif((n_signal+n_noise)^2,min = 0,max = noise_corr),n_signal+n_noise)
    
    Sigma = as.matrix(Matrix::bdiag(sigma_noise))
    
    diag(Sigma) = 1
    
    control_data = MASS::mvrnorm(n_control,mu = Mu,Sigma = sigma^2*Sigma)
    
    colnames(control_data) = stringr::str_c("V", 1001:(1000 + n_signal +
                                                           n_noise))
    
    return(control_data)
}
#### simulate control data
simulate_control = function(n_signal,
                            n_noise,
                            noise_mu,
                            Sigma,
                            n_control,
                            ...){
    
    mu_noise = rnorm(n_signal+n_noise,noise_mu,1)
    
    Mu = c(mu_noise)
    
    control_data = MASS::mvrnorm(n_control,mu = Mu,Sigma = Sigma)
    
    colnames(control_data) = stringr::str_c("V", 1001:(1000 + n_signal +
                                                           n_noise))
    
    return(control_data)
}
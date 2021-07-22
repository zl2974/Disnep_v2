### simulate case and control gene expression data

simulate_gene_expr_data = function(n_signal,
                         n_noise,
                         strong_signal_mu,
                         median_signal_mu,
                         noise_mu,
                         strong_corr,
                         median_corr,
                         noise_corr,
                         sigma) {
    
    source("R/simulate_case.R")
    source("R/simulate_control.R")
    
    n_case = 300
    
    n_control = 50
    
    ### Generate Sigma
    
    sigma_strong = runif((n_signal/2),min = 1/strong_corr,max = 2/strong_corr)
    
    sigma_median = runif((n_signal/2),min = 1/median_corr,max = 2/median_corr)
    
    sigma_noise = runif(n_noise,min = 1/noise_corr,max = 2/noise_corr)
    
    Sigma = as.matrix(Matrix::bdiag(diag(sigma_strong), diag(sigma_median), diag(sigma_noise)))
    
    Sigma = Sigma + matrix(1, (n_signal + n_noise), (n_signal + n_noise)) 
    
    Sigma = cov2cor(Sigma)
    
    Sigma = sigma^2*Sigma
    
    case_data = simulate_case(
        n_signal,
        n_noise,
        strong_signal_mu,
        median_signal_mu,
        noise_mu,
        Sigma,
        n_case
    )
    
    control_data = simulate_control(
        n_signal,
        n_noise,
        noise_mu,
        Sigma,
        n_control
    )
    
    
    return(list(case = case_data,
                control = control_data))
    
    
}

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
    
    case_data = simulate_case(
        n_signal,
        n_noise,
        strong_signal_mu,
        median_signal_mu,
        noise_mu,
        strong_corr,
        median_corr,
        noise_corr,
        sigma,
        n_case
    )
    
    control_data = simulate_control(n_signal,
                                    n_noise,
                                    noise_mu,
                                    noise_corr,
                                    sigma,
                                    n_control)
    
    
    return(list(case = case_data,
                control = control_data))
    
    
}

##### whole simulation process####
simulation = function() {
    library(MASS)
    library(tidyverse)
    
    
    id = 1
    n_signal = 50
    n_noise = 950
    strong_signal_mu = 10
    median_signal_mu = 9
    noise_mu         = 8
    strong_corr      = 0.1
    median_corr      = 0.02
    noise_corr       = 0.01
    sigma            = 1
    
    settings = tibble(
        n_signal = n_signal,
        n_noise = n_noise,
        strong_signal_mu = strong_signal_mu,
        median_signal_mu = median_signal_mu,
        noise_mu = noise_mu,
        strong_corr = strong_corr,
        median_corr = median_corr,
        noise_corr = noise_corr,
        sigma = sigma
    )
    
    #### generate a simulation data####
    
    source("R/simulate_gene_expr_data.R")
    
    gene_expr_data = pmap(settings, simulate_gene_expr_data)[[1]]
    
    
    #### perform t-test           #####
    
    source("R/t_test.R")
    
    signals = t_test(gene_expr_data$case, gene_expr_data$control)
    
    #### integrate case and control data
    
    gene_expr_data = do.call(rbind, gene_expr_data)
    
    
    #### simulate gene interaction data
    
    source("R/simulate_gene_int_data.R")
    
    gene_int_data = simulate_gene_int_data(n_signal,n_noise)
    
    rownames(gene_int_data) = rownames(gene_expr_data)
    colnames(gene_int_data) = colnames(gene_expr_data)
    
    
    ####  perform Score Prioritized ###
    
    ## GeneWander
    
    source("R/genewanderer.R")
    score_genewanderer = genewanderer(signals,gene_int_data)
    saveRDS(score_genewanderer,paste0("cache/",id,"genewanderer.Rdata"))
    
    
}
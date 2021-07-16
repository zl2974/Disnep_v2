##### whole simulation process####
simulation = function(id,
                      n_signal,
                      n_noise,
                      strong_signal_mu,
                      median_signal_mu,
                      noise_mu,
                      strong_corr,
                      median_corr,
                      noise_corr,
                      sigma,
                      ...) {
    
    if (F) {
        n_signal = 20
        n_noise = 980
        strong_signal_mu = 10
        median_signal_mu = 9
        noise_mu = 8
        strong_corr = 0.1
        median_corr = 0.02
        noise_corr = 0.01
        sigma = 1
    }
    
    
    
    library(MASS)
    library(tidyverse)
    
    system("mkdir cache/result cache/simulation -p")
    
    setting = tibble(
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
    
    gene_expr_data = pmap(setting, simulate_gene_expr_data)[[1]]
    
    
    #### perform t-test           #####
    
    source("R/t_test.R")
    
    signals = t_test(gene_expr_data$case, gene_expr_data$control, .return = "statistics")
    
    #### integrate case and control data
    
    gene_expr_data = do.call(rbind, gene_expr_data)
    
    
    #### simulate gene interaction data
    
    source("R/simulate_gene_int_data.R")
    
    gene_int_data = simulate_gene_int_data(n_signal, n_noise)
    
    if(ncol(gene_expr_data)!=ncol(gene_int_data)) stop("networks dimension not match")
    
    rownames(gene_int_data) = colnames(gene_expr_data)
    colnames(gene_int_data) = colnames(gene_expr_data)
    
    
    ####  perform Score Prioritized ###
    
    ## T-test
    score_ttest = signals %>%
        arrange(gene) %>% 
        nest(score = everything()) %>%
        bind_cols(setting) %>%
        mutate(id = id,
               method = "t-test")
    
    saveRDS(score_ttest,
            paste0("cache/result/", Sys.Date(), "_", id,  "_t_test.Rdata"))
    
    
    ## GeneWander
    
    source("R/genewanderer.R")
    score_genewanderer = genewanderer(signals, gene_int_data, type = "statistics") %>%
        arrange(gene) %>% 
        nest(score = everything()) %>%
        bind_cols(setting) %>%
        mutate(id = id,
               method = "Genewanderer")
    
    saveRDS(
        score_genewanderer,
        paste0(
            "cache/result/",
            Sys.Date(),
            "_",
            id,
            "_genewanderer.Rdata"
        )
    )
    
    
    ## DISNEP
    
    source("R/disnep.R")
    score_disnep = disnep(signals, gene_int_data, gene_expr_data, type = "statistics") %>%
        arrange(gene) %>% 
        nest(score = everything()) %>%
        bind_cols(setting) %>%
        mutate(id = id,
               method = "Disnep")
    
    saveRDS(score_disnep,
            paste0("cache/result/", Sys.Date(), "_", id, "_disnep.Rdata"))
    
    
    
    #### save simulated data with setting
    
    saveRDS(setting %>% 
                mutate(gene_int = map(NA,~gene_int_data),
                       gene_expr = map(NA,~gene_expr_data)),
            paste0("cache/simulation/", Sys.Date(), "_", id, "_simulate_data.Rdata")
            )

}
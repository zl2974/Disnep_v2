### t-test between data

t_test = function(data_1, data_2, .return = "pval") {
    
    if (ncol(data_1) != ncol(data_2))
        stop("dimension not match")
    
    score = c()
    
    for (i in 1:ncol(data_1)) {
        score[[i]] = t.test(data_1[, i], data_2[, i])$p.value
    }
    
    score = as.numeric(score)
    
    if (.return != "pval")
        score = qnorm(1 - score)
    
    t_score = tibble(gene = colnames(data_1),
                     score = score)
    
    return(t_score)
}

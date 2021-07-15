#### simulate gene enteraction data

simulate_gene_int_data = function(n_signal,
                                  n_noise,
                                  ...){
    
    library(tidyverse)
    
    disease = 
        read_tsv("data/C1378703_disease_gda_summary.tsv")$Gene
    
    protein_gene_name = readr::read_csv("data/String_Network_default_node_04.csv") %>% 
        select(gene = `display name` , protein = `stringdb::database identifier`) %>% 
        filter(!str_detect(gene,"ENSP"))
    
    ppi = read.csv("data/9606.protein.links.v11.0.txt",sep = "") %>% 
        left_join(protein_gene_name,by =c('protein1' = 'protein')) %>% 
        left_join(protein_gene_name,by =c('protein2' = 'protein')) %>% 
        select(gene_1 = gene.x,gene_2=gene.y,score = combined_score)
    
    gene_pool = unique(ppi$gene_1)
    
    disease_gene = sort(sample(intersect(gene_pool,disease),n_signal))
    
    noise_gene = sort(sample(gene_pool[!gene_pool %in% disease], n_noise))
    
    gene_name = c(disease_gene,noise_gene)
    
    gene_int = expand.grid(gene_1 = gene_name,
                           gene_2 = gene_name) %>%
        left_join(ppi,by = c("gene_1", "gene_2")) %>% 
        pivot_wider(names_from = gene_2,
                    values_from = score,
                    names_sort = T
                    ) %>% 
        arrange(gene_1) %>% 
        select(-gene_1) %>% 
        as.matrix()
    
    rownames(gene_int) = sort(gene_name)
    
    gene_int[is.na(gene_int)] = 0
    
    gene_int = gene_int[gene_name,gene_name] # make sure first n_signal is signal gene
    
    if(norm(gene_int-t(gene_int))>0) stop("gene interaction not symmetric")
    
    if (any(c(
        rownames(gene_int)[1:n_signal] %in% noise_gene,
        colnames(gene_int)[1:n_signal] %in% noise_gene,
        nrow(gene_int) != (n_signal + n_noise),
        ncol(gene_int) != (n_signal + n_noise)
    )))
        stop("gene interaction format not correct")
    
    return(gene_int)

}

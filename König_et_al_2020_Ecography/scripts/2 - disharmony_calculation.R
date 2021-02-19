library(tidyverse)

setwd(".") # Set WD here

load("../Data/source_matrix.RData")
load("../Data/raw_data/checklists.RData")
load("../Data/raw_data/geoentities_meta.RData")

############################
### Calculate disharmony ###
############################
islands = geoentities_meta$entity_ID[which(geoentities_meta$entity_class == "Island")]

D_family = plyr::adply(islands, 1, function(island){
  # Summarize island flora
  tmp_checklist = checklists[which(checklists$entity_ID == island),]
  tmp_comm = table(tmp_checklist[,"family"])
  n_i = nrow(tmp_checklist)
  
  # Summarize source regions and calculate source region probabilities
  tmp_sources = source_matrix[paste(island),]
  
  s_ij = data.frame(entity_ID = as.numeric(colnames(source_matrix)), s_ij = tmp_sources / sum(tmp_sources)) %>% # Source region probabilities
    filter(s_ij > 0)
  
  n_j = checklists %>% # Species richness of source regions
    filter(entity_ID %in% s_ij$entity_ID) %>% 
    group_by(entity_ID) %>% 
    dplyr::summarise(n_j = n())
  
  p_tj =  checklists %>% # Proportion of taxa in source regions
    inner_join(n_j, by = "entity_ID") %>% 
    group_by(family, entity_ID) %>% 
    dplyr::summarise(p_tj = n() / n_j[1])
  
  p_ti = p_tj %>% inner_join(s_ij, by = "entity_ID") %>% inner_join(n_j, by = "entity_ID") %>% # Expected proportion of taxa on island i
    group_by(family) %>% 
    dplyr::summarise(p_ti = sum(s_ij * p_tj))
  
  # Calculate observed number of species per taxon
  n_ti = tmp_checklist %>% 
    group_by(family) %>% 
    dplyr::summarise(n_ti = n())
  
  # Calculate family-wise 
  result = full_join(p_ti, n_ti, by = "family") %>% 
    replace(is.na(.), 0) %>% 
    mutate(entity_ID = island, 
           n_obs = n_ti,
           n_exp = p_ti * n_i,
           n_total = n_i,
           D = ifelse(n_exp < 1, NA, pbinom(q = n_obs-1, size = n_i, prob = p_ti) + dbinom(n_obs, size = n_i, prob = p_ti) * 0.5)) %>%  # removing cases with n_exp < 1 is necessary to prevent weird behaviour of p-values in descrete distributions when the expected value is small
    dplyr::select(family, entity_ID, n_obs, n_exp, n_total, D)
}, .id = NULL)

save(list = c("D_family"), file = "../Data/D_family.RData")
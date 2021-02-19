library(raster)
library(tidyverse)
library(ade4)

rm(list = ls())
setwd("") # Set working directory here

# Run PCA at different resolutions
load("Data/env_extract.RData")

env_1000 = inner_join(climtopo_extr_1000, lidar_extr_1000) %>% drop_na() %>% dplyr::select(-samplearea_id)
pca_1000 = dudi.pca(env_1000, center = T, scale = T, scannf = F, nf = ncol(env_1000))
env_500 = inner_join(climtopo_extr_500, lidar_extr_500) %>% drop_na() %>% dplyr::select(-samplearea_id_500)
pca_500 = dudi.pca(env_500, center = T, scale = T, scannf = F, nf = ncol(env_500))
env_250 = inner_join(climtopo_extr_250, lidar_extr_250) %>% drop_na() %>% dplyr::select(-samplearea_id_250)
pca_250 = dudi.pca(env_250, center = T, scale = T, scannf = F, nf = ncol(env_250))
env_125 = inner_join(climtopo_extr_125, lidar_extr_125) %>% drop_na() %>% dplyr::select(-samplearea_id_125)
pca_125 = dudi.pca(env_125, center = T, scale = T, scannf = F, nf = ncol(env_125))

# Quantify loadings and get most important vars per PC
get_var_ranking = function(pca){
  x = pca$c1
  variables = vector(mode = "character", length = ncol(x))
  for(i in 1:ncol(x)){
    variables[i] = rownames(x)[which.max(abs(x[,i]))]
    x = subset(x, rownames(x) != variables[i])
  }
  return(data.frame(var_name = variables, rank = 1:length(variables), stringsAsFactors = F))
}

rankings = list("1000" = get_var_ranking(pca_1000), "500" = get_var_ranking(pca_500),
                "250" = get_var_ranking(pca_250), "125" = get_var_ranking(pca_125))

# Reduce to uncorrelated variables using select07 function from Dorman et al 2013
select07 = function(ranking, X, threshold = 0.7, method = "spearman"){
  cm = cor(X, method = method)
  
  pairs = which(abs(cm) >= threshold, arr.ind = T) # identifies correlated variable pairs
  pairs = pairs[pairs[,1] != pairs[,2],] # removes diagonal entries
  
  exclude = NULL
  for(i in seq_len(length(ranking))){
    if(ranking[i] %in% row.names(pairs) & !(ranking[i] %in% exclude)){
      cv = cm[setdiff(row.names(cm), exclude), ranking[i]]
      cv = cv[setdiff(names(cv), ranking[1:i])]
      exclude = c(exclude, names(which((abs(cv) >= threshold)))) 
    }
  }
  return(ranking[!(ranking %in% exclude), drop=F])
}

var_select = list("1000" = select07(rankings$`125`$var_name, env_125), "500" = select07(rankings$`250`$var_name, env_250),
                  "250" = select07(rankings$`500`$var_name, env_500), "125" = select07(rankings$`1000`$var_name, env_1000))

var_select_final = reduce(var_select, intersect) # find shared vars
save(var_select_final, file = "Data/var_select_final.RData")

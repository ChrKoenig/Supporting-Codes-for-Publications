library(tidyverse)
library(Hmsc)
library(coda)
library(corrplot)
library(pROC)

rm(list = ls())
setwd("") # Set working directory here

load("Data/tm_points.RData")
load("Data/sample_areas.RData")
load("Data/var_select_final.RData")
load("Data/env_extract.RData")
load("Data/final_species.RData")
load("Models/jsdm_fit_1.RData") # list of fitted models
load("Models/jsdm_fit_5.RData") # list of fitted models

##################
# Create coda object
model_results_coda_1 = lapply(jsdm_model_results_1, convertToCodaObject)
model_results_coda_5 = lapply(jsdm_model_results_5, convertToCodaObject)

###################################
# Trace plots
traceplot_random = function(object_name, model_name, par_name, n_samples){
  tmp_object = get(object_name)
  tmp_coda = tmp_object[[model_name]][[par_name]]
  if(class(tmp_coda) != "mcmc.list"){
    for(i in 1:length(tmp_coda)){
      tmp_coda_i = tmp_coda[[i]]
      samples = sample(varnames(tmp_coda_i), n_samples, replace = F)
      tmp_coda_i = tmp_coda_i[,samples]
      plot(tmp_coda_i, auto.layout = T)
    }
  } else {
    samples = sample(varnames(tmp_coda), n_samples, replace = F)
    tmp_coda = tmp_coda[,samples]
    plot(tmp_coda, auto.layout = T)
  }
}

par(mar = c(1,2,1,1))
# Beta @ 1000m
traceplot_random("model_results_coda_1", "model_occ_1000", "Beta", 4)
traceplot_random("model_results_coda_1", "model_abund_1000", "Beta", 4)
traceplot_random("model_results_coda_5", "model_occ_1000", "Beta", 4)
traceplot_random("model_results_coda_5", "model_abund_1000", "Beta", 4)

# Beta @ 500m
traceplot_random("model_results_coda_1", "model_occ_500", "Beta", 4)
traceplot_random("model_results_coda_1", "model_abund_500", "Beta", 4)
traceplot_random("model_results_coda_5", "model_occ_500", "Beta", 4)
traceplot_random("model_results_coda_5", "model_abund_500", "Beta", 4)

# Beta @ 250m
traceplot_random("model_results_coda_1", "model_occ_250", "Beta", 4)
traceplot_random("model_results_coda_1", "model_abund_250", "Beta", 4)
traceplot_random("model_results_coda_5", "model_occ_250", "Beta", 4)
traceplot_random("model_results_coda_5", "model_abund_250", "Beta", 4)

# Beta @ 125m
traceplot_random("model_results_coda_1", "model_occ_125", "Beta", 4)
traceplot_random("model_results_coda_1", "model_abund_125", "Beta", 4)
traceplot_random("model_results_coda_5", "model_occ_125", "Beta", 4)
traceplot_random("model_results_coda_5", "model_abund_125", "Beta", 4)

# Omega @ 1000m
traceplot_random("model_results_coda_1", "model_occ_1000", "Omega", 4)
traceplot_random("model_results_coda_1", "model_abund_1000", "Omega", 4)
traceplot_random("model_results_coda_5", "model_occ_1000", "Omega", 4)
traceplot_random("model_results_coda_5", "model_abund_1000", "Omega", 4)

# Omega @ 500m
traceplot_random("model_results_coda_1", "model_occ_500", "Omega", 4)
traceplot_random("model_results_coda_1", "model_abund_500", "Omega", 4)
traceplot_random("model_results_coda_5", "model_occ_500", "Omega", 4)
traceplot_random("model_results_coda_5", "model_abund_500", "Omega", 4)

# Omega @ 250m
traceplot_random("model_results_coda_1", "model_occ_250", "Omega", 4)
traceplot_random("model_results_coda_1", "model_abund_250", "Omega", 4)
traceplot_random("model_results_coda_5", "model_occ_250", "Omega", 4)
traceplot_random("model_results_coda_5", "model_abund_250", "Omega", 4)

# Omega @ 125m
traceplot_random("model_results_coda_1", "model_occ_125", "Omega", 4)
traceplot_random("model_results_coda_1", "model_abund_125", "Omega", 4)
traceplot_random("model_results_coda_5", "model_occ_125", "Omega", 4)
traceplot_random("model_results_coda_5", "model_abund_125", "Omega", 4)

rm(model_results_coda_1, model_results_coda_5)

#######################
# Effective sample size
# 1 year
model_conv_ESS_1year = lapply(model_results_coda_1, function(x){  # Ideal: ESS = SS
  beta = effectiveSize(x$Beta)
  omega = lapply(x$Omega, effectiveSize)
  names(omega) = paste("omega", 1:length(omega), sep = "_")
  ESS = c(list("beta" = beta), omega)
  return(ESS)
})

pdf("Figures/ESS_beta_1year.pdf", width = 18, height = 10)
layout(matrix(1:length(model_conv_ESS_1year), nrow = 2))
for(i in 1:length(model_conv_ESS_1year)){
  hist(model_conv_ESS_1year[[i]]$beta, main = paste("Beta", names(model_conv_ESS_1year)[i]))
  abline(v = jsdm_model_results_1[[i]]$samples * length(jsdm_model_results_1[[i]]$repList), col = "red")
}
dev.off()

# 5 years
model_conv_ESS_5year = lapply(model_results_coda_5, function(x){  # Ideal: ESS = SS
  beta = effectiveSize(x$Beta)
  omega = lapply(x$Omega, effectiveSize)
  names(omega) = paste("omega", 1:length(omega), sep = "_")
  ESS = c(list("beta" = beta), omega)
  return(ESS)
})
save(model_conv_ESS_5year, file = "Data/model_conv_ESS_5year.RData")

pdf("Figures/ESS_beta_5years.pdf", width = 13.5, height = 10)
layout(matrix(1:length(model_conv_ESS_5year), nrow = 2))
for(i in 1:length(model_conv_ESS_5year)){
  hist(model_conv_ESS_5year[[i]]$beta, main = paste("Beta", names(model_conv_ESS_5year)[i]))
  abline(v = jsdm_model_results_5[[i]]$samples * length(jsdm_model_results_5[[i]]$repList), col = "red")
}
dev.off()

#####################
# Parameter estimates
# Betas (regression coefficients)
beta_estimates_1year = lapply(jsdm_model_results_1, getPostEstimate, parName = "Beta") 
beta_estimates_5year = lapply(jsdm_model_results_5, getPostEstimate, parName = "Beta") 

par(mar = c(4,4,1,1))
for(i in 1:length(jsdm_model_results_1)){
  tmp_name = names(jsdm_model_results_1)[i]
  png(filename = paste0("Figures/beta_", tmp_name, "_1year.png"), width = 1000, height = 700)
  plotBeta(jsdm_model_results_1[[i]], beta_estimates_1year[[i]], "Mean", support = 0.95)
  dev.off()
}

for(i in 1:length(jsdm_model_results_5)){
  tmp_name = names(jsdm_model_results_5)[i]
  png(filename = paste0("Figures/beta_", tmp_name, "_5years.png"), width = 1000, height = 700)
  plotBeta(jsdm_model_results_5[[i]], beta_estimates_5year[[i]], "Mean", support = 0.95)
  dev.off()
}

# Omegas (Residual correlations)
for(i in 1:length(jsdm_model_results_1)){
  tmp_name = names(jsdm_model_results_1)[i] # model name
  tmp_omega = computeAssociations(jsdm_model_results_1[[i]]) # Resiudals at each random effects
  tmp_omega = tmp_omega[[length(tmp_omega)]] # look only at sample-level random effect
  tmp_omega_plot = ((tmp_omega$support > 0.95) + (tmp_omega$support < (1 - 0.95)) > 0) * tmp_omega$mean
  png(filename = paste0("Figures/omega_", tmp_name, "_1year.png"), width = 1000, height = 700)
  corrplot(tmp_omega_plot, method = "color", diag = F, type = "lower", tl.col = "black",
           col = colorRampPalette(c("blue","white","red"))(200), mar=c(0,0,1,0))
  dev.off()
}

for(i in 1:length(jsdm_model_results_5)){
  tmp_name = names(jsdm_model_results_5)[i] # model name
  tmp_omega = computeAssociations(jsdm_model_results_5[[i]]) # Resiudals at each random effects
  tmp_omega = tmp_omega[[length(tmp_omega)-1]] # look only at sample-level random effect
  tmp_omega_plot = ((tmp_omega$support > 0.95) + (tmp_omega$support < (1 - 0.95)) > 0) * tmp_omega$mean
  png(filename = paste0("Figures/omega_", tmp_name, "_5year.png"), width = 1000, height = 700)
  corrplot(tmp_omega_plot, method = "color", diag = F, type = "lower", tl.col = "black",
           col = colorRampPalette(c("blue","white","red"))(200), mar=c(0,0,1,0))
  dev.off()
}

######################################
# Variance partitioning and model fits
# CAUTION: This requires a machine with a lot of RAM (> 100GB)
# 1-year models
jsdm_model_results_1 = subset(jsdm_model_results_1, names(jsdm_model_results_1) %in% c("model_occ_1000", "model_occ_500", "model_occ_250", "model_occ_125"))

# model fits
for(i in 1:length(jsdm_model_results_1)){
  print(i)
  tmp_name = names(jsdm_model_results_1)[i]
  tmp_model = jsdm_model_results_1[[i]]
  pred = computePredictedValues(tmp_model, thin = 100)
  tmp_fit = evaluateModelFit(tmp_model, pred)
  saveRDS(tmp_fit, file = paste0(tmp_name,"_fit_1.RDS"))
  gc()
} 

# variance partitioning
for(i in 1:length(jsdm_model_results_1)){
  print(i)
  tmp_name = names(jsdm_model_results_1)[i]
  tmp_model = jsdm_model_results_1[[i]]
  tmp_varPart = computeVariancePartitioning(tmp_model, group = rep(1, ncol(tmp_model$X)), groupnames = c("fixed"))
  saveRDS(tmp_varPart, file = paste0(tmp_name,"_varPart_1.RDS"))
  gc()
}

# 5-year models
jsdm_model_results_5 = subset(jsdm_model_results_5, names(jsdm_model_results_5) %in% c("model_occ_1000", "model_occ_500", "model_occ_250", "model_occ_125"))

# model fits
for(i in 1:length(jsdm_model_results_5)){
  print(i)
  tmp_name = names(jsdm_model_results_5)[i]
  tmp_model = jsdm_model_results_5[[i]]
  pred = computePredictedValues(tmp_model, thin = 100)
  tmp_fit = evaluateModelFit(tmp_model, pred)
  saveRDS(tmp_fit, file = paste0(tmp_name,"_fit_5.RDS"))
  gc()
} 

# variance partitioning
for(i in 1:length(jsdm_model_results_5)){
  print(i)
  tmp_name = names(jsdm_model_results_5)[i]
  tmp_model = jsdm_model_results_5[[i]]
  tmp_varPart = computeVariancePartitioning(tmp_model, group = rep(1, ncol(tmp_model$X)), groupnames = c("fixed"), )
  saveRDS(tmp_varPart, file = paste0(tmp_name, "_varPart_5.RDS"))
  gc()
}

######################################
#     Predictive Performance         #
######################################
# 1-year models
# Get unique site ids
sample_areas_1year = sample_areas %>% 
  rename(samplearea_id_1000 = samplearea_id) %>% 
  filter(samplearea_id_1000 %in% jsdm_model_results_1$model_occ_1000$studyDesign$samplearea_id) %>% 
  mutate_at(vars(contains("samplearea_id")), ~as.character(.))

# Get unique species ids
species_ids = jsdm_model_results_1$model_occ_1000$spNames

# Subset observations
holdout_obs = tm_points %>%
  filter(year %in% 2010:2012 & samplearea_id %in% sample_areas_1year$samplearea_id_1000 & species_id %in% final_species) %>%
  mutate(samplearea_id_1000 = factor(samplearea_id, levels = unique(sample_areas_1year$samplearea_id_1000)),
         samplearea_id_500 = factor(samplearea_id_500, levels = unique(sample_areas_1year$samplearea_id_500)),
         samplearea_id_250 = factor(samplearea_id_250, levels = unique(sample_areas_1year$samplearea_id_250)),
         samplearea_id_125 = factor(samplearea_id_125, levels = unique(sample_areas_1year$samplearea_id_125)))

pred_performance_1year = lapply(c(1000, 500, 250, 125), function(r){
  cat("\n", r)
  model_tmp = jsdm_model_results_1[[paste0("model_occ_", r)]]
  predictions_tmp = computePredictedValues(model_tmp, thin = 50)
  predictions_tmp_mean = apply(predictions_tmp, c(1,2), mean)
  rownames(predictions_tmp_mean) = gsub("_2009", "", rownames(predictions_tmp_mean))
  performance_r = lapply(2010:2012, function(y){
    cat("   ", y)
    observations_tmp = holdout_obs %>% 
      filter(year == y) %>% 
      group_by(species_id, samplearea_id = get(paste0("samplearea_id_", r))) %>% 
      tally() %>% 
      mutate(n = ifelse(n > 0, 1, 0)) %>% 
      complete(species_id, samplearea_id, fill = list(n = 0)) %>% 
      mutate(samplearea_id = as.character(samplearea_id)) %>% 
      pivot_wider(id_cols = samplearea_id, names_from = species_id, values_from = n) %>% 
      column_to_rownames("samplearea_id")
    observations_tmp = observations_tmp[rownames(predictions_tmp_mean),]
    
    AUC_y = suppressMessages(sapply(species_ids, function(x){
      pROC::auc(observations_tmp[,x] ~ predictions_tmp_mean[,x])
    }))
    
    R2_y = suppressMessages(sapply(species_ids, function(x){
      pred_pos = mean(predictions_tmp_mean[which(observations_tmp[,x] == 1), x])
      pred_neg = mean(predictions_tmp_mean[which(observations_tmp[,x] == 0), x])
      return(pred_pos - pred_neg)
    }))
    
    data.frame(species_id = c(names(AUC_y), names(R2_y)),
               year = y,
               metric = c(rep("AUC", length(AUC_y)), rep("TjurR2", length(R2_y))),
               value = c(AUC_y, R2_y))
  })
  performance_r = bind_rows(performance_r) %>% 
    mutate(grain_size = r) %>% 
    select(species_id, grain_size, year, metric, value)
})

pred_performance_1year = bind_rows(pred_performance_1year) %>% 
  mutate(temp_extent = "1", year = year-2009) %>% 
  pivot_wider(values_from = value, names_from = year, names_prefix = "pred_")

# 5-year models
# Get unique site ids
sample_areas_5years = sample_areas %>% st_drop_geometry() %>% 
  rename(samplearea_id_1000 = samplearea_id) %>% 
  filter(samplearea_id_1000 %in% jsdm_model_results_5$model_occ_1000$studyDesign$samplearea_id) %>% 
  mutate_at(vars(contains("samplearea_id")), ~as.character(.))

# Get unique species ids
species_ids = jsdm_model_results_5$model_occ_1000$spNames

# Subset observations
holdout_obs = tm_points %>%
  st_drop_geometry() %>% 
  filter(year %in% 2012:2014 & samplearea_id %in% sample_areas_5years$samplearea_id_1000 & species_id %in% species_ids) %>%
  mutate(samplearea_id_1000 = factor(samplearea_id, levels = unique(sample_areas_5years$samplearea_id_1000)),
         samplearea_id_500 = factor(samplearea_id_500, levels = unique(sample_areas_5years$samplearea_id_500)),
         samplearea_id_250 = factor(samplearea_id_250, levels = unique(sample_areas_5years$samplearea_id_250)),
         samplearea_id_125 = factor(samplearea_id_125, levels = unique(sample_areas_5years$samplearea_id_125)))

pred_performance_5years = lapply(c(1000, 500, 250, 125), function(r){
  cat("\n", r)
  model_tmp = jsdm_model_results_5[[paste0("model_occ_", r)]]
  predictions_tmp = computePredictedValues(model_tmp, thin = 50)
  predictions_tmp_mean = apply(predictions_tmp, c(1,2), mean)
  rownames(predictions_tmp_mean) = gsub("_.*$", "", rownames(predictions_tmp_mean))
  performance_r = lapply(2012:2014, function(y){
    cat("   ", y)
    observations_tmp = holdout_obs %>% 
      filter(year == y) %>% 
      group_by(species_id, samplearea_id = get(paste0("samplearea_id_", r))) %>% 
      tally() %>% 
      mutate(n = ifelse(n > 0, 1, 0)) %>% 
      complete(species_id, samplearea_id, fill = list(n = 0)) %>% 
      mutate(samplearea_id = as.character(samplearea_id)) %>% 
      pivot_wider(id_cols = samplearea_id, names_from = species_id, values_from = n) %>% 
      column_to_rownames("samplearea_id")
    observations_tmp = observations_tmp[rownames(predictions_tmp_mean),]
    
    AUC_y = suppressMessages(sapply(species_ids, function(x){
      tryCatch({
        pROC::auc(pull(observations_tmp,x) ~ predictions_tmp_mean[,x])
      }, error = function(e){
        return(NA)
      })
    }))
    
    R2_y = suppressMessages(sapply(species_ids, function(x){
      pred_pos = mean(predictions_tmp_mean[which(observations_tmp[,x] == 1), x])
      pred_neg = mean(predictions_tmp_mean[which(observations_tmp[,x] == 0), x])
      return(pred_pos - pred_neg)
    }))
    
    data.frame(species_id = c(names(AUC_y), names(R2_y)),
               year = y,
               metric = c(rep("AUC", length(AUC_y)), rep("TjurR2", length(R2_y))),
               value = c(AUC_y, R2_y))
  })
  performance_r = bind_rows(performance_r) %>% 
    mutate(grain_size = r) %>% 
    select(species_id, grain_size, year, metric, value)
})

pred_performance_5years = bind_rows(pred_performance_5years) %>% 
  mutate(temp_extent = "5", year = year-2011) %>% 
  pivot_wider(values_from = value, names_from = year, names_prefix = "pred_")

pred_performance = bind_rows(pred_performance_1year, pred_performance_5years)
save(pred_performance, file = "Data/pred_performance.RData")
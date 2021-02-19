library(tidyverse)
library(Hmsc)
library(raster)
library(rslurm)

rm(list = ls())
setwd("") # Set working directory here

load("Data/final_species.RData")
load("Data/tm_points.RData")
load("Data/sample_areas.RData")
load("Data/var_select_final.RData")
load("Data/env_extract.RData")
load("Data/forest_cells.RData")
load("Data/forest_species.RData")

#########################################
#        Prepare 1-year models         #
#########################################
sampling_period = 2009

# 1. 1000m
Y_1000_occ = tm_points %>%
  filter(year %in% sampling_period & samplearea_id %in% forest_cells & species_id %in% final_species) %>%
  group_by(species_id, samplearea_id, year) %>%
  tally() %>%
  mutate(present = 1, row_names = paste(samplearea_id, year, sep = "_")) %>%
  dplyr::select(-n) %>% 
  spread(key = species_id, value = present, fill = 0) %>%
  column_to_rownames(var = "row_names") 

Y_1000_abund = tm_points %>% 
  filter(year %in% sampling_period & samplearea_id %in% forest_cells & species_id %in% final_species) %>% 
  group_by(species_id, samplearea_id, year) %>% 
  tally() %>% 
  mutate(row_names = paste(samplearea_id, year, sep = "_")) %>% 
  spread(key = species_id, value = n, fill = 0) %>% 
  column_to_rownames(var = "row_names")

X_1000 = climtopo_extr_1000 %>%   
  left_join(lidar_extr_1000, by = "samplearea_id") %>% 
  filter(samplearea_id %in% Y_1000_abund$samplearea_id) %>% 
  dplyr::select(samplearea_id, var_select_final) %>% 
  mutate_at(vars(-samplearea_id), list(sq = ~.*.)) %>% # add quadratic terms
  mutate_at(vars(-samplearea_id), scale) %>% # center and recale variables
  mutate_at(vars(-samplearea_id), as.numeric) %>% #remove attributes from scaled vars
  right_join(Y_1000_abund[,c("samplearea_id", "year")], by = "samplearea_id") %>% 
  mutate(row_names = paste(samplearea_id, year, sep = "_")) %>% 
  dplyr::select(-samplearea_id, -year) %>% 
  column_to_rownames(var = "row_names")

all(rownames(X_1000) == rownames(Y_1000_occ))
all(rownames(X_1000) == rownames(Y_1000_abund))

Y_1000_occ = dplyr::select(Y_1000_occ, -samplearea_id, -year)
Y_1000_abund = dplyr::select(Y_1000_abund, -samplearea_id, -year)

study_design = str_split(rownames(X_1000), pattern = "_", simplify = T)[,1] %>% 
  as_tibble() %>% 
  mutate_each(funs = factor) %>% 
  set_names(c("samplearea_id"))
random_levels = list("samplearea_id" = HmscRandomLevel(units = unique(study_design$samplearea_id)))

model_occ_1000 = Hmsc(Y = as.matrix(Y_1000_occ), XData = X_1000, XFormula = ~., distr = "probit",
                      studyDesign = as.data.frame(study_design), ranLevels = random_levels)
model_abund_1000 = Hmsc(Y = as.matrix(Y_1000_abund), XData = X_1000, XFormula = ~., distr = "lognormal poisson",
                        studyDesign = as.data.frame(study_design), ranLevels = random_levels)

##########
# 2. 500m
Y_500_occ = tm_points %>%
  filter(year %in% sampling_period & samplearea_id %in% forest_cells & species_id %in% final_species) %>%
  group_by(species_id, samplearea_id_500, year) %>%
  tally() %>%
  mutate(present = 1, row_names = paste(samplearea_id_500, year, sep = "_")) %>%
  dplyr::select(-n) %>%
  spread(key = species_id, value = present, fill = 0) %>%
  column_to_rownames(var = "row_names")

Y_500_abund = tm_points %>% 
  filter(year %in% sampling_period & samplearea_id %in% forest_cells & species_id %in% final_species) %>% 
  group_by(species_id, samplearea_id_500, year) %>% 
  tally() %>% 
  mutate(row_names = paste(samplearea_id_500, year, sep = "_")) %>% 
  spread(key = species_id, value = n, fill = 0) %>% 
  column_to_rownames(var = "row_names")

X_500 = climtopo_extr_500 %>%   
  left_join(lidar_extr_500, by = "samplearea_id_500") %>% 
  filter(samplearea_id_500 %in% Y_500_occ$samplearea_id_500) %>% 
  dplyr::select(samplearea_id_500, var_select_final) %>% 
  mutate_at(vars(-samplearea_id_500), list(sq = ~.*.)) %>% # add quadratic terms
  mutate_at(vars(-samplearea_id_500), scale) %>% # center and recale variables
  mutate_at(vars(-samplearea_id_500), as.numeric) %>% # remove attributes from scaled vars
  right_join(Y_500_abund[,c("samplearea_id_500", "year")], by = "samplearea_id_500") %>% 
  mutate(row_names = paste(samplearea_id_500, year, sep = "_")) %>% 
  dplyr::select(-samplearea_id_500, -year) %>% 
  column_to_rownames(var = "row_names")

all(rownames(X_500) == rownames(Y_500_occ))
all(rownames(X_500) == rownames(Y_500_abund))

X_500 = X_500[complete.cases(X_500),]
Y_500_occ = dplyr::select(Y_500_occ, -samplearea_id_500, -year) %>% 
  filter(rownames(.) %in% rownames(X_500))
Y_500_abund = dplyr::select(Y_500_abund, -samplearea_id_500, -year) %>% 
  filter(rownames(.) %in% rownames(X_500))

study_design = str_split(rownames(X_500), pattern = "_", simplify = T)[,1] %>% 
  as_tibble() %>% 
  mutate_each(funs = as.integer) %>% 
  set_names(c("samplearea_id_500")) %>% 
  left_join(unique(sample_areas[,c("samplearea_id", "samplearea_id_500")]), by = "samplearea_id_500") %>% 
  dplyr::select(samplearea_id, samplearea_id_500) %>% 
  mutate_each(funs = as_factor)
random_levels = list("samplearea_id" = HmscRandomLevel(units = unique(study_design$samplearea_id)),
                     "samplearea_id_500" = HmscRandomLevel(units = unique(study_design$samplearea_id_500)))

model_occ_500 = Hmsc(Y = as.matrix(Y_500_occ), XData = X_500, XFormula = ~., distr = "probit",
                     studyDesign = as.data.frame(study_design), ranLevels = random_levels) 
model_abund_500 = Hmsc(Y = as.matrix(Y_500_abund), XData = X_500, XFormula = ~., distr = "lognormal poisson",
                       studyDesign = as.data.frame(study_design), ranLevels = random_levels) 

##########
# 3. 250m 
Y_250_occ = tm_points %>%
  filter(year %in% sampling_period & samplearea_id %in% forest_cells & species_id %in% final_species) %>%
  group_by(species_id, samplearea_id_250, year) %>%
  tally() %>%
  mutate(present = 1, row_names = paste(samplearea_id_250, year, sep = "_")) %>%
  dplyr::select(-n) %>%
  spread(key = species_id, value = present, fill = 0) %>%
  column_to_rownames(var = "row_names")

Y_250_abund = tm_points %>% 
  filter(year %in% sampling_period & samplearea_id %in% forest_cells & species_id %in% final_species) %>% 
  group_by(species_id, samplearea_id_250, year) %>% 
  tally() %>% 
  mutate(row_names = paste(samplearea_id_250, year, sep = "_")) %>% 
  spread(key = species_id, value = n, fill = 0) %>% 
  column_to_rownames(var = "row_names")

X_250 = climtopo_extr_250 %>%   
  left_join(lidar_extr_250, by = "samplearea_id_250") %>% 
  filter(samplearea_id_250 %in% Y_250_occ$samplearea_id_250) %>% 
  dplyr::select(samplearea_id_250, var_select_final) %>% 
  mutate_at(vars(-samplearea_id_250), list(sq = ~.*.)) %>% # add quadratic terms
  mutate_at(vars(-samplearea_id_250), scale) %>% # center and recale variables
  mutate_at(vars(-samplearea_id_250), as.numeric) %>% # remove attributes from scaled vars
  right_join(Y_250_abund[,c("samplearea_id_250", "year")], by = "samplearea_id_250") %>% 
  mutate(row_names = paste(samplearea_id_250, year, sep = "_")) %>% 
  dplyr::select(-samplearea_id_250, -year) %>% 
  column_to_rownames(var = "row_names")

all(rownames(X_250) == rownames(Y_250_occ))
all(rownames(X_250) == rownames(Y_250_abund))

X_250 = X_250[complete.cases(X_250),]
Y_250_occ = dplyr::select(Y_250_occ, -samplearea_id_250, -year) %>% 
  filter(rownames(.) %in% rownames(X_250))
Y_250_abund = dplyr::select(Y_250_abund, -samplearea_id_250, -year) %>% 
  filter(rownames(.) %in% rownames(X_250))

study_design = str_split(rownames(X_250), pattern = "_", simplify = T)[,1] %>% 
  as_tibble() %>% 
  mutate_each(funs = as.integer) %>% 
  set_names(c("samplearea_id_250")) %>% 
  left_join(unique(sample_areas[,c("samplearea_id", "samplearea_id_250")]), by = "samplearea_id_250") %>% 
  dplyr::select(samplearea_id, samplearea_id_250) %>% 
  mutate_each(funs = as_factor)
random_levels = list("samplearea_id" = HmscRandomLevel(units = unique(study_design$samplearea_id)),
                     "samplearea_id_250" = HmscRandomLevel(units = unique(study_design$samplearea_id_250)))

model_occ_250 = Hmsc(Y = as.matrix(Y_250_occ), XData = X_250, XFormula = ~., distr = "probit",
                     studyDesign = as.data.frame(study_design), ranLevels = random_levels) 
model_abund_250 = Hmsc(Y = as.matrix(Y_250_abund), XData = X_250, XFormula = ~., distr = "lognormal poisson",
                       studyDesign = as.data.frame(study_design), ranLevels = random_levels) 

#########
# 4. 125m
Y_125_occ = tm_points %>%
  filter(year %in% sampling_period & samplearea_id %in% forest_cells & species_id %in% final_species) %>%
  group_by(species_id, samplearea_id_125, year) %>%
  tally() %>%
  mutate(present = 1, row_names = paste(samplearea_id_125, year, sep = "_")) %>%
  dplyr::select(-n) %>%
  spread(key = species_id, value = present, fill = 0) %>%
  column_to_rownames(var = "row_names")

Y_125_abund = tm_points %>% 
  filter(year %in% sampling_period & samplearea_id %in% forest_cells & species_id %in% final_species) %>% 
  group_by(species_id, samplearea_id_125, year) %>% 
  tally() %>% 
  mutate(row_names = paste(samplearea_id_125, year, sep = "_")) %>% 
  spread(key = species_id, value = n, fill = 0) %>% 
  column_to_rownames(var = "row_names")

X_125 = climtopo_extr_125 %>%   
  left_join(lidar_extr_125, by = "samplearea_id_125") %>% 
  filter(samplearea_id_125 %in% Y_125_occ$samplearea_id_125) %>% 
  dplyr::select(samplearea_id_125, var_select_final) %>% 
  mutate_at(vars(-samplearea_id_125), list(sq = ~.*.)) %>% # add quadratic terms
  mutate_at(vars(-samplearea_id_125), scale) %>% # center and recale variables
  mutate_at(vars(-samplearea_id_125), as.numeric) %>% # remove attributes from scaled vars
  right_join(Y_125_abund[,c("samplearea_id_125", "year")], by = "samplearea_id_125") %>% 
  mutate(row_names = paste(samplearea_id_125, year, sep = "_")) %>% 
  dplyr::select(-samplearea_id_125, -year) %>% 
  column_to_rownames(var = "row_names")

all(rownames(X_125) == rownames(Y_125_occ))
all(rownames(X_125) == rownames(Y_125_abund))

X_125 = X_125[complete.cases(X_125),]
Y_125_occ = dplyr::select(Y_125_occ, -samplearea_id_125, -year) %>% 
  filter(rownames(.) %in% rownames(X_125))
Y_125_abund = dplyr::select(Y_125_abund, -samplearea_id_125, -year) %>% 
  filter(rownames(.) %in% rownames(X_125))

study_design = str_split(rownames(X_125), pattern = "_", simplify = T)[,1] %>% 
  as_tibble() %>% 
  mutate_each(funs = as.integer) %>% 
  set_names(c("samplearea_id_125")) %>% 
  left_join(unique(sample_areas[,c("samplearea_id", "samplearea_id_125")]), by = "samplearea_id_125") %>% 
  dplyr::select(samplearea_id, samplearea_id_125) %>% 
  mutate_each(funs = as_factor)
random_levels = list("samplearea_id" = HmscRandomLevel(units = unique(study_design$samplearea_id)),
                     "samplearea_id_125" = HmscRandomLevel(units = unique(study_design$samplearea_id_125)))

model_occ_125 = Hmsc(Y = as.matrix(Y_125_occ), XData = X_125, XFormula = ~., distr = "probit",
                     studyDesign = as.data.frame(study_design), ranLevels = random_levels) 
model_abund_125 = Hmsc(Y = as.matrix(Y_125_abund), XData = X_125, XFormula = ~., distr = "lognormal poisson",
                       studyDesign = as.data.frame(study_design), ranLevels = random_levels) 

####################################################
#           Fit 1-year models on cluster           #
####################################################
fit_model = function(model_name){
  tmp_model = model_list[[model_name]]
  tmp_model_fit = Hmsc::sampleMcmc(tmp_model, nChains = 8, thin = 50, samples = 500, transient = 5000, nParallel = 8)
  return(tmp_model_fit)
}
model_list = list("model_occ_1000" = model_occ_1000, "model_abund_1000" = model_abund_1000, 
                  "model_occ_500" = model_occ_500, "model_abund_500" = model_abund_500, 
                  "model_occ_250" = model_occ_250, "model_abund_250" = model_abund_250, 
                  "model_occ_125" = model_occ_125, "model_abund_125" = model_abund_125)
pars = data.frame(model_name = names(model_list), stringsAsFactors = F)

# CAUTION: This runs up to 20 days!
models_fit_1 <- slurm_apply(fit_model, # function to be parallelized
                            params = pars, # function parameters (one column = one parameter)
                            jobname = "fit_models_1year",
                            add_objects = c("model_list"),
                            slurm_options = list(qos = "long"),
                            nodes = 8, cpus_per_node = 1, # Requested resources 
                            submit = F) # Test rdc _un or submit?
# Submit = F because rslurm requests one core per row in 'pars', but we need 8 cores per row
# --> Change line "#SBATCH --cpus-per-task=1" to "#SBATCH --cpus-per-task=8" in submit.sh
# --> use terminal and navigate to newly created rslurm folder
# --> run command "sbatch submit.sh"

print_job_status(models_fit_1)

##############################################
#           Prepare 5-year models            #
##############################################
sampling_period = 2007:2011

# 1. 1000m
Y_1000_occ = tm_points %>%
  filter(year %in% sampling_period & samplearea_id %in% forest_cells & species_id %in% final_species) %>%
  group_by(species_id, samplearea_id, year) %>%
  tally() %>%
  mutate(present = 1, row_names = paste(samplearea_id, year, sep = "_")) %>%
  dplyr::select(-n) %>% 
  spread(key = species_id, value = present, fill = 0) %>%
  column_to_rownames(var = "row_names") 

Y_1000_abund = tm_points %>% 
  filter(year %in% sampling_period & samplearea_id %in% forest_cells & species_id %in% final_species) %>% 
  group_by(species_id, samplearea_id, year) %>% 
  tally() %>% 
  mutate(row_names = paste(samplearea_id, year, sep = "_")) %>% 
  spread(key = species_id, value = n, fill = 0) %>% 
  column_to_rownames(var = "row_names")

X_1000 = climtopo_extr_1000 %>%   
  left_join(lidar_extr_1000, by = "samplearea_id") %>% 
  filter(samplearea_id %in% Y_1000_abund$samplearea_id) %>% 
  dplyr::select(samplearea_id, var_select_final) %>% 
  mutate_at(vars(-samplearea_id), list(sq = ~.*.)) %>% # add quadratic terms
  mutate_at(vars(-samplearea_id), scale) %>% # center and recale variables
  mutate_at(vars(-samplearea_id), as.numeric) %>% #remove attributes from scaled vars
  right_join(Y_1000_abund[,c("samplearea_id", "year")], by = "samplearea_id") %>% 
  mutate(row_names = paste(samplearea_id, year, sep = "_")) %>% 
  dplyr::select(-samplearea_id, -year) %>% 
  column_to_rownames(var = "row_names")

all(rownames(X_1000) == rownames(Y_1000_occ))
all(rownames(X_1000) == rownames(Y_1000_abund))

Y_1000_occ = dplyr::select(Y_1000_occ, -samplearea_id, -year)
Y_1000_abund = dplyr::select(Y_1000_abund, -samplearea_id, -year)

study_design = str_split(rownames(X_1000), pattern = "_", simplify = T) %>% 
  as_data_frame() %>% 
  mutate_each(funs = factor) %>% 
  set_names(c("samplearea_id", "year"))
random_levels = list("samplearea_id" = HmscRandomLevel(units = unique(study_design$samplearea_id)),
                     "year" = HmscRandomLevel(units = unique(study_design$year)))

model_occ_1000 = Hmsc(Y = as.matrix(Y_1000_occ), XData = X_1000, XFormula = ~., distr = "probit",
                      studyDesign = as.data.frame(study_design), ranLevels = random_levels)
model_abund_1000 = Hmsc(Y = as.matrix(Y_1000_abund), XData = X_1000, XFormula = ~., distr = "lognormal poisson",
                        studyDesign = as.data.frame(study_design), ranLevels = random_levels)

##########
# 2. 500m
Y_500_occ = tm_points %>%
  filter(year %in% sampling_period & samplearea_id %in% forest_cells & species_id %in% final_species) %>%
  group_by(species_id, samplearea_id_500, year) %>%
  tally() %>%
  mutate(present = 1, row_names = paste(samplearea_id_500, year, sep = "_")) %>%
  dplyr::select(-n) %>%
  spread(key = species_id, value = present, fill = 0) %>%
  column_to_rownames(var = "row_names")

Y_500_abund = tm_points %>% 
  filter(year %in% sampling_period & samplearea_id %in% forest_cells & species_id %in% final_species) %>% 
  group_by(species_id, samplearea_id_500, year) %>% 
  tally() %>% 
  mutate(row_names = paste(samplearea_id_500, year, sep = "_")) %>% 
  spread(key = species_id, value = n, fill = 0) %>% 
  column_to_rownames(var = "row_names")

X_500 = climtopo_extr_500 %>%   
  left_join(lidar_extr_500, by = "samplearea_id_500") %>% 
  filter(samplearea_id_500 %in% Y_500_occ$samplearea_id_500) %>% 
  dplyr::select(samplearea_id_500, var_select_final) %>% 
  mutate_at(vars(-samplearea_id_500), list(sq = ~.*.)) %>% # add quadratic terms
  mutate_at(vars(-samplearea_id_500), scale) %>% # center and recale variables
  mutate_at(vars(-samplearea_id_500), as.numeric) %>% # remove attributes from scaled vars
  right_join(Y_500_abund[,c("samplearea_id_500", "year")], by = "samplearea_id_500") %>% 
  mutate(row_names = paste(samplearea_id_500, year, sep = "_")) %>% 
  dplyr::select(-samplearea_id_500, -year) %>% 
  column_to_rownames(var = "row_names")

all(rownames(X_500) == rownames(Y_500_occ))
all(rownames(X_500) == rownames(Y_500_abund))

X_500 = X_500[complete.cases(X_500),]
Y_500_occ = dplyr::select(Y_500_occ, -samplearea_id_500, -year) %>% 
  filter(rownames(.) %in% rownames(X_500))
Y_500_abund = dplyr::select(Y_500_abund, -samplearea_id_500, -year) %>% 
  filter(rownames(.) %in% rownames(X_500))

study_design = str_split(rownames(X_500), pattern = "_", simplify = T) %>% 
  as_tibble() %>% 
  mutate_each(funs = as.integer) %>% 
  set_names(c("samplearea_id_500", "year")) %>% 
  left_join(unique(sample_areas[,c("samplearea_id", "samplearea_id_500")]), by = "samplearea_id_500") %>% 
  dplyr::select(samplearea_id, samplearea_id_500, year) %>% 
  mutate_each(funs = as_factor)
random_levels = list("samplearea_id" = HmscRandomLevel(units = unique(study_design$samplearea_id)),
                     "samplearea_id_500" = HmscRandomLevel(units = unique(study_design$samplearea_id_500)),
                     "year" = HmscRandomLevel(units = unique(study_design$year)))

model_occ_500 = Hmsc(Y = as.matrix(Y_500_occ), XData = X_500, XFormula = ~., distr = "probit",
                     studyDesign = as.data.frame(study_design), ranLevels = random_levels) 
model_abund_500 = Hmsc(Y = as.matrix(Y_500_abund), XData = X_500, XFormula = ~., distr = "lognormal poisson",
                       studyDesign = as.data.frame(study_design), ranLevels = random_levels) 

##########
# 3. 250m 
Y_250_occ = tm_points %>%
  filter(year %in% sampling_period & samplearea_id %in% forest_cells & species_id %in% final_species) %>%
  group_by(species_id, samplearea_id_250, year) %>%
  tally() %>%
  mutate(present = 1, row_names = paste(samplearea_id_250, year, sep = "_")) %>%
  dplyr::select(-n) %>%
  spread(key = species_id, value = present, fill = 0) %>%
  column_to_rownames(var = "row_names")

Y_250_abund = tm_points %>% 
  filter(year %in% sampling_period & samplearea_id %in% forest_cells & species_id %in% final_species) %>% 
  group_by(species_id, samplearea_id_250, year) %>% 
  tally() %>% 
  mutate(row_names = paste(samplearea_id_250, year, sep = "_")) %>% 
  spread(key = species_id, value = n, fill = 0) %>% 
  column_to_rownames(var = "row_names")

X_250 = climtopo_extr_250 %>%   
  left_join(lidar_extr_250, by = "samplearea_id_250") %>% 
  filter(samplearea_id_250 %in% Y_250_occ$samplearea_id_250) %>% 
  dplyr::select(samplearea_id_250, var_select_final) %>% 
  mutate_at(vars(-samplearea_id_250), list(sq = ~.*.)) %>% # add quadratic terms
  mutate_at(vars(-samplearea_id_250), scale) %>% # center and recale variables
  mutate_at(vars(-samplearea_id_250), as.numeric) %>% # remove attributes from scaled vars
  right_join(Y_250_abund[,c("samplearea_id_250", "year")], by = "samplearea_id_250") %>% 
  mutate(row_names = paste(samplearea_id_250, year, sep = "_")) %>% 
  dplyr::select(-samplearea_id_250, -year) %>% 
  column_to_rownames(var = "row_names")

all(rownames(X_250) == rownames(Y_250_occ))
all(rownames(X_250) == rownames(Y_250_abund))

X_250 = X_250[complete.cases(X_250),]
Y_250_occ = dplyr::select(Y_250_occ, -samplearea_id_250, -year) %>% 
  filter(rownames(.) %in% rownames(X_250))
Y_250_abund = dplyr::select(Y_250_abund, -samplearea_id_250, -year) %>% 
  filter(rownames(.) %in% rownames(X_250))

study_design = str_split(rownames(X_250), pattern = "_", simplify = T) %>% 
  as_tibble() %>% 
  mutate_each(funs = as.integer) %>% 
  set_names(c("samplearea_id_250", "year")) %>% 
  left_join(unique(sample_areas[,c("samplearea_id", "samplearea_id_250")]), by = "samplearea_id_250") %>% 
  dplyr::select(samplearea_id, samplearea_id_250, year) %>% 
  mutate_each(funs = as_factor)
random_levels = list("samplearea_id" = HmscRandomLevel(units = unique(study_design$samplearea_id)),
                     "samplearea_id_250" = HmscRandomLevel(units = unique(study_design$samplearea_id_250)),
                     "year" = HmscRandomLevel(units = unique(study_design$year)))

model_occ_250 = Hmsc(Y = as.matrix(Y_250_occ), XData = X_250, XFormula = ~., distr = "probit",
                     studyDesign = as.data.frame(study_design), ranLevels = random_levels) 
model_abund_250 = Hmsc(Y = as.matrix(Y_250_abund), XData = X_250, XFormula = ~., distr = "lognormal poisson",
                       studyDesign = as.data.frame(study_design), ranLevels = random_levels) 

#########
# 4. 125m
Y_125_occ = tm_points %>%
  filter(year %in% sampling_period & samplearea_id %in% forest_cells & species_id %in% final_species) %>%
  group_by(species_id, samplearea_id_125, year) %>%
  tally() %>%
  mutate(present = 1, row_names = paste(samplearea_id_125, year, sep = "_")) %>%
  dplyr::select(-n) %>%
  spread(key = species_id, value = present, fill = 0) %>%
  column_to_rownames(var = "row_names")

Y_125_abund = tm_points %>% 
  filter(year %in% sampling_period & samplearea_id %in% forest_cells & species_id %in% final_species) %>% 
  group_by(species_id, samplearea_id_125, year) %>% 
  tally() %>% 
  mutate(row_names = paste(samplearea_id_125, year, sep = "_")) %>% 
  spread(key = species_id, value = n, fill = 0) %>% 
  column_to_rownames(var = "row_names")

X_125 = climtopo_extr_125 %>%   
  left_join(lidar_extr_125, by = "samplearea_id_125") %>% 
  filter(samplearea_id_125 %in% Y_125_occ$samplearea_id_125) %>% 
  dplyr::select(samplearea_id_125, var_select_final) %>% 
  mutate_at(vars(-samplearea_id_125), list(sq = ~.*.)) %>% # add quadratic terms
  mutate_at(vars(-samplearea_id_125), scale) %>% # center and recale variables
  mutate_at(vars(-samplearea_id_125), as.numeric) %>% # remove attributes from scaled vars
  right_join(Y_125_abund[,c("samplearea_id_125", "year")], by = "samplearea_id_125") %>% 
  mutate(row_names = paste(samplearea_id_125, year, sep = "_")) %>% 
  dplyr::select(-samplearea_id_125, -year) %>% 
  column_to_rownames(var = "row_names")

all(rownames(X_125) == rownames(Y_125_occ))
all(rownames(X_125) == rownames(Y_125_abund))

X_125 = X_125[complete.cases(X_125),]
Y_125_occ = dplyr::select(Y_125_occ, -samplearea_id_125, -year) %>% 
  filter(rownames(.) %in% rownames(X_125))
Y_125_abund = dplyr::select(Y_125_abund, -samplearea_id_125, -year) %>% 
  filter(rownames(.) %in% rownames(X_125))

study_design = str_split(rownames(X_125), pattern = "_", simplify = T) %>% 
  as_tibble() %>% 
  mutate_each(funs = as.integer) %>% 
  set_names(c("samplearea_id_125", "year")) %>% 
  left_join(unique(sample_areas[,c("samplearea_id", "samplearea_id_125")]), by = "samplearea_id_125") %>% 
  dplyr::select(samplearea_id, samplearea_id_125, year) %>% 
  mutate_each(funs = as_factor)
random_levels = list("samplearea_id" = HmscRandomLevel(units = unique(study_design$samplearea_id)),
                     "samplearea_id_125" = HmscRandomLevel(units = unique(study_design$samplearea_id_125)),
                     "year" = HmscRandomLevel(units = unique(study_design$year)))

model_occ_125 = Hmsc(Y = as.matrix(Y_125_occ), XData = X_125, XFormula = ~., distr = "probit",
                     studyDesign = as.data.frame(study_design), ranLevels = random_levels) 
model_abund_125 = Hmsc(Y = as.matrix(Y_125_abund), XData = X_125, XFormula = ~., distr = "lognormal poisson",
                       studyDesign = as.data.frame(study_design), ranLevels = random_levels) 

####################################################
#           Fit 5-year models on cluster           #
####################################################
fit_model = function(model_name){
  tmp_model = model_list[[model_name]]
  tmp_model_fit = Hmsc::sampleMcmc(tmp_model, nChains = 8, thin = 50, samples = 500, transient = 5000, nParallel = 8)
  return(tmp_model_fit)
}
model_list = list("model_occ_1000" = model_occ_1000, "model_abund_1000" = model_abund_1000, 
                  "model_occ_500" = model_occ_500, "model_abund_500" = model_abund_500, 
                  "model_occ_250" = model_occ_250, "model_abund_250" = model_abund_250, 
                  "model_occ_125" = model_occ_125, "model_abund_125" = model_abund_125)
pars = data.frame(model_name = names(model_list), stringsAsFactors = F)

# CAUTION: This runs up to 20 days!
models_fit_5 <- slurm_apply(fit_model, # function to be parallelized
                            params = pars, # function parameters (one column = one parameter)
                            jobname = "fit_models_5years",
                            add_objects = c("model_list"),
                            slurm_options = list(qos = "long", partition = "computehm"),
                            nodes = 8, cpus_per_node = 1, # Requested resources 
                            submit = F) # Test run or submit?
# Submit = F because rslurm requests one core per row in 'pars', but we need 8 cores per row
# --> Change line "#SBATCH --cpus-per-task=1" to "#SBATCH --cpus-per-task=8" in submit.sh
# --> use terminal and navigate to newly created rslurm folder
# --> run command "sbatch submit.sh"

print_job_status(models_fit_5)

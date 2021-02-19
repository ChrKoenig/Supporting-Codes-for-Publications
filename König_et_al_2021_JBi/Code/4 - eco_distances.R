library(tidyverse)
library(StatMatch)
library(corrplot)
library(fastDummies)

rm(list = ls())
setwd("") # Set working directory here

load("Data/tm_points.RData")
load("Data/final_species.RData")
load("Data/ID_lookup.RData")

# Prepare species trait data
traits_raw = read.delim("Data/Bird_traits_Storchova_2018.txt", header = T, stringsAsFactors = F) %>% 
  mutate(Species = case_when(Species == 'Poecile palustris' ~ 'Parus palustris',
                             Species == 'Cyanistes caeruleus' ~ 'Parus caeruleus',
                             Species == 'Periparus ater' ~ 'Parus ater',
                             Species == 'Lophophanes cristatus' ~ 'Parus cristatus',
                             Species == 'Poecile montanus' ~ 'Parus montanus',
                             Species == 'Regulus ignicapilla' ~ 'Regulus ignicapillus',
                             Species == 'Spinus spinus' ~ 'Carduelis spinus',
                             Species == 'Carduelis citrinella' ~ 'Serinus citrinella',
                             Species == 'Acanthis flammea' ~ 'Carduelis flammea',
                             Species == 'Lyrurus tetrix' ~ 'Tetrao tetrix',
                             TRUE ~ Species))

traits_final = ID_lookup %>% 
  left_join(traits_raw, by = c("latin" = "Species")) %>% 
  dplyr::filter(species_id %in% final_species) %>% 
  dplyr::select(species_id, bodymass = WeightU_MEAN, clutch_size = Clutch_MEAN, breeding_age = Age.of.first.breeding, 
                development_type = Young, nest_type = Nest.type, nesting_behaviour = Association.during.nesting, sedentary = Sedentary,
                Folivore_B:Omnivore_B) %>% 
  dummy_cols(select_columns = c("development_type", "nest_type", "nesting_behaviour", "sedentary"), remove_selected_columns = T) %>% 
  column_to_rownames("species_id")

variable_weights = c(1, # bodymass
                     1, # clutch_size
                     1, # breeding age
                     rep(1/9, 9), # Food type
                     rep(1/3, 3), # developmental type
                     rep(1/5, 5), # nest type
                     rep(1/4, 4), # nesting behaviour
                     rep(1/2, 2)) # sedentary

f_sim = 1 - gower.dist(traits_final, var.weights = variable_weights)
colnames(f_sim) = rownames(f_sim) = rownames(traits_final)
corrplot(f_sim, is.corr = F, order = "hclust")
save(f_sim, file = "Data/f_sim.Rdata")

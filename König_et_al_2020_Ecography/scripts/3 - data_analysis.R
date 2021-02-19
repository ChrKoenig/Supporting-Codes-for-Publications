library(brms)
library(car)
library(tidyverse)
library(broom)
library(dotwhisker)
library(grid)
library(gridExtra)

setwd(".") # Set WD here

load("../Data/raw_data/geoentities_meta.RData")
load("../Data/raw_data/geoentities.RData")
load("../Data/raw_data/checklists.RData")
load("../Data/raw_data/family_traits_final.RData")
load("../Data/D_family.RData")

########################################
### PLOT DATA (ADJUST TO YOUR NEEDS) ###
########################################
source("./plotting_functions.R")

# Plot dataset
png("../Plots/geoentities.png", 1600, 800, res = 110)
plot_geoentities()
dev.off()

# Plot sourcepools
# Hawaii
png("../Plots/Hawaii_source_pool.png", 1600, 600, res = 110) 
  plot_sourcepool(171, plot_histogram = T) 
dev.off()
png("../Plots/Hawaii_source_pool_agg.png", 1600, 600, res = 110) 
  plot_sourcepool_agg(171, plot_histogram = T) # Aggregated weights of available mainland floras
dev.off()

# La Reunion
png("../Plots/La_Reunion_source_pool.png", 1600, 600, res = 110) 
  plot_sourcepool(202, plot_histogram = T) 
dev.off()
png("../Plots/La_Reunion_source_pool_agg.png", 1600, 600, res = 110) 
  plot_sourcepool_agg(202, plot_histogram = T) # Aggregated weights of available mainland floras
dev.off()

# Tenerife
png("../Plots/Tenerife_source_pool.png", 1600, 600, res = 110) 
  plot_sourcepool(148, plot_histogram = T) 
dev.off()
png("../Plots/Tenerife_source_pool_agg.png", 1600, 600, res = 110) 
  plot_sourcepool_agg(148, plot_histogram = T) # Aggregated weights of available mainland floras
dev.off()

# Plot representation of families
plot_representational_disharmony("Brassicaceae")
family_disharmony = plot_representational_disharmony(family_names = c("Asteraceae", "Orchidaceae", "Polypodiaceae",
                                                                 "Primulaceae", "Amaranthaceae", "Pinaceae"), ncol = 2) 
ggsave("../Plots/representational_disharmony.pdf", family_disharmony, "pdf", width = 16.8, height = 10.5)
ggsave("../Plots/representational_disharmony.png", family_disharmony, "png", width = 16.8, height = 10.5, dpi = 300)

# Plot aggregated disharmony per island
island_disharmony = plot_compositional_disharmony()
ggsave("../Plots/compositional_disharmony.png", island_disharmony, width = 16.8, height = 6.7, dpi = 300, device = "png")

maps = ggpubr::ggarrange(family_disharmony, island_disharmony, ncol = 1, labels = c("A", "B"), widths = c(16.8, 16.8), heights = c(10.5,6.7), font.label = list(size = 20))
ggsave("../Plots/Fig3.png", maps, "png", width = 16.8, height = 10.5+6.7, dpi = 300)

# Plot trait effects
trait_plots = plot_trait_effects()
ggsave("../Plots/trait_plots.pdf", trait_plots, "pdf", width = 10, height = 10)
ggsave("../Plots/trait_plots.png", trait_plots, "png", width = 10, height = 10)

####################
### Analyze data ###
####################
# 1. Representational Disharmony
# Prepare data frame
family_df_model = D_family %>% 
  filter(!is.na(D)) %>% 
  left_join(family_traits_final %>% dplyr::select(-source) %>% spread(key = trait_name, value = trait_value), by = "family") %>% 
  mutate(woodiness = factor(woodiness, levels = c("non-woody", "woody"))) %>% 
  mutate(pollination_syndrome = factor(pollination_syndrome, levels = c("abiotic", "biotic"))) %>% 
  mutate(dispersal_syndrome = factor(dispersal_syndrome, levels = c("unspecialized", "autochorous", "anemochorous", "hydrochorous", "zoochorous"))) %>% 
  mutate_at(c("seed_mass", "plant_height", "sla"), list(~scale(log10(as.numeric(.))))) %>% 
  mutate(D = D * 0.99999999999) # Beta distribution only defined for 0 < x < 1 
  
# Fit models
brm_full = brm(D ~ woodiness + dispersal_syndrome + pollination_syndrome + plant_height + seed_mass + sla, family = Beta(), data = family_df_model)

# Inspect fit
(tidy(brm_full))
plot(brm_full)
summary(brm_full)
bayes_R2(brm_full)

# Plot effects
effect_size_df_families = tidy(brm_full) %>% 
  mutate(term = dplyr::recode(term, b_woodinesswoody = "woody", b_dispersal_syndromeautochorous = "autochorous", b_dispersal_syndromeanemochorous = "anemochorous", 
                              b_dispersal_syndromehydrochorous = "hydrochorous", b_dispersal_syndromezoochorous = "zoochorous", b_pollination_syndromebiotic = "biotic",
                              b_plant_height = "plant height", b_seed_mass = "seed mass", b_sla = "SLA")) %>% 
  filter(!term %in% c("b_Intercept", "phi", "lp__")) %>% 
  add_row(term = "non-woody", estimate = 0, std.error = 0.00000001, lower = 0.00000001, upper = 0.00000001) %>% 
  add_row(term = "unspecialized", estimate = 0, std.error = 0.00000001, lower = 0.00000001, upper = 0.00000001) %>% 
  add_row(term = "abiotic", estimate = 0, std.error = 0.00000001, lower = 0.00000001, upper = 0.00000001) %>% 
  mutate(variable = dplyr::recode(term, `non-woody` = "Woodiness", woody = "Woodiness", unspecialized = "Dispersal syndrome", autochorous = "Dispersal syndrome", 
                               anemochorous = "Dispersal syndrome", hydrochorous = "Dispersal syndrome", zoochorous = "Dispersal syndrome", 
                               abiotic = "Pollination syndrome", biotic = "Pollination syndrome", 
                               `plant height` = "Maximum plant height", `seed mass` = "Seed mass", `SLA` = "Specific Leaf Area"),
         term = factor(term, levels = rev(c("non-woody", "woody", "abiotic", "biotic", "unspecialized", "autochorous", "anemochorous", "hydrochorous", "zoochorous", 
                                            "plant height", "seed mass", "SLA"))))

seed_mass = plot_model_coefs(effect_size_df_families, "Seed mass")
plant_height = plot_model_coefs(effect_size_df_families, "Maximum plant height")
sla = plot_model_coefs(effect_size_df_families, "Specific Leaf Area") 
woodiness = plot_model_coefs(effect_size_df_families, "Woodiness", show_title = T)
pollination = plot_model_coefs(effect_size_df_families, "Pollination syndrome", show_title = T)
dispersal = plot_model_coefs(effect_size_df_families, "Dispersal syndrome", show_title = T, show_x_axis = T)

grobz1 = lapply(list(seed_mass, plant_height, sla, woodiness, pollination, dispersal), ggplotGrob)
grobz1 = rbind(grobz1[[1]], grobz1[[2]], grobz1[[3]], grobz1[[4]], grobz1[[5]], grobz1[[6]],  size = "first")

#png("../Plots/effect_sizes_fam.png", width = 12, height = 16, units = "cm", res = 300)
pdf("../Plots/effect_sizes_fam.pdf", width = 5, height = 5*1.25)
grid.newpage()
grid.draw(resize_heights(grobz1, c(1,1,1,2,2,5)))
dev.off()

# 2. Compositional disharmony
island_df_model =  D_family %>% 
  mutate(deviance = abs(D-0.5)) %>% 
  group_by(entity_ID) %>% 
  dplyr::summarize(D = median(deviance, na.rm = T)) %>% 
  inner_join(geoentities_meta, by = "entity_ID") %>% 
  mutate(area = scale(log10(area)), elev = scale(elev), dist = scale(dist), geology = recode_factor(geology, `6` = "volcanic", `1` = "atoll", `3` = "uplift"), 
         archipelago = as.factor(ifelse(is.na(arch_lvl_3), ifelse(is.na(arch_lvl_2), arch_lvl_1, arch_lvl_2), arch_lvl_3))) %>% 
  dplyr::select(entity_ID, geo_entity, longitude, latitude, area, dist, elev, geology, archipelago, D)

brm_islands_full = brm(D ~ area + dist + elev + geology, data = island_df_model)

# Inspect fit
(tidy(brm_islands_full))
plot(brm_islands_full)
summary(brm_islands_full)
bayes_R2(brm_islands_full)

# Plot effects
effect_size_df_islands = tidy(brm_islands_full) %>% 
  add_row(term = "Geology: volcanic", estimate = 0, std.error = 0.00000001, lower = 0.00000001, upper = 0.00000001) %>% 
  filter(!term %in% c("b_Intercept", "phi", "lp__", "sigma")) %>% 
  mutate(term = dplyr::recode(term, `b_area` = "Area", `b_dist` = "Distance to mainland", `b_elev` = "Elevation", `b_geologyatoll` = "Geology: atoll", `b_geologyuplift` = "Geology: uplift"),
         variable = c("Area", "Distance to mainland", "Elevation", "Geology", "Geology", "Geology"))

area = plot_model_coefs(effect_size_df_islands, "Area", limits = c(-0.04,0.04), breaks = seq(-0.04,0.04, 0.01))
dist = plot_model_coefs(effect_size_df_islands, "Distance to mainland", limits = c(-0.04,0.04), breaks = seq(-0.04,0.04, 0.01))
elev = plot_model_coefs(effect_size_df_islands, "Elevation",limits = c(-0.04,0.04), breaks = seq(-0.04,0.04, 0.01))
geology = plot_model_coefs(effect_size_df_islands, "Geology", limits = c(-0.04,0.04), breaks = seq(-0.04,0.04, 0.01), show_x_axis = T, show_title = T)

grobz2 = lapply(list(area, dist, elev, geology), ggplotGrob)
grobz2 = rbind(grobz2[[1]], grobz2[[2]], grobz2[[3]], grobz2[[4]],  size = "last")

pdf("../Plots/effect_sizes_isl.pdf", width = 5, height = 5 * 0.71666)
grid.newpage()
grid.draw(resize_heights(grobz2, c(1,1,1,3)))
dev.off()

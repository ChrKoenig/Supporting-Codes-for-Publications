library(tidyverse)
library(Hmsc)
library(vegan)

rm(list = ls())
setwd("") # Set working directory here

load("Data/tm_points.RData")
load("Data/sample_areas.RData")
load("Data/forest_cells.RData")
load("Data/pred_performance.RData")
load("Data/ID_lookup.RData")
load("Models/Model_runs_7/jsdm_fit_1.RData")
load("Models/Model_runs_7/jsdm_fit_5.RData")

#######################
# Model performance
model_performance = read_csv("../TableS2_model_performance.csv", 
                             col_types = c(species_name = "f", grain_size = "f", sampling_period = "f", 
                                           metric = "f", pred_type = "f", value = "n"))

model_performance_plot = ggplot(model_performance, aes(x = grain_size, y = value, color = pred_type)) +
  geom_jitter(position = position_jitterdodge(), alpha = 0.3) +
  geom_boxplot(outlier.shape = NA,  notch = T, alpha = 0.75, outlier.colour = "black") +
  facet_grid(cols = vars(metric), rows = vars(sampling_period), 
             labeller = labeller(sampling_period = c("1" = "1 year", "5" = "5 years"))) +
  scale_color_manual(labels = c("In-sample", "Out-of-sample"), values = c("#3498db", "#f39c12")) +
  labs(color = "Prediction type", x = "grain size") +
  theme_bw() +
  theme(plot.subtitle = element_text(hjust = 0.5),
        strip.background = element_rect(fill="black"),
        strip.text = element_text(colour = 'white', face = "bold"))
ggsave("./Figures/model_performance.png", plot = model_performance_plot, device = "png",  width = 16.4, height = 10, units = c("cm"), dpi = 300)

# Check relationship between model fit and prevalence
prevalence = tm_points %>% 
  filter(samplearea_id %in% forest_cells & year %in% 2007:2011) %>% 
  mutate(species_id = as.character(species_id)) %>% 
  group_by(species_id) %>% 
  summarise(count = n())

prevalence_fit_df = left_join(model_fits, prevalence)

model_fits_prevalence_plot = ggplot(filter(prevalence_fit_df, metric %in% c("AUC", "TjurR2")), aes(x = count, y = value, col = metric)) +
  geom_point() +
  stat_smooth(method = "glm", method.args = list(family = "binomial")) +
  scale_x_continuous(name = "Total observations (2007-2011)", trans = "log10") +
  labs(col = "Performance metric:") +
  facet_grid(rows = vars(temp_extent), cols = vars(grain_size)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill="black"),
        strip.text = element_text(colour = 'white', face = "bold"))
ggsave("./Figures/model_fits_vs_prevalence.png", plot = model_fits_prevalence_plot , width = 16.4, height = 11, units = c("cm"), dpi = 300)

#######################
# Variance partitioning
extract_var_part = function(x){
  object_name = gsub(".RDS", "", x)
  species_ids = jsdm_model_results_1$model_occ_1000$spNames # constant across models, so always take this one
  grain_size = gsub("(^.*occ_)([0-9]+)(_varPart.*$)", replacement = "\\2", object_name)
  temp_extent = gsub("(^.*)([0-9]$)", replacement = "\\2", object_name)
  var_part = readRDS(paste0("Data/varpart/", x))
  
  var_part_df = as.data.frame(t(var_part$vals)) %>% 
    rownames_to_column("species_id") %>% 
    mutate(grain_size = grain_size, temp_extent = temp_extent) %>% 
    rename_at(vars(ends_with(c("0","5"))), ~ paste0("sample_lvl")) 
}

var_parts = bind_rows(lapply(list.files("Data/varpart/", pattern = "varPart"), extract_var_part)) %>% 
  left_join(ID_lookup) %>% 
  select(species_name = latin, grain_size, sampling_period = temp_extent, fixed, plot_lvl = `Random: samplearea_id`, sample_lvl, year = `Random: year`) %>% 
  mutate(sample_lvl = ifelse(grain_size == "1000", plot_lvl, sample_lvl),
         plot_lvl = ifelse(grain_size == "1000", NA, plot_lvl)) %>% 
  replace_na(list(year = 0, sample_lvl = 0, plot_lvl = 0)) %>% 
  left_join(filter(model_performance, metric == "TjurR2"), by = c("species_name", "grain_size", "sampling_period")) %>% 
  select(species_name:sampling_period, total = fit, fixed, sample_lvl, plot_lvl, year) %>% 
  pivot_longer(cols = `total`:`year`, names_to = "term") %>% 
  mutate(grain_size = factor(grain_size, levels = c("125", "250", "500", "1000"))) %>% 
  remove_missing()

write.csv(var_parts, "Manuscript/Supporting Information/TableS3_variance_partitioning.csv", na = "",row.names = F)

var_parts_summary = var_parts %>% 
  group_by(grain_size, sampling_period, term) %>% 
  summarize(median_val = round(median(value), 2), 
            quant05 = round(quantile(value, .05), 2),
            quant95 = round(quantile(value, .95), 2)) %>% 
  mutate(val = paste0(median_val, " (", quant05, "-", quant95, ")")) %>% 
  pivot_wider(id_cols = c(grain_size, sampling_period), names_from = term, values_from = val) %>% 
  select(grain_size, sampling_period, total, fixed, sample_lvl, plot_lvl, year) %>% 
  arrange(grain_size, sampling_period)

ggplot(filter(var_parts, temp_extent == 5 & term == "sample_lvl"), aes(x = as.numeric(grain_size), y = value, group = species_id, color = species_id)) + 
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), se = F) +
  scale_x_continuous(breaks = 1:4, labels = c("125", "250", "500", "1000")) +
  ylab("Partial R²: Species associations") +
  xlab("Grain size (m)") +
  theme_bw() +
  theme(legend.position = "none")

##########################################################################
# Species associations
get_species_associations = function(x, y){
  # x: model, y: model name
  grain_size = paste0(parse_number(y), "m")
  data_type = ifelse(grepl("occ", y), "occurrence", "prevalence")
  vcv_all = computeAssociations(x) # species association matrices with all random effects 
  vcv_sampleUnit = vcv_all[[ifelse(grain_size == "1000m", 1, 2)]] # choose the correct random effect (at the level of the sampling unit)
  #species = factor(rownames(vcv_sampleUnit$mean))
  data.frame(spec_1 = rownames(vcv_sampleUnit$mean)[row(vcv_sampleUnit$mean)[upper.tri(vcv_sampleUnit$mean, diag = T)]], 
             spec_2 = colnames(vcv_sampleUnit$mean)[col(vcv_sampleUnit$mean)[upper.tri(vcv_sampleUnit$mean, diag = T)]], 
             grain_size = factor(grain_size, levels = c("125m", "250m", "500m", "1000m")),
             data_type = factor(data_type, levels = c("occurrence", "prevalence")),
             estimate = vcv_sampleUnit$mean[upper.tri(vcv_sampleUnit$mean, diag = T)],
             support = vcv_sampleUnit$support[upper.tri(vcv_sampleUnit$support, diag = T)],
             stringsAsFactors = F)
}

# Produce flat VCV matrices
vcv_flat_1 = imap_dfr(jsdm_model_results_1, .f = get_species_associations) %>% 
  mutate(temp_extent = "1 year") 
vcv_flat_5 = imap_dfr(jsdm_model_results_5, .f = get_species_associations) %>% 
  mutate(temp_extent = "5 years")

vcv_flat_all = bind_rows(vcv_flat_1, vcv_flat_5) %>% 
  filter(data_type == "occurrence") %>% # prevalence models had convergence problems
  drop_na() %>% 
  mutate(sgnf = (support > 0.975 | support < 0.025))
save(vcv_flat_all, file = "Data/vcv_flat_all.RData")

# Order species based on clustering (use 5 year, 125m model)
reference = computeAssociations(jsdm_model_results_5$model_occ_250)[2][[1]][[1]] # get sample level VCV matrix
plot_order = colnames(reference)[hclust(d = dist(reference), method = "average")$order]

vcv_plot = vcv_flat_all %>% # reorder spec1-spec2
  mutate(spec_1 = fct_relevel(spec_1, plot_order),
         spec_2 = fct_relevel(spec_2, plot_order)) %>% 
  mutate(spec_1_ind = ifelse(as.numeric(spec_1) < as.numeric(spec_2), spec_1, spec_2),
         spec_2_ind = ifelse(as.numeric(spec_1) < as.numeric(spec_2), spec_2, spec_1),
         spec_1 = factor(levels(spec_1)[spec_1_ind], levels = levels(spec_1)),
         spec_2 = factor(levels(spec_2)[spec_2_ind], levels(spec_2))) %>% 
  select(-spec_1_ind, -spec_2_ind) 

correlation_matrices_plot = ggplot(data = vcv_plot, aes(x = spec_1, y = spec_2, fill = estimate)) +
  geom_tile() +
  geom_tile(data = filter(vcv_plot, sgnf == T), aes(x = spec_1, y = spec_2), col = "black") +
  scale_fill_gradient2(name = "Residual correlation", low = "#C10000", mid = "grey90", high = "#0000C1", limits = c(-1,1),
                       guide = guide_colorbar(direction = "horizontal", title.position = "top", title.hjust = 0.5, barwidth = 20, barheight = .6)) +
  scale_y_discrete(limits = rev(levels(vcv_plot$spec_2))) +
  xlab("Species 1") + ylab("Species 2") +
  facet_grid(cols = vars(grain_size), rows = vars(temp_extent)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.position = "bottom", legend.justification = c(0.53,1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        strip.background = element_rect(fill="black"),
        strip.text = element_text(colour = 'white', face = "bold"))
ggsave("./Figures/correlation_matrices.png", plot = correlation_matrices_plot, width = 16.4, height = 11, units = c("cm"), dpi = 300)

# Consistency of residual correlations
vcv_grain_contrast = vcv_flat_all %>% 
  left_join(filter(vcv_flat_all, grain_size == "125m"), by = c("spec_1", "spec_2", "temp_extent")) %>% 
  mutate(significant = case_when(
    (support.x > 0.975 | support.x < 0.025) & (support.y < 0.975 & support.y > 0.025) ~ "coarse",
    (support.y > 0.975 | support.y < 0.025) & (support.x < 0.975 & support.x > 0.025) ~ "fine",
    (support.x > 0.975 | support.x < 0.025) & (support.y > 0.975 | support.y < 0.025) ~ "both",
    (support.x < 0.975 & support.x > 0.025) & (support.y < 0.975 & support.y > 0.025) ~ "none")) %>% 
  dplyr::select(spec_1, spec_2, grain_size = grain_size.x, temp_extent, estimate_125 = estimate.y, estimate = estimate.x, significant) %>% 
  remove_missing()

grain_contrast_plot = ggplot(arrange(vcv_grain_contrast, desc(significant)), aes(y = estimate_125, x = estimate, color = significant)) +
  geom_point(color = "grey80") +
  geom_point(fill = NA, shape = 21) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  xlim(-1, 1) + xlab("residual correlation") +
  ylim(-1, 1) + ylab("residual correlation (125 m)") +
  facet_grid(cols = vars(grain_size), rows = vars(temp_extent)) +
  scale_color_manual(values = c(fine = "#e6a129", both = "black", coarse = "#bd29e6", none = "grey80"),
                     name = "Significant at grain size:",
                     labels=c("Fine & Coarse", "Coarse", "Fine", "None"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.position = "bottom",
        strip.background = element_rect(fill="black"),
        strip.text = element_text(colour = 'white', face = "bold"))

ggsave("./Figures/residuals_grain_contrast.png", plot = grain_contrast_plot, width = 16.4, height = 10, units = c("cm"), dpi = 300)

# Relationship between species prevalence and residual associations
vcv_prevalence = vcv_flat_all %>% 
  pivot_longer(cols = c(spec_1, spec_2), values_to = "species_id") %>% 
  filter(data_type == "occurrence", sgnf == T) %>% 
  dplyr::select(species_id, grain_size, temp_extent, estimate) %>% 
  group_by(species_id, grain_size, temp_extent) %>% 
  summarize(mean_estimate = mean(estimate)) %>% 
  left_join(prevalence)

vcv_prevalence_plot = ggplot(vcv_prevalence, aes(x = count, y = mean_estimate)) +
  geom_point() +
  stat_smooth(method = "glm") +
  scale_y_continuous("Mean of significant associations", limits = c(-1,1)) +
  scale_x_continuous(name = "Total observations (2007-2011)", trans = "log10") +
  facet_grid(rows = vars(temp_extent), cols = vars(grain_size)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill="black"),
        strip.text = element_text(colour = 'white', face = "bold"))
ggsave("./Figures/associations_vs_prevalence.png", plot = vcv_prevalence_plot, width = 16.4, units = c("cm"), dpi = 300)

##########################################################################
# Look at species pairs with inconsistent associations
vcv_wide_all = vcv_flat_all %>% 
  filter(sgnf == T) %>% 
  pivot_wider(id_cols = c(spec_1, spec_2), names_from = c(grain_size, temp_extent), values_from = c(estimate))

# Negative at fine grain, positive at coarse
switch_lowerright = vcv_wide_all %>% 
  filter(`125m_5 years` < 0 & `1000m_5 years` > 0) %>% 
  left_join(prevalence, by = c("spec_1" = "species_id")) %>% 
  left_join(prevalence, by = c("spec_2" = "species_id")) %>% 
  left_join(ID_lookup[,c("species_id", "latin")], by = c("spec_1" = "species_id")) %>% 
  left_join(ID_lookup[,c("species_id", "latin")], by = c("spec_2" = "species_id")) %>% 
  select(spec_1 =  latin.x, prevalence_1 = count.x, 
         spec_2 = latin.y, prevalence_2 = count.y, 
         omega_1000 = `1000m_5 years`, omega_125 = `125m_5 years`)

# Positive at fine grain, negative at coarse
switch_upperleft = vcv_wide_all %>% 
  filter(`125m_5 years` > 0 & `1000m_5 years` < 0) %>% 
  left_join(prevalence, by = c("spec_1" = "species_id")) %>% 
  left_join(prevalence, by = c("spec_2" = "species_id")) %>% 
  left_join(ID_lookup[,c("species_id", "latin")], by = c("spec_1" = "species_id")) %>% 
  left_join(ID_lookup[,c("species_id", "latin")], by = c("spec_2" = "species_id")) %>% 
  select(spec_1 =  latin.x, prevalence_1 = count.x, 
         spec_2 = latin.y, prevalence_2 = count.y, 
         omega_1000 = `1000m_5 years`, omega_125 = `125m_5 years`)

##############################
# Mantel correlation between species associations and trait similarity
load("Data/f_sim.Rdata")

# Get VCV-Matrices
species_associations = list("vcv_1000_1" = computeAssociations(jsdm_model_results_1$model_occ_1000)[[1]]$mean,
                            "vcv_500_1" = computeAssociations(jsdm_model_results_1$model_occ_500)[[2]]$mean,
                            "vcv_250_1" = computeAssociations(jsdm_model_results_1$model_occ_250)[[2]]$mean,
                            "vcv_125_1" = computeAssociations(jsdm_model_results_1$model_occ_125)[[2]]$mean,
                            "vcv_1000_5" = computeAssociations(jsdm_model_results_5$model_occ_1000)[[2]]$mean,
                            "vcv_500_5" = computeAssociations(jsdm_model_results_5$model_occ_500)[[2]]$mean,
                            "vcv_250_5" = computeAssociations(jsdm_model_results_5$model_occ_250)[[2]]$mean,
                            "vcv_125_5" = computeAssociations(jsdm_model_results_5$model_occ_125)[[2]]$mean)


mantel_correlations = lapply(species_associations, function(x){
  m1 = x
  m2 = f_sim
  sp_names = intersect(colnames(m1), colnames(m2))
  mantel(m1[sp_names, sp_names], m2[sp_names, sp_names])
})

sapply(mantel_correlations, function(x) {x$statistic})

vcv_flat_fsim = vcv_flat_all %>% 
  mutate(f_sim = apply(select(., spec_1, spec_2), MARGIN = 1, function(x){f_sim[x[1], x[2]]})) 

correlations = vcv_flat_fsim %>%
  group_by(temp_extent, grain_size) %>% 
  summarize(R = round(cor(estimate, f_sim), 2))

correlation_plot = ggplot(vcv_flat_fsim, aes(x = estimate, y = f_sim)) +
  geom_point(colour = "gray 80") +
  geom_point(data = filter(vcv_flat_fsim, sgnf == T), shape = 21, colour = "gray 40") +
  geom_smooth(data = vcv_flat_fsim, aes(x = estimate, y = f_sim), method = "lm", col = "white", fill = "red") +
  geom_text(data = correlations, x = 0, y = 0.475, aes(label = paste("R = ", correlations$R)), inherit.aes = F, hjust = 0.5, size = 2) +
  xlab("residual correlation") +
  ylab("functional similarity") +
  xlim(-1,1) +
  facet_grid(rows = vars(temp_extent), cols = vars(grain_size)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.position = "bottom",
        strip.background = element_rect(fill="black"),
        strip.text = element_text(colour = 'white', face = "bold"))

ggsave("./Figures/fsim_correlations.png", plot = correlation_plot , width = 16.4, height = 9, units = c("cm"), dpi = 300)

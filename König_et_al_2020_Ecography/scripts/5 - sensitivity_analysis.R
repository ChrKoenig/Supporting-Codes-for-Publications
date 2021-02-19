library(tidyverse)

############################################
### Sensitivity to sampling completeness ###
############################################
# Create artificial source pool with 100 families and lognormal distribution of species richness 
set.seed(1)
source_pool = rep(paste0("fam_", 1:100), ceiling(rlnorm(100, 3, 3)))

calc_D_sampled = function(focal_family){
  # Calculate disharmony at different levels of sample completeness
  #   n_i: total species richness of focal island
  #   n_ti: species richness of focal taxon on focal island
  #   p_ti: prevalence of focal taxon in the source pool of the focal island
  #   r_tmp: representation in sample relative to true representation in the source pool
  p_ti = length(which(source_pool == focal_family)) / length(source_pool) 
  D_sampled = bind_rows(
    lapply(c(0.1, 0.5, 0.75, 0.9, 1, 1.1, 1.25, 1.5, 1.9), FUN = function(r_tmp){
      D_r = bind_rows(
        lapply(seq(1000, 50, -50), FUN = function(n_i){ # vary sample size from 1000 to 50
          D_n = replicate(100, expr = {
            sample_tmp = sample(source_pool, size = n_i, replace = F)
            n_ti = round(length(which(sample_tmp == focal_family)) * r_tmp, 0)
            n_exp = round(n_i * p_ti, 0)
            D = ifelse(n_exp < 1, NA, pbinom(q = n_ti-1, size = n_i, prob = p_ti) + dbinom(n_ti, size = n_i, prob = p_ti) * 0.5)
          }) 
          D_n = tibble(sample_size = n_i, 
                       D = D_n,
                       r_tmp = as.character(r_tmp))
          return(D_n)
        }))
      return(D_r)
    }))
}

focal_family = "fam_61" # observe most frequent family - 27035 species
D_sampled = calc_D_sampled(focal_family) 
ggplot(D_sampled, aes(x = sample_size, y = D)) +
  geom_point(alpha = 0.05) +
  geom_smooth(method="glm", method.args = list(family="binomial"), se = F) +
  facet_wrap("r_tmp") +
  xlab("Assemblage size") +
  ylab("Representational disharmony") +
  theme_bw()

focal_family = "fam_11" # observe medium-sized family - 1874 species
D_sampled = calc_D_sampled(focal_family)
ggplot(D_sampled, aes(x = sample_size, y = D)) +
  geom_point(alpha = 0.05) +
  geom_smooth(method="glm", method.args = list(family="binomial"), se = F) +
  facet_wrap("r_tmp") +
  xlab("Assemblage size") +
  ylab("Representational disharmony") +
  theme_bw()

focal_family = "fam_20" # observe medium-sized family - 120 species
D_sampled = calc_D_sampled(focal_family)
ggplot(D_sampled, aes(x = sample_size, y = D)) +
  geom_point(alpha = 0.05) +
  geom_smooth(method="glm", method.args = list(family="binomial"), se = F) +
  facet_wrap("r_tmp") +
  xlab("Assemblage size") +
  ylab("Representational disharmony") +
  theme_bw()
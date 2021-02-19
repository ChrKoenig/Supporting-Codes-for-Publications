library(plyr)
library(dplyr)
library(rgdal)
library(vegan)
library(gdm)
library(geosphere)

setwd(".") # Set WD here
load("../Data/raw_data/checklists.RData")
load("../Data/raw_data/geoentities_meta.RData")

###########################
### Calculate distances ###
###########################
# 1. Compositional dissimilarity (Species level)
spec_tab = unclass(table(checklists[,c("entity_ID", "work_ID")])) # look at species composition
spec_dist = betadiver(spec_tab, "sim") #   Turnover / Simpsons Beta diversity
spec_dist = as.matrix(spec_dist) # Convert to matrix

# 2. Geographical distance
geo_dist = round(distm(as.matrix(geoentities_meta[, c("longitude", "latitude")]), fun = distHaversine)/1000,1)+0.1
colnames(geo_dist) = geoentities_meta$entity_ID
rownames(geo_dist) = geoentities_meta$entity_ID

save(list = c("geo_dist", "spec_dist"), file = "../Data/distances.RData")
gc()

############################################
###  Generalized dissimilarity modeling  ###
############################################
# We fit a model to reflect the default behaviour of floristic similarity on the mainland 
# using distance + climatic variables. Given new data, we can predict the expected number of shared species
# relative to a global equal area grid in order to localize potential source pools for islands. 
spec_dist = cbind("entity_ID" = as.numeric(rownames(spec_dist)), spec_dist)
geo_dist = cbind("entity_ID" = as.numeric(rownames(geo_dist)), geo_dist)
mainlands = geoentities_meta$entity_ID[geoentities_meta$entity_class %in% c("Mainland", "Island/Mainland")]
islands = geoentities_meta$entity_ID[geoentities_meta$entity_class == "Island"]

### 1. Model fitting
spec_dist_ml = spec_dist[paste(mainlands), c("entity_ID", paste(mainlands))]
geo_dist_ml = geo_dist[paste(mainlands), c("entity_ID", paste(mainlands))]
env_data_ml = geoentities_meta %>% filter(entity_ID %in% mainlands) %>% dplyr::select(entity_ID, longitude, latitude, T_mean, T_var, P_mean, P_var)
gdm_tab_ml = formatsitepair(bioData = spec_dist_ml, bioFormat = 3, siteColumn = "entity_ID", abundance = F, XColumn = "longitude", YColumn = "latitude", predData = env_data_ml, distPreds = list("Geo_Distance" = geo_dist_ml))
gdm_fit_ml = gdm(gdm_tab_ml, geo = F)
summary(gdm_fit_ml)

save(gdm_fit_ml, file = "../Data/gdm_fit_ml.RData")

### 2. Prepare prediction Data
# load global equal area grid and prepared environmental variables for each grid cell
load(file = "../Data/raw_data/grid.RData")
load(file = "../Data/raw_data/env_grid.RData")

# Combine environmental data of grid cells and islands 
env_data_is = geoentities_meta %>% filter(entity_ID %in% islands) %>% dplyr::select(entity_ID, longitude, latitude, T_mean, T_var, P_mean, P_var)
env_grid = rbind(env_data_is, env_grid)

# Prepare geographic distance matrix for gdm function 
geo_dist_pred <- sapply(1:nrow(env_grid), function(i){
  distGeo(as.matrix(env_grid[i ,c("longitude","latitude")]), env_grid[,c("longitude","latitude")])
})/1000
geo_dist_pred = cbind(env_grid$entity_ID, geo_dist_pred)
colnames(geo_dist_pred) = c("entity_ID", env_grid$entity_ID)

# Prepare empty compositional dissimilarity matrix (gdm::formatsitepair requires this)
spec_dist_pred = cbind(env_grid$entity_ID, matrix(0, nrow = nrow(env_grid), ncol = nrow(env_grid)))
colnames(spec_dist_pred) = c("entity_ID", env_grid$entity_ID)

# Prepare gdm table
gdm_tab_pred = formatsitepair(bioData = spec_dist_pred, bioFormat = 3, siteColumn = "entity_ID", abundance = F,
                              XColumn = "longitude", YColumn = "latitude", predData = env_grid, 
                              distPreds = list("Geo_Distance" = geo_dist_pred))
save(gdm_tab_pred, file = "../Data/gdm_tab_pred.RData")

### 3. Predict
gdm_pred = predict.gdm(gdm_fit_ml, gdm_tab_pred, geo = F, time = F)
gdm_pred_mat = matrix(NA, nrow = nrow(env_grid), ncol = nrow(env_grid))
gdm_pred_mat[which(lower.tri(gdm_pred_mat))] = gdm_pred
gdm_pred_mat = as.matrix(as.dist(gdm_pred_mat))
entity_names = sort(env_grid$entity_ID)
colnames(gdm_pred_mat) = entity_names
rownames(gdm_pred_mat) = entity_names

save(gdm_pred_mat, file = "../Data/gdm_pred_mat.RData")

#############################
### Aggregate predictions ###
#############ä###############
load("../Data/raw_data/geoentities.RData")

# Since the predictions are for grid cells, not the actual mainland units, we have to aggregate the predicted values
# Get cell IDs per mainland entity
mainland_gridcells = lapply(mainlands, function(ml_ID){
  print(ml_ID)
  tmp_entity = geoentities[which(geoentities@data$entity_ID == ml_ID),]
  tryCatch({
    tmp_grid = raster::intersect(grid, tmp_entity)
    list(ID = tmp_grid@data$ID,
         area = area(tmp_grid)/1000000)},
    error = function(e){NA})
})

names(mainland_gridcells) = mainlands

# Prepare source matirx
source_matrix = matrix(NA, nrow = length(islands), ncol = length(mainlands), byrow = T)
rownames(source_matrix) = islands
colnames(source_matrix) = mainlands

# Fill source matrix
for(island in islands){
  print(island)
  island_richness = length(which(checklists$entity_ID == island))
  tmp_sources = sapply(mainlands, function(x){
    tmp_grid = mainland_gridcells[[paste(x)]]
    if(is.na(tmp_grid) || length(tmp_grid) == 0){return(0)}
    grid_subset = which(paste("cell_", tmp_grid$ID, sep = "") %in%  colnames(gdm_pred_mat))
    grid_IDs = tmp_grid$ID[grid_subset]
    grid_areas = tmp_grid$area[grid_subset]
    if(is.na(grid_IDs) || length(grid_IDs) == 0){return(0)}
    grid_index = which(colnames(gdm_pred_mat) %in% paste("cell_", grid_IDs, sep = ""))
    grid_index = grid_index[which(!is.na(grid_index))]
    grid_values = gdm_pred_mat[paste(island), grid_index]
    p = 1 - weighted.mean(grid_values, grid_areas)
  })
  names(tmp_sources) = mainlands
  
  # Remove all mainland units beyond "knee" (see Supporting Text 2 for justification)
  y = 1- sort(tmp_sources, decreasing = T)
  y_norm = (y-min(y)) / (max(y)-min(y)) # Normalize y to (0,1) to make slopes comparable
  spline_fit = smooth.spline(1:length(y), y_norm, spar = 0.75) # Best for this purpose, because smooth and cont. increasing
  spline_prime = diff(spline_fit$y) # Derivative, i.e. what is the slope of 'spline_fit'?
  spline_runlength = rle(spline_prime < 0.001) # Consecutive slope values smaller/greater than 0.001
  if(tail(spline_runlength$values, n = 1) == T){ # Is the last run < 0.001, i.e. the function flattens out?
    cutoff = length(y) - tail(spline_runlength$lengths, n = 1) # cutoff value
  } else {
    cutoff = length(y) - sum(tail(spline_runlength$lengths, n = 2))
  }
  tmp_sources[names(y[cutoff:length(y)])] = 0 # set all mainlands beyond cutoff to zero
  source_matrix[paste(island),] = tmp_sources # write in matrix
}

save(source_matrix, file = "../Data/source_matrix.RData")

library(rgeos)
library(spdep)
library(raster)
library(maptools)
library(scales)
library(tidyverse)
library(RColorBrewer)
library(ape)

#########################################################
### Plotting functions for visualizing source regions ###
#########################################################
plot_geoentities = function(geoentities, area_threshold = 25000, proj_string = "+proj=longlat +datum=WGS84 +ellps=WGS84"){
  
  load("../Data/raw_data/geoentities_meta.RData")
  load("../Data/raw_data/continents.RData")
  load("../Data/raw_data/geoentities.RData")
  
  continents = raster::crop(continents, raster::extent(-180, 180, -62, 90))
  continents = spTransform(continents, CRS(proj_string))
  geoentities = spTransform(geoentities, CRS(proj_string))
  
  large_entities = subset(geoentities, area > area_threshold)
  par(mar = c(0,0,0,0))  
  plot(continents, border = "lightgrey", col = "lightgrey")
  for(i in 1:length(large_entities)){
    entity = large_entities[i,]
    print(entity$entity_ID)
    if(entity@data$entity_class == "Island"){
      tryCatch(expr = {
        plot(gSimplify(entity, tol = 0.01), add = T, col= "#0099FF66", border =  "#0099FF66")}, 
        warning = function(w){warning(w)}, 
        error = function(e){plot(entity, add = T, col= "#0099FF66", border =  "#0099FF66")})
    } else {
      tryCatch(expr = {
        plot(gSimplify(entity, tol = 0.01), add = T, col= "#F7595966", border = "#FF354466")},
        warning = function(w){warning(w)}, 
        error = function(e){plot(entity, add = T, col= "#F7595966", border = "#FF354466")})
    }
    small_entities = subset(geoentities@data, area < area_threshold)
    coords = SpatialPoints(small_entities[,c("point_x","point_y")], proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
    coords = spTransform(coords, CRS(proj_string))
    colors =  ifelse(small_entities$entity_class == "Island", "#5465FF44", "#FF322444")
    plot(coords, pch = 20, col = colors, cex = 2, add = T)
  }
}

# Visualize the grid-level source region predictions for a given island
plot_sourcepool = function(entity_ID, plot_histogram = T, main = NULL){
  # entity_ID = focal island
  # plot_histogram = plot a histogram next to the map?
  # main = plot title
  
  load("../Data/raw_data/geoentities_meta.RData")
  load("../Data/raw_data/grid.RData")
  load("../Data/raw_data/continents.RData")
  load("../Data/gdm_pred_mat.RData")
  
  if(plot_histogram){layout(matrix(1:2, ncol = 2), widths = c(10,3))}
  par(mar = c(0,0,0,0))
  if(!is.null(main)){par(mar = c(0,0,3,0))}
  
  tmp_entity_orig = gdm_pred_mat[paste(entity_ID), grep("cell", colnames(gdm_pred_mat))]
  names(tmp_entity_orig) = gsub("cell_", "", names(tmp_entity_orig))
  tmp_entity_orig = tmp_entity_orig[order(as.numeric(names(tmp_entity_orig)))]
  tmp_entity_orig = 1-tmp_entity_orig # Dissimilarity --> similarity
  tmp_entity = tmp_entity_orig / max(tmp_entity_orig)
  tmp_colors = apply(colorRamp(c("black", "yellow", "red"))(tmp_entity), 1, function(x){rgb(x[1], x[2], x[3], maxColorValue = 255)})
  names(tmp_colors) = names(tmp_entity)
  
  names(tmp_colors) = gsub("cell_", "", names(tmp_colors))
  grid = subset(grid, ID %in% names(tmp_colors))
  tmp_colors = tmp_colors[match(grid@data$ID, names(tmp_colors))]
  
  continents = raster::crop(continents, raster::extent(-180, 180, -62, 90))
  plot(continents, col = "grey80", border = NA)
  points(geoentities_meta[which(geoentities_meta$entity_ID == entity_ID), c("longitude", "latitude")], pch = 10, cex = 4, lwd = 2)
  for(cell in names(tmp_colors)){
    cell_ID = gsub("cell_", "", cell)
    plot(grid[which(grid@data$ID == cell_ID),], col = tmp_colors[cell], border = NA, add = T)
  }
  
  if(!is.null(main)){
    mtext(main, 3, 1, cex = 1.5)
  }
  if(plot_histogram){
    par(mar = c(5,4,1,1))
    hist_colors = apply(colorRamp(c("black", "yellow", "red"))(seq(0, 1, length.out = 20)), 1, function(x){
      rgb(x[1], x[2], x[3], maxColorValue = 255)
    })
    hist(tmp_entity_orig, main = NA, col = hist_colors, border = NA,  las = 1,  ylab = NA, cex.lab = 1,
         xlab = expression(bold('Source likelihood [1-Beta'[sim]*']')), breaks = seq(0, max(tmp_entity_orig), length.out = 21))
  }
}

# Visualize the aggregated source region values for a given island
plot_sourcepool_agg = function(entity_ID, plot_histogram = F, main = NULL){
  # entity_ID = focal island
  # plot_histogram = plot a histogram next to the map?
  
  load("../Data/raw_data/geoentities.RData")
  load("../Data/raw_data/geoentities_meta.RData")
  load("../Data/raw_data/continents.RData")
  load("../Data/source_matrix.RData")
  
  if(plot_histogram){layout(matrix(1:2, ncol = 2), widths = c(10,3))}
  par(mar = c(0,0,0,0))
  if(!is.null(main)){par(mar = c(0,0,3,0))}
  
  tmp_probs = source_matrix[paste(entity_ID),]
  tmp_probs_std = (tmp_probs-min(tmp_probs))/(max(tmp_probs)-min(tmp_probs))
  
  color_table = colorRamp(c("black", "yellow", "red"))(tmp_probs_std)
  color_table[is.na(color_table)] = 255
  tmp_colors = apply(color_table, 1, function(x){rgb(x[1], x[2], x[3], maxColorValue = 255)})
  tmp_colors[tmp_probs == 0] = "black"
  names(tmp_colors) = names(tmp_probs)
  tmp_geoentities = geoentities[which(geoentities$entity_ID %in% names(tmp_colors)),]
  
  continents = raster::crop(continents, raster::extent(-180, 180, -62, 90))
  plot(continents, col = "grey80", border = NA)
  points(geoentities_meta[which(geoentities_meta$entity_ID == entity_ID), c("longitude", "latitude")], pch = 10, cex = 4, lwd = 2)
  plot(tmp_geoentities, col = tmp_colors[match(tmp_geoentities$entity_ID, names(tmp_colors))], border = NA, add = T)
  
  if(plot_histogram){
    par(mar = c(5,4,1,1))
    hist_colors = apply(colorRamp(c("black", "yellow", "red"))(seq(0, 1, length.out = 20)), 1, function(x){
      rgb(x[1], x[2], x[3], maxColorValue = 255)
    })
    hist(tmp_probs, main = NA, col = hist_colors, border = NA,  las = 1,  ylab = NA, cex.lab = 1,
         xlab = expression(bold('Mean source likelihood [1-Beta'[sim]*']')), breaks = seq(0, max(tmp_probs), length.out = 21))
  }
}

# Plot the over- or under-representation of a given family
plot_representational_disharmony = function(family_names, ncol = 2){
  # family_names: character vector of family names
  
  load("../Data/raw_data/geoentities_meta.RData")
  load("../Data/raw_data/continents.RData")
  load("../Data/raw_data/checklists.RData")
  load("../Data/D_family.RData")
  
  continents = raster::crop(continents, raster::extent(-180, 180, -62, 90))  
  geoentities_plot = D_family %>% 
    filter(family %in% family_names & n_exp >= 1) %>% 
    inner_join(geoentities_meta, by = "entity_ID") %>% 
    mutate(family = factor(family, levels = family_names), present = as.factor(ifelse(n_obs > 0, "present", "absent"))) %>% 
    mutate(present = factor(present, levels = c("present", "absent"))) %>% 
    arrange(family, !is.na(D))
  
  
  ggplot() + 
    geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = "gray80")  +
    geom_point(data = geoentities_plot, mapping = aes(x = longitude, y = latitude, color = D, size = area, shape = present), stroke = 1.5) +
    scale_size_continuous(breaks = c(0.1,1,10,100,1000,10000,100000), labels = c(0.1,1,10,100,1000,10000,100000), trans = "sqrt", guide = "legend", 
                          limits = c(1, 10000000), range = c(2,30), name = expression(bold(paste('Area [', km^2, ']')))) +
    scale_color_gradient2(low = "#ff553f", mid = "gray30", high = "#1c7bff", limits=c(0, 1), midpoint = 0.5, na.value = "black", 
                          breaks = c(1,0.5,0), labels = c("1.0 (over-represented)", "0.5 (harmonic)", "0.0 (under-represented)"))  +
    scale_shape_manual(values = c(1, 4), name = "Occurrence", drop = F) +
    labs(x = "Longitude", y = "Latitude") +
    facet_wrap(~ family, ncol = ncol) + 
    guides(colour = guide_colorbar(order = 1, title = "Representational \ndisharmony"),
           shape = guide_legend(order = 2),
           size = guide_legend(order = 3)) +
    theme_bw() + 
    theme(strip.text = element_text(size=15),
          axis.text = element_text(size=11),
          axis.title = element_text(size=15),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 15, face = "bold"))
}

# Plot the the overall bias of island floras compared to their source regions
plot_compositional_disharmony = function(){
  load("../Data/raw_data/geoentities_meta.RData")
  load("../Data/raw_data/continents.RData")
  load("../Data/raw_data/checklists.RData")
  load("../Data/D_family.RData")
  
  continents = raster::crop(continents, raster::extent(-180, 180, -62, 90))  
  geoentities_plot = D_family %>% 
    inner_join(geoentities_meta, by = "entity_ID") %>% 
    mutate(deviance = abs(D-0.5)) %>% 
    group_by(entity_ID) %>% 
    summarize(deviance_avg = median(deviance, na.rm = T), longitude = mean(longitude), latitude = mean(latitude), area = mean(area), geology = mean(geology)) %>% 
    mutate(geology = recode_factor(geology, `6` = "volcanic", `1` = "atoll", `3` = "tectonic uplift                 "))
  
  ggplot() +
    geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = "gray80")  +
    geom_point(data = geoentities_plot, mapping = aes(x = longitude, y = latitude, color = deviance_avg, size = area, shape = geology), stroke = 1.5) +
    scale_size_continuous(breaks = c(0.1,1,10,100,1000,10000,100000), labels = c(0.1,1,10,100,1000,10000,100000), trans = "sqrt", guide = "legend", 
                          limits = c(1, 10000000), range =  c(2,30), name = expression(bold(paste('Area [', km^2, ']')))) +
    scale_color_gradient2(low = "#c1cbff", mid = "#ffe5a8", high = "#ff0000", midpoint = 0.25, limits=c(0, 0.5)) +
    scale_shape_manual(values = c(17, 1, 15), name = "Island type", drop = F) +
    labs(x = "Longitude", y = "Latitude") +
    guides(colour = guide_colorbar(order = 1, title = "Compositional \ndisharmony"),
           shape = guide_legend(order = 2),
           size = guide_legend(order = 3)) +
    theme_bw() + 
    theme(strip.text = element_text(size=15),
          axis.text = element_text(size=11),
          axis.title = element_text(size=15),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 15, face = "bold"))
}

# Plot dot-whisker plots of standardized model coefficients
plot_model_coefs = function(results_table, var_name, show_x_axis = F, show_title = F, limits = c(-1,1), breaks = seq(-1, 1, 0.25)){
  effect_size_tmp = filter(results_table, variable == var_name)
  plot_tmp = ggplot(effect_size_tmp, aes(x = term, y = estimate)) +
    geom_hline(yintercept = 0, color = "gray50", linetype = 2, size = 0.5) +
    geom_pointrange(mapping = aes(ymin = lower, ymax = upper), size = 1.3, fatten=1.3, colour = "red") +
    labs(x="" , y="Standardized coefficient") +
    coord_flip(clip = 'off') +
    scale_y_continuous(limits = limits, breaks = breaks) +
    theme_light() 
  
  if(show_title == T){
    plot_tmp = plot_tmp + ggtitle(var_name)
  }
  
  if(show_x_axis == F) {
    plot_tmp = plot_tmp +
      theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
            axis.title.x=element_blank(),
            axis.text.x=element_blank())
  } else {
    plot_tmp = plot_tmp +
      theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))
  }
  return(plot_tmp)
}

resize_heights <- function(g, heights = rep(1, length(idpanels))){
  idpanels <- unique(g$layout[grepl("panel",g$layout$name), "t"])
  g$heights <- grid:::unit.list(g$heights)
  g$heights[idpanels] <- unit.c(do.call(unit, list(heights, 'null')))
  g
}

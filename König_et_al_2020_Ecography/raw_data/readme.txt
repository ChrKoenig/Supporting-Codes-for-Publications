Data description for König et al. 2020 - Source pools and disharmony of the world's island floras. Ecography (https://doi.org/10.1111/ecog.05174)

============================================================================================

checklists.RData - Species checklists
	entity_ID - unique geographical region ID
	work_ID - unique species ID
	family - plant family 
	
continents.RData - Shape file of continents

env_grid.RData - Extracted environmental information for global equal-area grid
	entity_ID - unique geographical region ID
	longitude - longitude in degrees
	latitude - latitude in degrees
	T_mean - Mean annual temperature 
	T_var - Temperature seasonality
	P_mean - Mean annual precipitation
	P_var - Precipitation seasonality
	
family_traits_final.RData - Family traits used in the analysis of representational disharmony
	supergroup - Taxonomic supergroup (Angiospermae/Gymnospermae/Pteridophyta)
	family - plant family
	trait_name - name of the trait_name
	trait_value - value of the trait
	source - source of information (botanical literature/aggregation from species-level data)
	
geoentities.RData - Shape file of geographical regions

geoentities_meta.RData - Additional information on geographical regions
	entity_ID - unique geographical region ID
	geo_entity - name of the geographical region
	entity_class - category of the geographical region (Mainland/Island/both)
	longitude - longitude in degrees
	latitude - latitude in degrees
	area - area in km²
	dist - distance to the nearest mainland in km
	elev - maximum elevation in m
	T_mean - Mean annual temperature 
	T_var - Temperature seasonality
	P_mean - Mean annual precipitation
	P_var - Precipitation seasonality
	arch_lvl_1 - archipelago ID (coarse resolution)
	arch_lvl_2 - archipelago ID (medium resolution)
	arch_lvl_3 - archipelago ID (fine resolution)
	geology - geological origin (atoll, volcanic, tectonic uplift)

grid.RData - Shape file of global equal-area grid
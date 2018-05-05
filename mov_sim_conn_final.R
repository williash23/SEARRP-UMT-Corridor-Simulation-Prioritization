###############################################################################
#  Connectivity through movement simulation for bookend scenarios.
#  January 19, 2018
#  Script to generate simulated movement paths for "bookend" spp, incorporating landscape
#   characteristics, and evaluating probability of cell use as a measure of connectivity.
#  Sara Williams
###############################################################################



# =============================================================================
#  Notes on model implementation, assumptions, and areas for improvement.
# =============================================================================

	# ----------------------
	#  Initial locations (where an animal starts) is constrained to a protected area 
	
	# ----------------------
	#  Barriers (e.g., main paved roads) incorporated with relatively high resistance. Integrated 
	#   within forest cover and resistance layer.


	# ----------------------
	#  Scenarios:
	#   All combinations of maximum and minimum of 4 parameters:
	#   1. Directional persistence: low, high 
	#   2. Step length: short, long
	#   3. Permeability: resistance low, resistance high
	#   4. Perceptual range: short, long



# =============================================================================
#  Load packages.
# =============================================================================
library(sf)
library(sp)
library(raster)
library(rgdal)
library(rgeos)
library(dplyr)
library(SiMRiv)
library(rasterVis)


# =============================================================================
#  Specify parameters for simulations and CRW
# =============================================================================

	# ----------------------
	#  Number of simulations, individuals, and steps per individual
	nsim <- 1000
	ni <- 1
	nstep <- 3000

	# ----------------------
	#  Concentration of turning angles between steps (i.e., directional persistence), which
	#   ranges from 0 to 1; 1 = straight line (most persistent).
	#   This will be a vector of length = nscen
	#turn_conc <- c(0.2, 0.9)
	turn_conc <- c(0.8)
	
	# ----------------------
	#  Perceptual range of spp in meters.
	#   This will be a vector of length = nscen
	# percep_range <- c(910, 2730) # meters
	percep_range <- c(500) # meters
	
	# ----------------------
	#  Step length in meters, which should be less than perceptual range.
	step_l <- c(250)

	# ----------------------
	#  Permeability: resistance value for each cell based on forest cover quality (ACD value, or
	#   reclassified ACD value)
	# 0 = non-forest
	# 1 = forest
	# 2 = roads
	# 3 = ocean
	# 4 = PA (assumption highest quality forest)
	resist_l <- c(0.4, 0.1, 0.25, 0.50, 0)
	resist_h <- c(0.98, 0.05, 0.75, 0.50, 0.0)
	#resist_val <- as.data.frame(cbind(resist_l, resist_h))
	resist_val <- as.data.frame(resist_h)
	
	
	
# =============================================================================
#  Load base layer input data.
# =============================================================================

	# ----------------------
	#  Boundaries
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/trans_crop_proj/border_sabah_d.Rdata")
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/trans_crop_proj/border_sarawak_d.Rdata")
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/processed_spat_data/trans_crop_proj/border_kali_d.Rdata")
	load("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/border_sabah_sf.Rdata")
	border_sabah_sp  <- as(border_sabah_sf, "Spatial")
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/main_sabah_sf.Rdata")
	main_sabah_sp <- as(main_sabah_sf, "Spatial")
	
	# ----------------------
	# Forest cover
	#for_cov <- raster("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/for_cov_conn.grd")
	# ----------------------
	# Forest cover - larger buffer size
	for_cov <- raster("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/for_cov_conn_lg.grd")
	# 0-1
	
	# ----------------------
	#  Exiting protected areas, clustered (within 5m of each other become single PA)
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/ssk_pa_near_sf_clust.Rdata")
	
	# ----------------------
	#  Select only PA cluster with area of 5000 ha.
	pa_sf <- ssk_pa_near_sf_clust %>%
		mutate(size = as.integer(st_area(.)*0.0001)) %>%
		dplyr::filter(size > 5000) %>%
		dplyr::select(size)
	pa_sp  <- as(pa_sf, "Spatial")
	
	

# =============================================================================
#  Create empty object for storage of loop outputs.
# =============================================================================
	
	# ----------------------
	#  Create raster following template of study area (cell values are empty)
	r_mat <- matrix(0, nrow(for_cov), ncol(for_cov))
	r_template <- raster(r_mat)
	extent(r_template) <- extent(for_cov)
	projection(r_template) <- CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0") 


	
# =============================================================================
#  Function to normalize raster to 0 to 1 scale.
# =============================================================================
	
	# ----------------------	
	#  Normalize to 0 to 1 scale.
	range01 <- function(x){(x-x_min)/(x_max - x_min)}
	#  Must always also set up min and max for each raster before using function, e.g.,:
	# x_min <- move_smooth@data@min
	# x_max <-move_smooth@data@max
	
	

# =============================================================================
#  Function for calculation of relative probability of initial location.
# =============================================================================
	
	# ----------------------
	# Function for use in selecting initial location of individual.
	probsel <- function(probrast, ni){
		x <- getValues(probrast)
		x[is.na(x)] <- 0
		vec_cells <- seq(1:length(probrast))
		samp <- sample(vec_cells,  ni, replace = F, prob=x)
		samprast <- raster(probrast)
		samprast[samp] <- 1 
		samp_pts <- rasterToPoints(samprast, fun = function(x){x > 0})
		samp_pts <- SpatialPoints(samp_pts)
		crs(samp_pts) <- CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
		return(samp_pts)
		}

	# ----------------------
	#  Rasterize existing protected areas according to this raster template.
	pa_r <- rasterize(pa_sp, r_template, field = pa_sp@data$size)
	
	# ----------------------
	#  Assign relative probability of initial occurrence to each grid cell from relative predicted 
	#   density using function assigned above (probsel()) in matrix format.
	#  Initial occurrence only from within protected areas, probability based on size.
	probrast <- pa_r/sum(getValues(pa_r), na.rm = T)



# =============================================================================
#  Run simulation for all selected spp using parameters specified above.
# =============================================================================

	# ----------------------
	#  Portion to generate resistance matrix.
	#   Nested loops run over every possible combination of parameters.
	
	for(i in 1:length(turn_conc)){
		for(j in 1:length(percep_range)){
			for(k in 1:length(step_l)){
				for(q in 1:ncol(resist_val)){
		
				# ----------------------
				#  Reclassify forest cover into a resistance map using values specified for 
				#   each spp permeability.
				rc_val_for <- c(0, resist_val[1,q],
					1, resist_val[2,q],
					2, resist_val[3,q],
					3, resist_val[4,q], 
					4, resist_val[5,q])
				rc_mat_for <- matrix(rc_val_for, 
					ncol = 2, 
					byrow = TRUE)
				resist <- raster::reclassify(for_cov, rc_mat_for)
				
				
					# ----------------------
					#  Portion to simulate movement paths for each individual.
					#   Loop w runs over the number of desired simulated movement paths (i.e., an individual).
					
						# ----------------------
						#  Run animal movement simulations for each spp and store output in a list.
						sim_out_list <- list()
						for(w in 1:nsim){
							
							# ----------------------
							#  Randomly select starting location for individuals, which is constrained to protected areas 
							#   and the spp range. 
							init_loc <- probsel(probrast, ni)
							init_loc_df <- as.data.frame(init_loc) 
							init_loc_m <- as.matrix(init_loc_df[,1:2])

							# ----------------------
							#  Use functions from SiMRiv packge to generate simulated paths of CRW for each individual.
							#  Generate the spp' CRW using specified parameters.
							corr_rw <- species(state.CRW(turn_conc[i]) * percep_range[j] + step_l[k])
							
							# ----------------------
							#  Simulate each individual movement path from that spp' CRW, specified number of steps, its 
							#   random starting location, and permemability.
							sim <- simulate(corr_rw, nstep, resist = resist, coords = c(init_loc_m[1], init_loc_m[2]))

							# ----------------------			
							#   Convert each simulated path to data frame then to an sf object.
							sim_move_df <- as.data.frame(sim[,1:2])
							colnames(sim_move_df) =  c("x", "y")
							coordinates(sim_move_df) = ~ x + y
							crs(sim_move_df) <- CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
							sim_move_sf_tmp <- st_as_sf(sim_move_df)

							# ----------------------			
							#   Format column names of simulated movement path sf object.			
							sim_move_sf <- sim_move_sf_tmp %>%
								mutate(sim_path = w) %>%
								mutate(step_num = 1:n()) 
								
							# ----------------------			
							#  Assign simulated movement path sf object to specific name for each spp.
							nam1 <- paste("sim_move_sf", "crw", w, sep = "_") 
							assign(nam1, sim_move_sf)
							
							# ----------------------			
							#  Fill in simulation list with each spp' simulated movement path sf object.
							sim_out_list[[w]] <- as.name(nam1)
							}

					# ----------------------			
					#  Bind the list of simulated paths together to make a collection of total = nsim paths of 
					#   animal movement for each spp.
					sim_move_bind <- do.call(rbind, sim_out_list)
					sim_move_bind_sp <- as(sim_move_bind, "Spatial")
					
					# ----------------------	
					#  Tabulate the number of steps (i.e., point locations from simulated movement path) that occur
					#   in each grid cell as a measure of each cell's use by the spp.
					sim_move_r <- r_template
					tab <- table(cellFromXY(sim_move_r, sim_move_bind_sp))
					sim_move_r[as.numeric(names(tab))] <- tab
					
					# ----------------------	
					#  Save outputs for each spp.
					nam2 <- paste("sim_move_bind", i, j, k, q, sep = "_")
					assign(nam2, sim_move_bind)
					nam3 <- paste(nam2, "Rdata", sep = ".")
					save(nam2, file = paste("C:/Users/saraw/Documents/SEARRP_Analyses/move_sim/moves_in", nam3, sep = "/"))
					nam3 <- paste("sim_move_r", i, j, k, q, sep = "_")
					writeRaster(sim_move_r, paste("C:/Users/saraw/Documents/SEARRP_Analyses/move_sim/moves_in", nam3, sep = "/"))
					
					}	
				}
			}
		}
	
	
	
# =============================================================================
#  Combine outputs from each spp (all simulations per spp).
# =============================================================================

	# ----------------------	
	#  Stack simulated movement rasters for each spp. 
	# setwd("C:/Users/saraw/Documents/SEARRP_Analyses/move_sim/directed_mov")
	# img <- list.files(pattern='\\.grd$')
	# all_sim_move <- stack(img)
	
	# ----------------------	
	#  Sum values over cells across all RasterLayers in stack to give a cumulative value of how cells 
	#   are used across simulated movement.
	# all_sim_move_sum <- raster::calc(all_sim_move, sum, na.rm = T)
	
	# ----------------------	
	#  Or, if using only output from the loop, assign it to name for further processing.
	all_sim_move_sum <- sim_move_r
	
	# ----------------------	
	#  Remove already protected areas and crop to Sabah.
	all_sim_move_sum_m_pa <- raster::mask(all_sim_move_sum, pa_sp, inverse = TRUE, updatevalue = 0)
	all_sim_move_sum_c <- raster::crop(all_sim_move_sum_m_pa, main_sabah_sp)
	all_sim_move_sum_m <- raster::mask(all_sim_move_sum_c, main_sabah_sp)
	
	# ----------------------	
	#  Normalize to 0 to 1 scale,  if wanted for use as a conservation feature input.
	x_min <- all_sim_move_sum_m@data@min
	x_max <- all_sim_move_sum_m@data@max
	moves_in <- raster::calc(all_sim_move_sum_m, range01)
	
	# ----------------------	
	#  Save outputs of aggregating all spp and all simulations.
	# writeRaster(moves_in, "C:/Users/saraw/Documents/SEARRP_Analyses/move_sim/moves_in.grd")
	
	# ----------------------	
	#  Smooth using a focal window
	moves_sm <- focal(all_sim_move_sum_m,  w = matrix(1, 3, 3))
	
	# ----------------------	
	#  Normalize to 0 to 1 scale,  if wanted for use as a conservation feature input.
	x_min <- moves_sm@data@min
	x_max <-moves_sm@data@max
	moves_in_sm <- raster::calc(moves_sm, range01)
	
	# ----------------------	
	#  Save outputs of aggregating all spp and all simulations.
	# writeRaster(moves_in_sm, "C:/Users/saraw/Documents/SEARRP_Analyses/move_sim/moves_in_sm.grd")
	

	
# =============================================================================
#  Plot using RasterVis.
# =============================================================================

	# # ----------------------	
	# #   Set color palette theme
	r_theme <- rasterTheme(brewer.pal(10, "Greens"))
	r_theme$regions$col <- c("grey","grey","grey","grey","grey", 
		"lightgreen", "lightgreen", "lightgreen", "lightgreen", "lightgreen",
		"black", "black", "black", "black", "black", 
		"lightblue", "lightblue", "lightblue", "lightblue", "lightblue", 
		"darkgreen", "darkgreen", "darkgreen", "darkgreen", "darkgreen")
	r_theme <- rasterTheme(brewer.pal(5, "Reds"))
	r_theme$panel.background$col <- "white"
	
	
	# ----------------------	
	plot_move_r <- levelplot(moves_in, par.settings = r_theme, 
		#main = "Probability of use during species movement \n",
		xlab= "Longitude (UTM)",
		ylab="Latitude (UTM)",
		margin = FALSE) +
		layer(sp.polygons(border_sabah_d, lwd = 0.8, col = 'grey40')) + 
		layer(sp.polygons(border_sarawak_d, lwd = 0.8, col = 'grey40', fill = 'gray40')) +
		layer(sp.polygons(border_kali_d, lwd = 0.8, col = 'grey40', fill = 'gray40')) #+
	#	layer(sp.polygons(pa_sp, lwd = 0.8, col = 'grey40', fill = 'darkorange2')) 
	plot_move_r
	
	

# =============================================================================	
###############################################################################
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
#library(rasterVis)
library(doParallel)


# =============================================================================
#  Specify parameters for simulations and CRW
# =============================================================================

	# ----------------------
	#  Number of simulations, individuals, and steps per individual
	nsim <- 1000
	nrep <- 100
	ni <- 1
	nstep <- 250

	# ----------------------
	#  Concentration of turning angles between steps (i.e., directional persistence), which
	#   ranges from 0 to 1; 1 = straight line (most persistent).
	#   This will be a vector of length = nscen
	turn_conc <- c(0.92)
	
	# ----------------------
	#  Perceptual range of spp in meters.
	#   This will be a vector of length = nscen
	percep_range <- c(1500) # meters
	
	# ----------------------
	#  Step length in meters, which should be less than perceptual range.
	step_l <- 500

	
	
# =============================================================================
#  Load base layer input data.
# =============================================================================

	# ----------------------
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/main_sabah_sf.Rdata")
	main_sabah_sp <- as(main_sabah_sf, "Spatial")
	
	# ----------------------
	# Forest cover - larger buffer size
	for_cov <- raster("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/ext_for_cov_rds.grd")
	
	# ----------------------
	#  Exiting protected areas, clustered (within 5m of each other become single PA)
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/pa_sp_move_sim.Rdata")
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/pa_sf_move_sim.Rdata")
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/pa_sp_move_sim_smooth.Rdata")
	
	
	
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
	pa_r_tmp <- rasterize(pa_sp, r_template, field = pa_sp@data$size)
	pa_r <- mask(pa_r_tmp, pa_sp_smooth, inverse = TRUE)
	
	# ----------------------
	#  Assign relative probability of initial occurrence to each grid cell from relative predicted 
	#   density using function assigned above (probsel()) in matrix format.
	#  Initial occurrence only from within protected areas, probability based on size.
	probrast <- pa_r/sum(getValues(pa_r), na.rm = T)
	probrast <- raster::crop(probrast, extent(231556, 751756, 400000, 783428))
	
	
	
# =============================================================================
#  Set up resistance matrix.
# =============================================================================

	# ----------------------
	#  Permeability: 
	resist_l <- c(0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.001, 0.5)
	resist_val <- as.data.frame(resist_l)

	# ----------------------
	#  Reclassify forest cover into a resistance map using values specified for 
	#   each spp permeability.
	rc_val_for <- c(1, resist_val[1,1],
		2, resist_val[2,1],
		3, resist_val[3,1],
		4, resist_val[4,1],
		5, resist_val[5,1],
		6, resist_val[6,1],
		7, resist_val[7,1],
		8, resist_val[8,1], 
		9, resist_val[9,1])
		rc_mat_for <- matrix(rc_val_for, 
		ncol = 2, 
		byrow = TRUE)
	resist <- raster::reclassify(for_cov, rc_mat_for)

	
# =============================================================================
#  Run simulation.
# =============================================================================

	# ----------------------
	#  Initiate cluster for parallel processing.
	cl <- makeCluster(detectCores()-1)
	registerDoParallel(cl)

	# ----------------------
	#  Initiatialize stack to hold output
	out_stack <- stack()
	
	for(q in 1:nrep){
		
		move_sim_paths <- foreach(w = 1:nsim,  .packages=c("dplyr", "sf", "sp", "raster", "SiMRiv")) %dopar% {
			
			# ----------------------
			#  Randomly select starting location for individuals, which is constrained to protected areas 
			#   and the spp range. 
			init_loc <- probsel(probrast, ni)
			init_loc_df <- as.data.frame(init_loc) 
			init_loc_m <- as.matrix(init_loc_df[,1:2])

			# ----------------------
			#  Generate the CRW using specified parameters.
			spp_crw <- SiMRiv::species(state(turn_conc, perceptualRange("cir", percep_range), step_l, "CorrelatedRW"))

			# ----------------------
			#  Simulate a single movement path from CRW model
			sim <- SiMRiv::simulate(spp_crw, nstep, resist = resist, coords = c(init_loc_m[1], init_loc_m[2]))

			# ----------------------			
			#   Convert each simulated path to data frame then to an sf object.
			sim_move_df <- as.data.frame(sim[,1:2])
			colnames(sim_move_df) =  c("x", "y")
			coordinates(sim_move_df) = ~ x + y
			crs(sim_move_df) <- CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
			
			# ----------------------			
			#   Format column names of simulated movement path sf object.			
			sim_move_sf_tmp <- st_as_sf(sim_move_df) %>%
				mutate(sim_path = w) %>%
				mutate(step_num = 1:n()) 

			first_loc <- sim_move_sf_tmp %>%
				slice(1)
			first_pa <- as.numeric(raster::extract(pa_r, as(first_loc, "Spatial")))
			all_pa <- raster::extract(pa_r, as(sim_move_sf_tmp, "Spatial"))
			
			sim_move_sf_tmp2 <- cbind(sim_move_sf_tmp, all_pa) %>%
				mutate(sim_path = w) %>%
				mutate(step_num = 1:n()) %>%
				mutate(start_pa = first_pa) 
			sim_move_sf_tmp2[is.na(sim_move_sf_tmp2)] <- first_pa
			
			sim_move_df_tmp <- sim_move_sf_tmp2 %>%
				rowwise() %>%
				mutate(cont_path = ifelse(all_pa == start_pa, 1, 0)) %>%
				as.data.frame() 
			
			end_step_tmp <- which(sim_move_df_tmp$cont_path == 0)[1]
			
			if(is.na(end_step_tmp)){
				end_step <- nstep}else{
					end_step <- end_step_tmp
					}
			
				# ----------------------
				#  Select only the last step for dispersers ending location
				sim_move_df_fin <- sim_move_df_tmp[1:end_step,]
				last_loc_n <- nrow(sim_move_df_fin)
				last_loc <- sim_move_df_fin[last_loc_n,]
				stayed <- ifelse(last_loc$start_pa - last_loc$all_pa != 0, 0, 5)
				
				if(stayed < 5){
					sim_move_sf <- NA}else{
						sim_move_sf <- st_as_sf(sim_move_df_fin, sf_column_name = "geometry") 
						}
				
			return(sim_move_sf)
			}	

		# ----------------------			
		#  Remove "null" paths that only stayed within the same PA	
		good_paths <- na.omit(move_sim_paths)
		length_idx <- length(good_paths)
		
		if(length_idx  == 0){
			sim_move_r <- r_template
			sim_move_r[] <- 0 } else {
				
				# ----------------------			
				#  Bind the list of simulated paths together to make a collection of total = nsim paths of 
				#   animal movement for each spp.
				paths_bind <- do.call(rbind, good_paths)
				paths_bind_sp <- as(paths_bind, "Spatial")
				
				# ----------------------	
				#  Tabulate the number of steps (i.e., point locations from simulated movement path) that occur
				#   in each grid cell as a measure of each cell's use by the spp.
				sim_move_r <- r_template
				tab <- table(cellFromXY(sim_move_r, paths_bind_sp))
				sim_move_r[as.numeric(names(tab))] <- tab
				}
	
		out_stack <- stack(out_stack, sim_move_r)
			
	}

	
	
	out_sum <- sum(out_stack)

	out_moves_tmp <- stack(out1_sum, corr_out) 
	out_moves <- sum(out_moves_tmp, na.rm = TRUE)
	

	
	
	moves_in_c <- raster::crop(out_moves, main_sabah_sp)
	moves_in_m <- raster::mask(moves_in_c, main_sabah_sp)
	moves_in_m2 <- raster::mask(moves_in_m, pa_sp, inverse = TRUE)
	
	plot(moves_in_m2)
		
	m <- matrix(1, ncol=3, nrow=3)
	moves_sm <- focal(moves_in_m2, m, fun="mean", na.rm=TRUE, NAonly=TRUE, pad=TRUE) 

	

	#moves_in_m2 <- raster::mask(moves_in_m, mal_dan_bound_sp, inverse = TRUE)
	
	
	
	out_sum[is.na(out_sum)] <- 0
	
	x_min <- corr_r@data@min
	x_max <- corr_r@data@max
	corr_fac <- raster::calc(corr_r, range01)
	corr_fac[corr_fac == 0] <- NA
	plot(corr_fac)
	
	 
	
	##################################
	# ----------------------	
	#  Use only southwest corner of Sabah
	load(file = "C:/Users/saraw/Desktop/planning_unit_grids/sw_sabah_grid.Rdata")
	sw <- st_buffer(st_union(sw_sabah_grid), 5000)
	area_thresh <- units::set_units(5000, km^2)
	sw_fill <- fill_holes(sw, threshold = area_thresh)
	sw_sp <- as(sw_fill, "Spatial")
	
	# ----------------------	
	#  Normalize to 0 to 1 scale,  if wanted for use as a conservation feature input.
	sim_move_r_c <- raster::crop(sim_move_r, sw_sp, updatevalue = 0)
	sim_move_r_m <- raster::mask(sim_move_r_c, sw_sp, updatevalue = 0)
	
	#moves_in_sm_c <- raster::crop(sim_move_r_c, main_sabah_sp)
	#moves_in_sm_m <- raster::mask(moves_in_sm_c, main_sabah_sp)
	
	x_min <- sim_move_r_m@data@min
	x_max <- sim_move_r_m@data@max
	moves_in <- raster::calc(sim_move_r_m, range01)
	moves_in[moves_in == 0] <- NA
	

	# ----------------------	
	#  Save outputs of aggregating all spp and all simulations.
	# writeRaster(moves_in, "C:/Users/saraw/Desktop/moves_in.grd")
	
	# ----------------------	
	#  Smooth using a focal window
	moves_sm <- focal(all_sim_move_sum_m,  w = matrix(1, 3, 3))
	

	
# =============================================================================
#  Plot using RasterVis.
# =============================================================================

	# # ----------------------	
	# #   Set color palette theme
	r_theme <- rasterTheme(brewer.pal(10, "Reds"))
	r_theme$regions$col <- c("grey","grey","grey","grey","grey", 
		"lightgreen", "lightgreen", "lightgreen", "lightgreen", "lightgreen",
		"black", "black", "black", "black", "black", 
		"lightblue", "lightblue", "lightblue", "lightblue", "lightblue", 
		"darkgreen", "darkgreen", "darkgreen", "darkgreen", "darkgreen")
	r_theme <- rasterTheme(brewer.pal(5, "Reds"))
	r_theme$panel.background$col <- "white"
	
	
	# ----------------------	
	plot_move_r <- levelplot(probrast, par.settings = r_theme, 
		#main = "Probability of use during species movement \n",
		#xlab= "Longitude (UTM)",
		#ylab="Latitude (UTM)",
		margin = FALSE) #+
		#layer(sp.polygons(border_sabah_d, lwd = 0.8, col = 'grey40')) #+ 
		#layer(sp.polygons(border_sarawak_d, lwd = 0.8, col = 'grey40', fill = 'gray40')) +
		#layer(sp.polygons(border_kali_d, lwd = 0.8, col = 'grey40', fill = 'gray40')) #+
	#	layer(sp.polygons(pa_sp, lwd = 0.8, col = 'grey40', fill = 'darkorange2')) 
	plot_move_r
	
	

# =============================================================================	
###############################################################################
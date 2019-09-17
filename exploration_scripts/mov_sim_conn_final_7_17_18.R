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
library(doParallel)


# =============================================================================
#  Specify parameters for simulations and CRW
# =============================================================================

	# ----------------------
	#  Number of simulations, individuals, and steps per individual
	nsim <- 5000
	ni <- 1
	nstep <- 500

	# ----------------------
	#  Concentration of turning angles between steps (i.e., directional persistence), which
	#   ranges from 0 to 1; 1 = straight line (most persistent).
	#   This will be a vector of length = nscen
	#turn_conc <- c(0.2, 0.9)
	turn_conc <- c(0.92)
	
	# ----------------------
	#  Perceptual range of spp in meters.
	#   This will be a vector of length = nscen
	# percep_range <- c(910, 2730) # meters
	percep_range <- c(1000) # meters
	
	# ----------------------
	#  Step length in meters, which should be less than perceptual range.
	step_l <- c(500)

	
	
# =============================================================================
#  Load base layer input data.
# =============================================================================

	# ----------------------
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/main_sabah_sf.Rdata")
	main_sabah_sp <- as(main_sabah_sf, "Spatial")
	
	# ----------------------
	# Forest cover - larger buffer size
	for_cov_tmp <- raster("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/for_cov_conn_lg.grd")
	for_cov_tmp[for_cov_tmp == 3] <- NA
	rc_val <- c(0, 1,
		1, 2,
		2, 3,
		4, 4)
	rc_mat <- matrix(rc_val, 
		ncol = 2, 
		byrow = TRUE)
	for_cov <- raster::reclassify(for_cov_tmp, rc_mat)
	# 1 = non-forest
	# 2 = forest
	# 3 = roads
	# 4 = PA (assumption highest quality forest)
	
	# ----------------------
	#  Permeability: 
	resist_l <- c(0.2, 0.02, 0.2, 0)
	resist_val <- as.data.frame(resist_l)
	
	# ----------------------
	#  Exiting protected areas, clustered (within 5m of each other become single PA)
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/ssk_pa_near_sf_clust.Rdata")
	
	# ----------------------
	#  Select only PA cluster with area of 5000 ha.
	pa_tmp1 <- st_intersection(ssk_pa_near_sf_clust, 
		st_set_crs(st_as_sf(as(raster::extent(231556, 751756, 400000, 783428), 
		"SpatialPolygons")), st_crs(ssk_pa_near_sf_clust)))
	pa_sf <- pa_tmp1 %>%
		mutate(size = as.integer(st_area(.)*0.0001)) %>%
		#dplyr::filter(size > 5000) %>%
		dplyr::select(size)
	pa_sp  <- as(pa_sf, "Spatial")
	
	
	
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/mal_dan_bound.Rdata")
	
	tmp1 <- st_buffer(mal_dan_bound, 100)
	area_thresh <- units::set_units(400, km^2)
	tmp2 <- fill_holes(tmp1, area_thresh)
	mal_dan_bound_sp <- as(tmp2, "Spatial")

	
	

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
	probrast <- raster::crop(probrast, extent(231556, 751756, 400000, 783428))
	
	
# =============================================================================
#  Run simulation for all selected spp using parameters specified above.
# =============================================================================

	# ----------------------
	#  Reclassify forest cover into a resistance map using values specified for 
	#   each spp permeability.
	rc_val_for <- c(1, resist_val[1,1],
	2, resist_val[2,1],
	3, resist_val[3,1],
	4, resist_val[4,1])
	rc_mat_for <- matrix(rc_val_for, 
	ncol = 2, 
	byrow = TRUE)
	resist <- raster::reclassify(for_cov, rc_mat_for)


	# ----------------------
	#  Initiate cluster for parallel processing.
	cl <- makeCluster(detectCores()-1)
	registerDoParallel(cl)

	moves_rep <- foreach(q = 1:nrep,  .packages=c("dplyr", "sf", "sp", "raster", "SiMRiv", "doParallel")) %dopar% {

		move_sim_conn_locs <- foreach(w = 1:nsim,  .packages=c("dplyr", "sf", "sp", "raster", "SiMRiv")) %dopar% {
		
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

			first_loc <-  sim_move_sf_tmp %>%
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
				end_step <- nstep} else{
					end_step <- end_step_tmp
					}

			# ----------------------
			#  Select only the last step for dispersers ending location
			sim_move_df_fin <- sim_move_df_tmp[1:end_step,]
			sim_move_sf <- st_as_sf(sim_move_df_fin, sf_column_name = "geometry") 

			return(sim_move_sf)
			}	


		# ----------------------			
		#  Bind the list of simulated paths together to make a collection of total = nsim paths of 
		#   animal movement for each spp.
		sim_move_bind <- do.call(rbind, move_sim_conn_locs)
		sim_move_bind_sp <- as(sim_move_bind, "Spatial")
		
		# ----------------------	
		#  Tabulate the number of steps (i.e., point locations from simulated movement path) that occur
		#   in each grid cell as a measure of each cell's use by the spp.
		r_template_small <- disaggregate(r_template, fact = 3)
		sim_move_r <- r_template
		tab <- table(cellFromXY(sim_move_r, sim_move_bind_sp))
		sim_move_r[as.numeric(names(tab))] <- tab
		
		# ----------------------	
		#  Smooth using a focal window
		moves_sm <- focal(sim_move_r,  w = matrix(1, 9, 9))
		
		moves_in_c <- raster::crop(moves_sm, main_sabah_sp)
		moves_in_m <- raster::mask(moves_in_c, main_sabah_sp)
		moves_in_m2 <- raster::mask(moves_in_m, mal_dan_bound_sp, inverse = TRUE)
		moves_in_m3 <- raster::mask(moves_in_m2, pa_sp, inverse = TRUE)
		
	return(moves_in_m3)
	}

	
	tst <- stack(moves_rep)
	tst[is.na(tst)] <- 0
	tst2 <- mean(tst, na.rm = TRUE)
	
	
	x_min <- moves_in_m3@data@min
	x_max <- moves_in_m3@data@max
	moves_in <- raster::calc(moves_in_m3, range01)
	moves_in[moves_in == 0] <- NA
	
	
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
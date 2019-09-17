###############################################################################
#  Connectivity through movement simulation for bookend scenarios.
#  January 19, 2019; last updated June 23, 2019
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
#library(doParallel)
library(fasterize)
library(gurobi)
library(prioritizr)

setwd("C:/Users/saraw/Documents/Prioritization/")



# =============================================================================
#  Specify parameters for simulations and CRW
# =============================================================================

	# ----------------------
	#  Number of simulations, individuals, and steps per individual
	nsim <- 5000
	nrep <- 10
	ni <- 1
	nstep <- 500

	# ----------------------
	#  Concentration of turning angles between steps (i.e., directional persistence), which
	#   ranges from 0 to 1; 1 = straight line (most persistent).
	#   This will be a vector of length = nscen
	turn_conc <- c(0.96)
	
	# ----------------------
	#  Perceptual range of spp in meters.
	#   This will be a vector of length = nscen
	percep_range <- c(2000) # meters
	
	# ----------------------
	#  Step length in meters, which should be less than perceptual range.
	step_l <- 500

	# ----------------------
	#  Corridor/location weighting values
	dist_wt <- 1
	smaller_pa_wt <- 2
	larger_pa_wt <- 5



# =============================================================================
#  Load base layer input data.
# =============================================================================

	# ----------------------
	#  Mainland Sabah boundary and planning units
	load(file = "study_area_boundaries_and_pas/main_sabah_sf.Rdata")
	main_sabah_sp <- as(main_sabah_sf, "Spatial")
	pu_r <- raster("mainland_sabah_planning_units.tif")
	
	# ----------------------
	# Forest cover - larger buffer size
	for_cov <- raster("ext_for_cov_rds.tif")
	
	# ----------------------
	#  Existing protected areas, clustered (within 5m of each other become single PA)
	load(file = "pa_corr_sim.Rdata")

	prj <- "+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
	
# =============================================================================
#  Create empty object for storage of loop outputs.
# =============================================================================
	
	# ----------------------
	#  Create raster following template of study area (cell values are empty)
	r_mat <- matrix(0, nrow(for_cov), ncol(for_cov))
	r_template <- raster(r_mat)
	extent(r_template) <- extent(for_cov)
	projection(r_template) <- CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0") 

	# ----------------------
	#  Rasterize existing protected areas according to this raster template.
	pa_sf <- pa_corr_sim %>%
		dplyr::mutate(pa_id = 1:n())%>%
		dplyr::select(-group)
	pa_sp <- as(st_union(st_buffer(pa_sf, 0.1), by_feature = TRUE), "Spatial")
	pa_r_tmp <- rasterize(pa_sp, r_template, field = pa_sp@data$size)
	
	
	
# =============================================================================
#  Functions 
# =============================================================================
	
	# ----------------------	
	#  Normalize raster to 0 to 1 scale.
	range01 <- function(x){(x-x_min)/(x_max - x_min)}
	#  Must always also set up min and max for each raster before using function, e.g.,:
	# x_min <- move_smooth@data@min
	# x_max <- move_smooth@data@max
	
	# ----------------------	
	#  Erase for sf objects
	st_erase = function(x, y) st_difference(x, st_union(st_combine(y)))

	# ----------------------	
	#  Function to determine if a vector output is empty.
	isEmpty <- function(x){
		return(length(x)==0)
		}
		
	# ----------------------
	# Function for calculation of relative probability of initial location.
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

	
	
# =============================================================================
#  Set up inital starting location and resistance matrices.
# =============================================================================

	# ----------------------
	#  Assign relative probability of initial occurrence to each grid cell from relative predicted 
	#   density using function assigned above (probsel()) in matrix format.
	#  Initial occurrence only from within protected areas, probability based on size.
	probrast <- pa_r_tmp/sum(getValues(pa_r_tmp), na.rm = T)
	x_min <- cellStats(probrast, min)
	x_max <- cellStats(probrast, max)
	probrast <- range01(probrast)
	
	# ----------------------
	#  Permeability: 
	resist_l <- c(0.75, 0.5, 0.3, 0.05, 0.05, 0.05, 0.05, 0.05, 0.5)
	resist_val <- as.data.frame(resist_l)
	#  Forest cover values:
	# NA - No data
	# 1 - 0 to 40 MT ACD/Non-forest Gaveau
	# 2  - 41 to 80 MT ACD/Regrowth forest Gaveau
	# 3 - 81 to 120 MT ACD/Logged forest Gaveau
	# 4 - 121 to 160 MT ACD
	# 5 - 161 to 200 MT ACD/Intact forest Gaveau
	# 6 - 201 to 240 MT ACD
	# 7 - 241 to 280 MT ACD
	# 8 - 281+ MT ACD
	# 9 - Roads
	
	
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
	#cl <- makeCluster(2)
	#registerDoParallel(cl)

	# ----------------------
	#  Initiatialize stack to hold output
	#out_stack_ls <- foreach(q = 1:nrep, .errorhandling="pass", .packages=c("dplyr", "sf", "sp", "raster", "SiMRiv", "foreach", "doParallel")) %dopar% {
		
		#move_sim_paths <- foreach(w = 1:nsim, .errorhandling="pass", .packages=c("dplyr", "sf", "sp", "raster", "SiMRiv")) %dopar% {
		#}
		
		
		out_ls <- list()
		
		for(w in 1:nsim){	
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
			first_pa <- as.numeric(as.character(st_intersects(first_loc, pa_sf)))
			all_pa <- as.numeric(as.character(st_intersects(sim_move_sf_tmp, pa_sf)))
		
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
			sim_move_sf <- st_as_sf(sim_move_df_fin, sf_column_name = "geometry")
			st_crs(sim_move_sf) = 32650
			sim_move_sf[is.na(sim_move_sf)] <- 0	
				
			out_ls[[w]] <- sim_move_sf
			}
		

		
		paths_bind <- do.call(rbind, out_ls) 
		paths <- paths_bind %>%
			dplyr::select(-all_pa, -start_pa, -cont_path) %>%
			group_by(sim_path) %>%
			mutate(lead = geometry[row_number() + 1],
				dist = st_distance(geometry, lead, by_element = T),
				path_dist = sum(dist)) %>%
			as.data.frame()
		
		paths_sf <- st_as_sf(paths, sf_column_name = "geometry")
		st_crs(paths_sf) = 32650
		
		start_point <- paths_sf %>%
			group_by(sim_path) %>%
			slice(1) %>%
			st_join(pa_sf) %>%
			as.data.frame() %>%			
			dplyr::select(sim_path, size, start_pa = pa_id)
			
		end_point <- paths_sf %>%
			group_by(sim_path) %>%
			slice(n()) %>%
			st_join(pa_sf) %>%
			as.data.frame() %>%
			dplyr::select(sim_path, size, end_pa = pa_id)
			
		paths_wt <- paths %>% #paths_sf %>%
			left_join(start_point, by = c("sim_path" = "sim_path")) %>%
			dplyr::rename(start_size = size) %>%
			left_join(end_point, by = c("sim_path" = "sim_path")) %>%
			dplyr::rename(end_size = size) %>%
			rowwise(.) %>%
			mutate(path_wt = (dist_wt * (1/log(path_dist))) + 
				(smaller_pa_wt * (log(min(start_size, end_size, na.rm = TRUE)))) + 
				(larger_pa_wt * (log(max(start_size, end_size,  na.rm = TRUE))))) %>%
			as.data.frame() %>%
			dplyr::select(-lead) %>%
			dplyr::filter(!is.na(end_pa)) %>%
			dplyr::filter(start_pa != end_pa) 
		
		paths_sf_wt <- st_as_sf(paths_wt, sf_column_name = "geometry")  #%>%
			#dplyr::select(path_wt)
		st_crs(paths_sf_wt) = 32650
		
		#save(paths_sf_wt, file = "paths_sf_wt.Rdata")
		#load(file = "paths_sf_wt5000.Rdata")
	
		large_pas <- c(1,3,4,5,9,39,50,52,57) # PA ids of PAs over 60000 ha
		paths_sf_wt <- paths_sf_wt %>% dplyr::filter(end_pa %in% large_pas)
				
		sim_ln <- paths_sf_wt %>%
			group_by(sim_path) %>%
			dplyr::summarise(do_union = FALSE) %>%
			st_cast("LINESTRING") %>%
			mutate(path_len = as.numeric(st_length(.))/1000) %>%
			as.data.frame()
		sim_ln <- st_as_sf(sim_ln, sf_column_name = "geometry") 
		paths_ln_wt <- as.data.frame(paths_sf_wt) %>%
			left_join(sim_ln, by = c("sim_path")) %>%
			group_by(sim_path) %>%
			slice(1) %>%
			rowwise(.) %>%
			mutate(path_wt_ln = (dist_wt * (1/log(path_len))) + 
				(smaller_pa_wt * (log(min(start_size, end_size, na.rm = TRUE)))) + 
				(larger_pa_wt * (log(max(start_size, end_size,  na.rm = TRUE))))) %>%
			dplyr::select(-dist, -geometry.x, -path_dist) %>%
			as.data.frame()
		paths_ln_sf_wt <- st_as_sf(paths_ln_wt, sf_column_name = "geometry.y")
		st_crs(paths_ln_sf_wt) = 32650
		paths_po_sf_wt <- st_buffer(paths_ln_sf_wt, 1000)
		
		corr_sim_r <- fasterize(paths_po_sf_wt, r_template, field = "path_wt_ln", fun="sum")
		plot(corr_sim_r)		
	
		m <- matrix(1, ncol=5, nrow=5)
		corr_sim_sm <- focal(corr_sim_m, m, fun="mean", na.rm=TRUE, pad=TRUE) 
		plot(corr_sim_sm)
	
		corr_sim_c <- raster::crop(corr_sim_sm, pu_r)
		corr_sim_rs <- raster::resample(corr_sim_c, pu_r)
		corr_sim_m_tmp <- raster::mask(corr_sim_rs, pu_r)
		corr_sim_m <- raster::mask(corr_sim_m_tmp, pa_sp, inverse = TRUE)
		plot(corr_sim_m)
		#writeRaster(corr_sim_m, file = "corr_sim_m.tif")
		
		
		
# =============================================================================
#  Plot using RasterVis.
# =============================================================================

	# ----------------------	
	#   Set color palette theme
	corr_sim_m[corr_sim_m == 0] <- NA
	r_theme <- rasterTheme(brewer.pal(9, "YlOrRd"))
	r_theme$panel.background$col = 'white'
	
	# ----------------------	
	plot_corr_sim <- levelplot(corr_sim_m, par.settings = r_theme, margin = FALSE) +
		layer(sp.polygons(pa_sp, lwd = 1, col = 'palegreen4', fill = 'palegreen4')) +
		layer(sp.polygons(main_sabah_sp, lwd = 0.8, col = 'grey40', fill = 'transparent')) 
	plot_corr_sim
	
	
	
# =============================================================================
#  Prioritization set up.
# =============================================================================

	# ----------------------
	#  Create raster following template of smaller study area (mainland Sabah)
	r_mat <- matrix(0, nrow(corr_sim_m), ncol(corr_sim_m))
	r_template <- raster(r_mat)
	extent(r_template) <- extent(corr_sim_m)
	projection(r_template) <- CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
	
	# ----------------------
	#  Create possible starting locations of corridors - can only be next to a PA.
	pa_sf <- pa_corr_sim %>%
		#dplyr::filter(size > 60000) %>%
		dplyr::mutate(pa_id = 1:n())%>%
		dplyr::select(-group)
		
	pa_bound_tmp <- st_buffer(pa_sf, 5000)
	pa_bound <- st_erase(pa_bound_tmp, st_buffer(pa_sf, 1))
	pu_start <- raster::mask(pu_r, pa_bound)
	pu_agg <- raster::aggregate(pu_r, fact = 3)
	pu_agg[pu_agg  > 0] <- 1
	pu_p <- rasterToPolygons(pu_r)
	
	#pa_sp <- as(pa_sf, "Spatial")
	pa_tst <- rasterize(st_buffer(pa_sf, 2000), r_template, field = pa_sp@data$pa_id)
	

	total_budget <- 200000
	single_pu_ha <-  as.numeric(st_area(st_as_sf(pu_p[1,])))/10000

	

# =============================================================================
#  Initiate lists to hold loop outputs
# =============================================================================

	# ----------------------
	#  Loop output holders
	whole_solution_out <- list()
	n_sol_steps <- vector("numeric", 50L)
	start_pa_id <- list()
	end_pa_id <- list()
	conn_pa <- vector("numeric", 15L)
	conn_start_pa <- vector("numeric", 15L)

	
	
# =============================================================================
#  Initiate prioritization problem and solve first iteration
# =============================================================================

		# ----------------------
		#  Initial problem
		p_init <- problem(x = pu_p, features = corr_sim_m, cost_column = "mainland_sabah_planning_units") %>%
			add_max_utility_objective(4) %>%
			add_feature_contiguity_constraints() %>%
			add_boundary_penalties(penalty = 0.1) %>%
			add_binary_decisions() %>%
			add_gurobi_solver(time_limit = 30)

		# ----------------------
		#  Solve initial problem.
		s_1 <- solve(p_init)
		
		# ----------------------
		#  Set up locked in areas (previous solution) for next round of selection.
		s_1_sf <- st_as_sf(s_1) %>% dplyr::filter(solution_1 == 1) %>% dplyr::select(solution_1)
		s_1_sp <- as(s_1_sf, "Spatial")
		
		plot(pu_p, col = "dodgerblue")
		plot(s_1_sp, add = T, col = "red")
		plot(pa_sp, add = T, col = "black")

		# ----------------------
		#  Check if the solution resulted in a full corridor from one planning unit to another.
		start_pu_sp <- as(st_buffer(s_1_sf, 1000), "Spatial")
		start_pa_tmp <- raster::extract(pa_tst, start_pu_sp) 
		start_pa <- unlist(lapply(start_pa_tmp, function(x) if (!is.null(x)) mean(x, na.rm=TRUE) else NA ))
		start_pa_id[[1]] <- ifelse(start_pa == "NaN", NA, start_pa)
		
	
		for(i in 2:50){
			
			# ----------------------
			#  Specify that planning unit selected in initiation or previous iteration is locked in.
			nam1 <- paste("s", i-1, "sp", sep = "_")
			nam1_sf <- paste("s", i-1, "sf", sep = "_")
			locked_in <- get(nam1)
			locked_in_area <- nrow(get(nam1_sf))
			
			nam2 <- paste("s", i-1, "sf", sep = "_")
			tmp <- get(nam2)
			new_corr <- raster::mask(corr_sim_m, st_buffer(tmp, 25000))
			new_pu <- raster::mask(pu_r, st_buffer(tmp, 25000))
			
			# ----------------------
			#  Set up problem to add planning units to first planning unit selected.
			#p <- problem(x = new_pu, features = new_corr) %>%
			p <-  problem(x = pu_p, features = new_corr, cost_column = "mainland_sabah_planning_units") %>%
				add_max_utility_objective(budget = locked_in_area + 1) %>%
				add_locked_in_constraints(locked_in) %>%
				add_feature_contiguity_constraints() %>%
				#add_boundary_penalties(penalty = 0.01) %>% #i =7, 12,14, 18, 27, 28, 31
				#add_boundary_penalties(penalty = 6.226575e-13) %>%
				add_binary_decisions() %>%
				add_gurobi_solver(time_limit = 30)

			# ----------------------
			#  Solve iteration problem.
			tmp3 <- solve(p)
			nam3 <- paste("s", i, sep = "_")
			assign(nam3, tmp3)
			
			# ----------------------
			#  Set up locked in areas (previous solution) for next round of selection.
			tmp4 <- st_as_sf(tmp3) %>% dplyr::filter(solution_1 == 1) %>% dplyr::select(solution_1)
			tmp5 <- as(tmp4, "Spatial")

			# ----------------------
			#  Assign new solution name.
			nam4 <- paste("s", i, "sf", sep = "_")
			assign(nam4, tmp4)
			nam5 <- paste("s", i, "sp", sep = "_")
			assign(nam5, tmp5)
	
			plot(pu_p, col = "dodgerblue")
			#plot(new_corr, add = TRUE)
			plot(tmp5, add = T, col = "red")
			plot(pa_sp, add = T, col = "black")
			
			# ----------------------
			#  Check if the solution resulted in a full corridor from one planning unit to another.
			#   Using start PA from first iteration of this scenario.
			start_pa_id[[i]] <- start_pa_id[[i-1]]
		
			end_pu_sp <- as(st_buffer(tmp4, 1000), "Spatial")
			end_pa_tmp <- raster::extract(pa_tst, end_pu_sp) 
			end_pa <- unlist(lapply(end_pa_tmp, function(x) if (!is.null(x)) mean(x, na.rm=TRUE) else NA ))
			end_pa_id[[i]] <- end_pa
		
			diff_pa_tst <- ifelse(isEmpty(setdiff(end_pa_id[[i]], start_pa_id[[i]])), 0, 
				setdiff(end_pa_id[[i]], start_pa_id[[i]])) 
			diff_pa_tst[is.nan(diff_pa_tst)] <- 0

			#if(any(diff_pa_tst > 0) == TRUE) {break}
			any(diff_pa_tst > 0)
			#if(cellStats(locked_in, sum) == cellStats(tmp3, sum)) {break}
			#cellStats(locked_in, sum) == cellStats(tmp3, sum)
			#if(corr_overlap > 100000) {break}
			
			}
			
		# stopped at i = 8 (one iteration past reaching two PAs)		

		j = 1
		n_sol_steps[j] <- i		
		
		conn_tst <- unlist(end_pa_tmp)
		conn_idx <- first(which((start_pa_id[[1]] - conn_tst) != 0))
		conn_pa[j] <- conn_tst[conn_idx]
		conn_start_pa[j] <- start_pa_id[[1]]
			
		last_sol <- get(paste("s", i, "sf", sep = "_"))
		last_sol_sf <- s_1_sf %>%
			mutate(sel_ord = 1) %>%
			mutate(corridor = j)
		 for(m in 2:n_sol_steps[1]){
			sol_nam <- paste("s", m, "sf", sep = "_")
			sol_tmp <- get(sol_nam) %>%
			mutate(sel_ord = m) %>%
			mutate(corridor = j)
			last_sol_sf <- rbind(last_sol_sf, sol_tmp)
			}
		whole_solution_out[[j]] <- last_sol_sf
		tot_area_h <- sum(as.numeric(st_area(last_sol_sf)))/10000
		


		#while(tot_area_h < 1000000){
			j = j + 1
	
		
			# ----------------------
			#  Re-initiate loop output holders
			start_pa_id_loop <- list()
			end_pa_id_loop <- list()
			
			# ----------------------
			#  Specify that planning units selected in initiation problem solution that should become 
			#  "locked out" and not available for selection.
			locked_out <- as(last_sol_sf, "Spatial")
			
			# ----------------------
			#  Initial problem
			p_init_loop <-  problem(x = pu_p, features = corr_sim_m, cost_column = "mainland_sabah_planning_units") %>%
				add_max_utility_objective(5) %>%
				add_locked_out_constraints(locked_out) %>%
				add_contiguity_constraints() %>%
				#add_boundary_penalties(penalty = 0.1) %>%
				add_binary_decisions() %>%
				add_gurobi_solver(time_limit = 30)

			# ----------------------
			#  Solve initial problem.
			s_1_loop <- solve(p_init_loop)
			
			# ----------------------
			#  Set up locked in areas (previous solution) for next round of selection.
			s_1_loop_sf <- st_as_sf(s_1_loop) %>% dplyr::filter(solution_1 == 1) %>% dplyr::select(solution_1)
			s_1_loop_sp <- as(s_1_loop_sf, "Spatial")
			
			plot(pu_p, col = "dodgerblue")
			plot(locked_out, add = TRUE, col = "magenta4")
			plot(s_1_loop_sp, add = T, col = "red")
			plot(pa_sp, add = T, col = "black")

			# ----------------------
			#  Check if the solution resulted in a full corridor from one planning unit to another.
			start_pu_sp <- as(st_buffer(s_1_loop_sf, 1000), "Spatial")
			start_pa_tmp <- raster::extract(pa_tst, start_pu_sp) 
			start_pa <- unlist(lapply(start_pa_tmp, function(x) if (!is.null(x)) mean(x, na.rm=TRUE) else NA ))
			start_pa_id_loop[[1]] <- ifelse(start_pa == "NaN", NA, start_pa)
			
			#for(i in 2:25){
				
				# ----------------------
				#  Specify that planning unit selected in initiation or previous iteration is locked in.
				nam1 <- paste("s", i-1, "loop", "sp", sep = "_")
				nam1_sf <- paste("s", i-1, "loop", "sf", sep = "_")
				locked_in <- get(nam1)
				locked_in_area <- nrow(get(nam1_sf))
				
				tmp <- get(nam1_sf)
				new_corr <- raster::mask(corr_sim_m, st_buffer(tmp, 25000))
				
				# ----------------------
				#  Set up problem to add planning units to first planning unit selected.
				p_loop <-  problem(x = pu_p, features = new_corr, cost_column = "mainland_sabah_planning_units") %>%
					add_max_utility_objective(budget = locked_in_area + 5) %>%
					add_locked_in_constraints(locked_in) %>%
					add_locked_out_constraints(locked_out) %>%
					add_contiguity_constraints() %>%
					#add_boundary_penalties(penalty = 0.05) %>%
					add_binary_decisions() %>%
					add_gurobi_solver(time_limit = 30)
	
				# ----------------------
				#  Solve iteration problem.
				tmp3 <- solve(p_loop)
				nam3 <- paste("s", i, "loop", sep = "_")
				assign(nam3, tmp3)
				
				# ----------------------
				#  Set up locked in areas (previous solution) for next round of selection.
				tmp4 <- st_as_sf(tmp3) %>% dplyr::filter(solution_1 == 1) %>% dplyr::select(solution_1)
				tmp5 <- as(tmp4, "Spatial")

				# ----------------------
				#  Assign new solution name.
				nam4 <- paste("s", i, "loop", "sf", sep = "_")
				assign(nam4, tmp4)
				nam5 <- paste("s", i, "loop", "sp", sep = "_")
				assign(nam5, tmp5)
	
				plot(pu_p, col = "dodgerblue")
				plot(locked_out, add = TRUE, col = "magenta4")
				plot(tmp5, add = T, col = "red")
				plot(pa_sp, add = T, col = "black")
			
				# ----------------------
				#  Check if the solution resulted in a full corridor from one planning unit to another.
				#   Using start PA from first iteration of this scenario.
				start_pa_id_loop[[i]] <- start_pa_id_loop[[i-1]]

				end_pu_sp <- as(st_buffer(tmp4, 1000), "Spatial")
				end_pa_tmp <- raster::extract(pa_tst, end_pu_sp) 
				end_pa <- unlist(lapply(end_pa_tmp, function(x) if (!is.null(x)) mean(x, na.rm=TRUE) else NA ))
				end_pa_id_loop[[i]] <- end_pa
	
				diff_pa_tst <- ifelse(isEmpty(setdiff(end_pa_id_loop[[i]], start_pa_id_loop[[i]])), 0, 
					setdiff(end_pa_id_loop[[i]], start_pa_id_loop[[i]])) 
				diff_pa_tst[is.nan(diff_pa_tst)] <- 0
				
				corr_overlap <- as.numeric(st_area(st_intersection(st_buffer(locked_out_sf, 10), st_buffer(tmp4, 10))))
				if(length(corr_overlap) == 0){
					corr_overlap <- 0
					} else {
						corr_overlap <- corr_overlap
						}
						
				#if(any(diff_pa_tst > 0) == TRUE) {break}
				#if(cellStats(locked_in, sum) == cellStats(tmp3, sum)) {break}
				#if(corr_overlap > 100000) {break}
			
				}
			
			# stop at i = 6 for corridor 2
			# stop at i = 8 for corridor 3
			# stop at i = 5 for corridor 4
			# stop at i = 2 for corridor 5
			# stop at i = 5 for corridor 6
			
			n_sol_steps[j] <- i		
			
			last_sol <- get(paste("s", i, "loop", "sf", sep = "_"))
			last_sol_sf <- s_1_loop_sf %>%
				mutate(sel_ord = 1) %>%
				mutate(corridor = j)
			 for(m in 2:n_sol_steps[j]){
				sol_nam <- paste("s", m, "loop", "sf", sep = "_")
				sol_tmp <- get(sol_nam) %>%
				mutate(sel_ord = m) %>%
				mutate(corridor = j)
				last_sol_sf <- rbind(last_sol_sf, sol_tmp)
				}
			whole_solution_out[[j]] <- last_sol_sf
			last_sol_sf <- rbind(st_as_sf(locked_out), last_sol_sf)
			tot_area_h <- sum(as.numeric(st_area(last_sol_sf)))/10000
	
			
			}

		all_sol <- do.call(rbind, whole_solution_out)

		# ----------------------
		#  Specify that planning unit selected in initiation or previous iteration is locked in.
		locked_in_sf <- distinct(all_sol)
		locked_in <- as(locked_in_sf, "Spatial")
		new_corr <- raster::mask(corr_sim_m, st_buffer(all_sol, 25000))
		locked_in_area <- sum(as.numeric(st_area(locked_in_sf)))/10000
		rem_pu <- round((total_budget - locked_in_area)/single_pu_ha)
		
		# ----------------------
		#  Set up problem to add planning units to first planning unit selected.
		p_last <-  problem(x = pu_p, features = corr_sim_m, cost_column = "mainland_sabah_planning_units") %>%
			add_max_utility_objective(budget = nrow(locked_in_sf) + rem_pu) %>%
			add_locked_in_constraints(locked_in) %>%
			add_boundary_penalties(penalty  = 0.005) %>%
			add_binary_decisions() %>%
			add_gurobi_solver(time_limit = 30)

		# ----------------------
		#  Solve iteration problem.
		s_last <- solve(p_last)
	
		# ----------------------
		#  Set up locked in areas (previous solution) for next round of selection.
		s_last_sf <- st_as_sf(s_last) %>% dplyr::filter(solution_1 == 1) %>% dplyr::select(solution_1)
		s_last_sp <- as(s_last_sf, "Spatial")

		plot(pu_p, col = "dodgerblue")
		plot(s_last_sp, add = T, col = "red")
		plot(locked_in, add = T, col = "orange")
		plot(pa_sp, add = T, col = "black")
	
		
		
		
		all_sol_xy <- all_sol %>%
			st_centroid() %>%
			st_coordinates() %>%
			as.data.frame()
			
		
		all_sol_arr <- all_sol %>%
			cbind(all_sol_xy$X) %>%
			cbind(all_sol_xy$Y) %>%
			arrange(corridor, sel_ord) %>%
			group_by(all_sol_xy.X, all_sol_xy.Y) %>%
			slice(1) %>%
			as.data.frame() %>%
			dplyr::select(-all_sol_xy.X, -all_sol_xy.Y, -solution_1) %>%
			st_as_sf(sf_column_name = "geometry")
	
		max_corr <- max(all_sol_arr$corridor)

		complete_sol <- s_last_sf %>%
			st_erase(all_sol_arr) %>%
			dplyr::select(-solution_1) %>%
			mutate(sel_ord = 1) %>%
			mutate(corridor = max_corr + 1) %>%
			rbind(all_sol_arr) %>%
			arrange(corridor, sel_ord) %>%
			mutate(all_order = 1:n())



	colr <- RColorBrewer::brewer.pal(max(complete_sol$corridor) + 1, rev("YlOrRd"))
	colr <- colr[-1]
	
	p <- ggplot() +
			geom_sf(data = main_sabah_sf, colour = "grey80", fill = "grey80", size = 1) +
			geom_sf(data = pa_sf, colour = "darkseagreen4", fill = "darkseagreen4") +
			xlim(300000, 800000) +
			ylim(400000, 800000) +
			theme_bw() 
			
		nam <- paste("plots/p_", 0, ".png", sep = "")	
		ggsave(nam, p, dpi = 300)
	
	for(g in 1:max(complete_sol$corridor)){
	
		dat1 <- complete_sol %>%
			dplyr::filter(corridor == g) 
			
		for(i in 1:length(unique(dat1$sel_ord))){
		
			dat2 <- dat1 %>%
				dplyr::filter(sel_ord == i)
			
			tmp_colr <- colr[dat2[1,]$corridor]
			
			p <- p + geom_sf(data = dat2, colour = tmp_colr, fill = tmp_colr) 
			
			nam <- paste("plots/p_", g, i, ".png", sep = "")	
			ggsave(nam, p, dpi = 300)
		
			}
		}








		
# =============================================================================	
###############################################################################



	'%!in%' <- function(x,y)!('%in%'(x,y))

###############################################################################
#  Do not allow corridors to start and end at same PAs as a
#   previous round
###############################################################################




# =============================================================================
#  Initiate lists to hold loop outputs
# =============================================================================

	# ----------------------
	#  Loop output holders
	whole_solution_out <- list()
	n_sol_steps <- vector("numeric", 50L)
	start_pa_id <- list()
	end_pa_id <- list()
	#conn_pa <- vector("numeric", 15L)
	conn_start_pa <- vector("numeric", 15L)

	# ----------------------
	#  Create possible starting locations of corridors - can only be next to a PA.
	pa_sf <- pa_corr_sim %>%
		#dplyr::filter(size > 60000) %>%
		dplyr::mutate(pa_id = 1:n())%>%
		dplyr::select(-group)
		
	pa_bound_tmp <- st_buffer(pa_sf, 5000)
	pa_bound <- st_erase(pa_bound_tmp, st_buffer(pa_sf, 1))
	pa_bound_sp <- as(pa_bound, "Spatial")
	pu_start <- raster::mask(pu_r, pa_bound_sp)

	
	pa_sp <- as(pa_sf, "Spatial")
	pa_tst <- rasterize(pa_sp, r_template, field = pa_sp@data$pa_id)
	
	
	
# =============================================================================
#  Initiate prioritization problem and solve first iteration
# =============================================================================

		# ----------------------
		#  Initial problem
		p_init <- problem(x = pu_start, features = corr_sim_m) %>%
			add_max_utility_objective(5005) %>%
			add_contiguity_constraints() %>%
			add_binary_decisions() %>%
			add_gurobi_solver(time_limit = 600)


		# ----------------------
		#  Solve initial problem.
		s_1 <- solve(p_init)
		plot(pu_r, col = "dodgerblue")
		plot(s_1, add = T, col = c("transparent", "red"))
									
		# ----------------------
		#  Set up locked in areas (previous solution) for next round of selection.
		s_1_sp <- rasterToPolygons(s_1, fun = function(x){x>0}, dissolve = TRUE)
		s_1_sf <- st_as_sf(s_1_sp)
		
		# ----------------------
		#  Check if the solution resulted in a full corridor from one planning unit to another.
		start_pu_sp <- as(st_buffer(s_1_sf, 1000), "Spatial")
		start_pa_tmp <- raster::extract(pa_tst, start_pu_sp) 
		start_pa <- unlist(lapply(start_pa_tmp, function(x) if (!is.null(x)) mean(x, na.rm=TRUE) else NA ))
		start_pa_id[[1]] <- ifelse(start_pa == "NaN", NA, start_pa)
		
		
		for(i in 2:150){
				
			# ----------------------
			#  Specify that planning unit selected in initiation or previous iteration is locked in.
			nam1 <- paste("s", i-1, sep = "_")
			locked_in <- get(nam1)
			locked_in_area <- cellStats(locked_in, sum) * 100
			
			# ----------------------
			#  Set up problem to add planning units to first planning unit selected.
			p <- problem(x = pu_r, features = corr_sim_m) %>%
				add_max_utility_objective(budget = locked_in_area + 5001) %>%
				add_locked_in_constraints(locked_in) %>%
				#add_contiguity_constraints() %>%
				add_boundary_penalties(penalty = 50) %>%
				add_binary_decisions() %>%
				add_gurobi_solver(time_limit = 600)

			
			# ----------------------
			#  Solve iteration problem.
			tmp3 <- solve(p)
			nam3 <- paste("s", i, sep = "_")
			assign(nam3, tmp3)
			
			print(cellStats(locked_in, sum))
			print(cellStats(tmp3, sum))
			plot(pu_r, col = "dodgerblue")
			plot(tmp3, add = T, col = c("transparent", "red"))

			# ----------------------
			#  Set up locked in areas (previous solution) for next round of selection.
			tmp4 <- rasterToPolygons(tmp3, fun = function(x){x>0}, dissolve = TRUE)
			tmp5 <- st_as_sf(tmp4)
				
			# ----------------------
			#  Assign new solution name.
			nam4 <- paste("s", i, "sp", sep = "_")
			assign(nam4, tmp4)
			nam5 <- paste("s", i, "sf", sep = "_")
			assign(nam5, tmp5)

			# ----------------------
			#  Check if the solution resulted in a full corridor from one planning unit to another.
			#   Using start PA from first iteration of this scenario.
			start_pa_id[[i]] <- start_pa_id[[i-1]]
							
			end_pu_sp <- as(st_buffer(tmp5, 1000), "Spatial")
			end_pa_tmp <- raster::extract(pa_tst, end_pu_sp) 
			end_pa <- unlist(lapply(end_pa_tmp, function(x) if (!is.null(x)) mean(x, na.rm=TRUE) else NA ))
			end_pa_id[[i]] <- end_pa
																
			diff_pa_tst <- ifelse(isEmpty(setdiff(end_pa_id[[i]], start_pa_id[[i]])), 0, 
				setdiff(end_pa_id[[i]], start_pa_id[[i]])) 
			diff_pa_tst[is.nan(diff_pa_tst)] <- 0
							
			if(any(diff_pa_tst > 0) == TRUE) {break}
			if(cellStats(locked_in, sum) == cellStats(tmp3, sum)) {break}
			}
			
			i
			
			for(j in i+1:150){
				
				# ----------------------
				#  Specify that planning unit selected in initiation or previous iteration is locked in.
				nam1 <- paste("s", j-1, sep = "_")
				locked_in <- get(nam1)
				locked_in_area <- cellStats(locked_in, sum) * 100
				
				# ----------------------
				#  Set up problem to add planning units to first planning unit selected.
				p <- problem(x = pu_r, features = corr_sim_m) %>%
					add_max_utility_objective(budget = locked_in_area + 5001) %>%
					add_locked_in_constraints(locked_in) %>%
					#add_contiguity_constraints() %>%
					add_boundary_penalties(penalty = 50) %>%
					add_binary_decisions() %>%
					add_gurobi_solver(time_limit = 600)
	
				
				# ----------------------
				#  Solve iteration problem.
				tmp3 <- solve(p)
				nam3 <- paste("s", j, sep = "_")
				assign(nam3, tmp3)
				
				print(cellStats(locked_in, sum))
				print(cellStats(tmp3, sum))
				plot(pu_r, col = "dodgerblue")
				plot(tmp3, add = T, col = c("transparent", "red"))
	
				# ----------------------
				#  Set up locked in areas (previous solution) for next round of selection.
				tmp4 <- rasterToPolygons(tmp3, fun = function(x){x>0}, dissolve = TRUE)
				tmp5 <- st_as_sf(tmp4)
					
				# ----------------------
				#  Assign new solution name.
				nam4 <- paste("s", j, "sp", sep = "_")
				assign(nam4, tmp4)
				nam5 <- paste("s", j, "sf", sep = "_")
				assign(nam5, tmp5)
	
				# ----------------------
				#  Check if the solution resulted in a full corridor from one planning unit to another.
				#   Using start PA from first iteration of this scenario.
				start_pa_id[[i]] <- start_pa_id[[j-1]]
								
				end_pu_sp <- as(st_buffer(tmp5, 1000), "Spatial")
				end_pa_tmp <- raster::extract(pa_tst, end_pu_sp) 
				end_pa <- unlist(lapply(end_pa_tmp, function(x) if (!is.null(x)) mean(x, na.rm=TRUE) else NA ))
				end_pa_id[[j]] <- end_pa
																	
				diff_pa_tst <- ifelse(isEmpty(setdiff(end_pa_id[[j]], start_pa_id[[j]])), 0, 
					setdiff(end_pa_id[[j]], start_pa_id[[j]])) 
				diff_pa_tst[is.nan(diff_pa_tst)] <- 0
								
				if(any(diff_pa_tst > 0) == TRUE) {break}
				if(cellStats(locked_in, sum) == cellStats(tmp3, sum)) {break}
				}
				
			
			
			
			
			
			
			
			
			
		j = 1
		n_sol_steps[j] <- i		
		
		#conn_tst <- unlist(end_pa_tmp)
		#conn_idx <- first(which((start_pa_id[[1]] - conn_tst) != 0))
		#conn_pa[j] <- conn_tst[conn_idx]
		if(cellStats(locked_in, sum) == cellStats(tmp3, sum)) {
			conn_start_pa[j] <- 0
			} else {
			conn_start_pa[j] <- start_pa_id[[1]]
			}
			
		last_sol <- get(paste("s", i, sep = "_"))
		last_sol_sf <- s_1_sf %>%
			mutate(sel_ord = 1) %>%
			mutate(corridor = j)
		 for(m in 2:n_sol_steps[1]){
			sol_nam <- paste("s", m, "sf", sep = "_")
			sol_tmp <- get(sol_nam) %>%
			mutate(sel_ord = m) %>%
			mutate(corridor = j)
			last_sol_sf <- rbind(last_sol_sf, sol_tmp)
			}
		whole_solution_out[[j]] <- last_sol_sf
		tot_area_h <- cellStats(last_sol, sum) * 100
		





		#for(j in 2:20){
		#while(tot_area_h < 1000000){
		j = j + 1
		
			pa_start_new_sf <- pa_corr_sim %>%
				#dplyr::filter(size > 60000) %>%
				dplyr::mutate(pa_id = 1:n())%>%
				dplyr::select(-group) %>%
				dplyr::filter(pa_id %!in% conn_start_pa)
				
			pa_bound_tmp <- st_buffer(pa_start_new_sf, 5000)
			pa_bound <- st_erase(pa_bound_tmp, st_buffer(pa_start_new_sf, 1))
			pa_bound_sp <- as(pa_bound, "Spatial")
			pu_start_new <- raster::mask(pu_r, pa_bound_sp)
		

			# ----------------------
			#  Re-initiate loop output holders
			start_pa_id_loop <- list()
			end_pa_id_loop <- list()
			i = 1
			
			# ----------------------
			#  Specify that planning units selected in initiation problem solution that should become 
			#  "locked out" and not available for selection.
			locked_out_start <- raster::mask(last_sol, pu_start_new)
			locked_out <- last_sol
			locked_out_h <- cellStats(locked_out_start, sum, nam.rm = TRUE)
			locked_out_sp <- rasterToPolygons(locked_out, fun = function(x){x>0}, dissolve = TRUE)
			locked_out_sf <- st_as_sf(locked_out_sp)
			
		
			# ----------------------
			#  Initial problem
			if(locked_out_h < 1){
				
				p_init_loop <- problem(x = pu_start_new, features = corr_sim_m) %>%
					add_max_utility_objective(5005) %>%
					#add_locked_out_constraints(locked_out_start) %>%
					add_contiguity_constraints() %>%
					add_binary_decisions() %>%
					add_gurobi_solver(time_limit = 600)
				} else {
				
				p_init_loop <- problem(x = pu_start_new, features = corr_sim_m) %>%
					add_max_utility_objective(5005) %>%
					add_locked_out_constraints(locked_out_start) %>%
					add_contiguity_constraints() %>%
					add_binary_decisions() %>%
					add_gurobi_solver(time_limit = 600)
				
				}
	
			# ----------------------
			#  Solve initial problem.
			s_1_loop <- solve(p_init_loop)
			plot(pu_r, col = "dodgerblue")
			plot(pa_sp, add = T)
			plot(locked_out, add = T, col = c("transparent", "purple"))
			plot(s_1_loop, add = T, col = c("transparent", "red"))
										
			# ----------------------
			#  Set up locked in areas (previous solution) for next round of selection.
			s_1_sp_loop <- rasterToPolygons(s_1_loop, fun = function(x){x>0}, dissolve = TRUE)
			s_1_sf_loop <- st_as_sf(s_1_sp_loop)
			
			# ----------------------
			#  Check if the solution resulted in a full corridor from one planning unit to another.
			start_pu_sp <- as(st_buffer(s_1_sf_loop, 1000), "Spatial")
			start_pa_tmp <- raster::extract(pa_tst, start_pu_sp) 
			start_pa <- unlist(lapply(start_pa_tmp, function(x) if (!is.null(x)) mean(x, na.rm=TRUE) else NA ))
			start_pa_id_loop[[1]] <- ifelse(start_pa == "NaN", NA, start_pa)
			
			corr_overlap <- as.numeric(st_area(st_intersection(st_buffer(locked_out_sf, 10), st_buffer(s_1_sf_loop, 10))))
			if(length(corr_overlap) == 0){
				corr_overlap <- 0
				} else {
				corr_overlap <- corr_overlap
				}
			
			if(corr_overlap > 10000){
				print("STOP")
				
				n_sol_steps[j] <- i		
			
				conn_start_pa[j] <- start_pa_id_loop[[1]]
				
				last_sol <- get(paste("s", i, "loop", sep = "_"))
				last_sol_sf <- s_1_sf_loop %>%
					mutate(sel_ord = 1) %>%
					mutate(corridor = j)
				
				whole_solution_out[[j]] <- last_sol_sf
				
				last_sol <- stack(last_sol, locked_out)
				last_sol <- raster::calc(last_sol, sum, na.rm = TRUE)
				
				tot_area_h <- cellStats(last_sol, sum) * 100
				
				
				} else {
					print("CONTINUE")
				}
					
					
					
			for(i in 2:50){
					
				# ----------------------
				#  Specify that planning unit selected in initiation or previous iteration is locked in.
				nam1 <- paste("s", i-1, "loop", sep = "_")
				locked_in <-get(nam1)
				locked_in_area <- cellStats(locked_in, sum) * 100
				
				# ----------------------
				#  Set up problem to add planning units to first planning unit selected.
				p_loop <- problem(x = pu_r, features = corr_sim_m) %>%
					add_max_utility_objective(budget = locked_in_area + 5001) %>%
					add_locked_out_constraints(locked_out) %>%
					add_locked_in_constraints(locked_in) %>%
					#add_contiguity_constraints() %>%
					add_boundary_penalties(penalty = 50) %>%
					add_binary_decisions() %>%
					add_gurobi_solver(time_limit = 600)
				
				# ----------------------
				#  Solve iteration problem.
				tmp3 <- solve(p_loop)
				nam3 <- paste("s", i, "loop", sep = "_")
				assign(nam3, tmp3)
				
				cellStats(locked_in, sum)
				cellStats(tmp3, sum)
				plot(pu_r, col = "dodgerblue")
				plot(pa_sp, add = T)
				plot(locked_out, add = T, col = c("transparent", "purple"))
				plot(tmp3, add = T, col = c("transparent", "red"))
				
	
				# ----------------------
				#  Set up locked in areas (previous solution) for next round of selection.
				tmp4 <- rasterToPolygons(tmp3, fun = function(x){x>0}, dissolve = TRUE)
				tmp5 <- st_as_sf(tmp4)
					
				# ----------------------
				#  Assign new solution name.
				nam4 <- paste("s", i, "sp", "loop", sep = "_")
				assign(nam4, tmp4)
				nam5 <- paste("s", i, "sf", "loop", sep = "_")
				assign(nam5, tmp5)
	
				# ----------------------
				#  Check if the solution resulted in a full corridor from one planning unit to another.
				#   Using start PA from first iteration of this scenario.
				start_pa_id_loop[[i]] <- start_pa_id_loop[[i-1]]
								
				end_pu_sp <- as(st_buffer(tmp5, 1000), "Spatial")
				end_pa_tmp <- raster::extract(pa_tst, end_pu_sp) 
				end_pa <- unlist(lapply(end_pa_tmp, function(x) if (!is.null(x)) mean(x, na.rm=TRUE) else NA ))
				end_pa_id_loop[[i]] <- end_pa
								
				diff_pa_tst <- ifelse(isEmpty(setdiff(end_pa_id_loop[[i]], start_pa_id_loop[[i]])), 0, 
					setdiff(end_pa_id_loop[[i]], start_pa_id_loop[[i]])) 
				diff_pa_tst[is.nan(diff_pa_tst)] <- 0
				
				corr_overlap <- as.numeric(st_area(st_intersection(st_buffer(locked_out_sf, 10), st_buffer(tmp5, 10))))
				if(length(corr_overlap) == 0){
					corr_overlap <- 0
					} else {
					corr_overlap <- corr_overlap
					}			
								
				if(any(diff_pa_tst > 0) == TRUE) {break}
				if(cellStats(locked_in, sum) == cellStats(tmp3, sum)) {break}
				if(corr_overlap > 10000) {break}
				}
				
				
				n_sol_steps[j] <- i		
				
				#conn_tst <- unlist(end_pa_tmp)
				#conn_idx <- first(which((start_pa_id_loop[[1]] - conn_tst) != 0))
				#conn_pa[j] <- conn_tst[conn_idx]
				if(cellStats(locked_in, sum) == cellStats(tmp3, sum)) {
					conn_start_pa[j] <- 0
					} else {
					conn_start_pa[j] <- start_pa_id_loop[[1]]
					}
			
				#conn_start_pa[j] <- start_pa_id_loop[[1]]
				
				last_sol <- get(paste("s", i, "loop", sep = "_"))
				last_sol_sf <- s_1_sf_loop %>%
					mutate(sel_ord = 1) %>%
					mutate(corridor = j)
				 for(m in 2:n_sol_steps[j]){
					sol_nam <- paste("s", m, "sf","loop", sep = "_")
					sol_tmp <- get(sol_nam) %>%
					mutate(sel_ord = m) %>%
					mutate(corridor = j)
					last_sol_sf <- rbind(last_sol_sf, sol_tmp)
					}
				whole_solution_out[[j]] <- last_sol_sf
				
				last_sol <- stack(last_sol, locked_out)
				last_sol <- raster::calc(last_sol, sum, na.rm = TRUE)
				
				tot_area_h <- cellStats(last_sol, sum) * 100
			#}



			all_sol <- do.call(rbind, whole_solution_out)


r <- stack()
for(i in 1:nrow(all_sol_arr)){
	sf_tmp <- all_sol_arr[i,]
	r_tmp <- fasterize(sf_tmp, r_template, field = "tot_ord")
	r <- stack(r, r_tmp)
	}
r_sum <- raster::calc(r, sum, na.rm = TRUE)
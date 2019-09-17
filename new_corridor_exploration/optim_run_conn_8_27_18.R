###############################################################################
#  Run optimization of connectivity input iteratively with 'Prioritizr' package.
#  September 29, 2017; last updated April 9, 2018
#  Script to generate optimize selected land units for protection of a variety of input 
#   conservation features.
#  Sara Williams
###############################################################################



# =============================================================================
#  Notes on model implementation, assumptions, and areas for improvement.
# =============================================================================

	# ----------------------
	#  Inputs:
	#   1. Conservation *features*. These may be rasters or spatial polygons. I am thinking I will use rasters, where 
	#   the extent is the entirety of Sabah, and each cell value is a 1 or 0, for if the spp range covers that cell or not. These 
	#   will take the format of a RasterStack object with nLayers = number of conservation features.
	#   2. Planning units. This may be raster or polygon. If raster, then each pixel value is the cost information. If polygon, 
	#   than should have 3 variables: cost, locked_in, locked_out.

	
	
# =============================================================================
#  Load packages.
# =============================================================================
library(sf)
library(sp)
library(raster)
library(dplyr)
library(tidyr)
# install.packages("C:/gurobi801/win64/R/gurobi_8.0-1.zip", repos = NULL)
library(gurobi)
library(units)
#devtools::install_github("prioritizr/priortizr")
library(prioritizr)
library(ggplot2)
#library(rasterVis)
#library(beepr)

#options(digits = 10)

setwd("C:/Users/saraw/Documents/Prioritization/")



# =============================================================================
#  Load data.
# =============================================================================		

	# ----------------------
	#  Boundaries
	load(file = "study_area_boundaries_and_pas/main_sabah_sf.Rdata")
	load(file = "study_area_boundaries_and_pas/main_ssk_sf.Rdata")
	main_ssk_sp <- as(main_ssk_sf, "Spatial")
	main_sabah_sp <- as(main_sabah_sf, "Spatial")


	
	# ----------------------
	# Existing protected areas prepared for movement simulation corridor prioritization.
	load(file = "conn_move_sim/pa_sf_move_sim.Rdata")
	#load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/ssk_pa_near_sf_clust.Rdata")
	
	# ----------------------
	#  Load study area grid that holds planning units for connectivity optimization.
	load(file = "planning_unit_grids/sa_grid.Rdata")
	
	# ----------------------
	#  Load connectivity information from movemnt simulation and weighting of shortest paths
	#  between protected areas.
	moves_sm <- raster("C:/Users/saraw/Documents/SEARRP_Analyses/conn_move_sim/moves_sm.grd")
	
	
	
# =============================================================================
#  Generate needed spatial layers and functions needed for assessment.
# =============================================================================

	# ----------------------
	#  Create raster following template of study area (cell values are empty)
	r_mat <- matrix(0, nrow(moves_sm), ncol(moves_sm))
	r_template <- raster(r_mat)
	extent(r_template) <- extent(moves_sm)
	projection(r_template) <- CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0") 

	# ----------------------
	#  Rasterize existing protected areas according to this raster template.
	pa_sf <- pa_sf_move_sim
	# pa_sf <- ssk_pa_near_sf_clust %>%
		# mutate(area_tmp = as.numeric(st_area(.) * 0.0001)) %>%
		# dplyr::filter(area_tmp > 10000) %>%
		# mutate(id = 1:n())
	# pa_buff_sf <- st_buffer(pa_sf, 1100)
	# pa_sp <- as(pa_sf, "Spatial")
	# pa_buff_sp <- as(pa_buff_sf, "Spatial")
	# pa_r <- rasterize(pa_buff_sp, r_template, field = pa_buff_sp$id)
	
	# ----------------------	
	#  Function to determine if a vector output is empty.
	isEmpty <- function(x){
		return(length(x)==0)
		}
	st_erase = function(x, y) st_difference(x, st_union(st_combine(y)))
	
	# ----------------------	
	#  Function to cluster nearby units..
	clusterSF <- function(sfpolys, thresh){

		dmat = st_distance(sfpolys)
		hc = hclust(as.dist(dmat>thresh), method="single")
		groups = cutree(hc, h=0.5)
		
		d = st_sf(
			geom = do.call(c,
				lapply(1:max(groups), function(g){
					st_union(sfpolys[groups==g,])
				})
				)
			)
		d$group = 1:nrow(d)
		d
	}
	
	# ----------------------	
	#  Function to normalize rasters to 0 to 1 scale.
	range01 <- function(x){(x-x_min)/(x_max - x_min)}
	
	
	
# =============================================================================
#  Define planning units.
# =============================================================================

	# ----------------------
	#  Set up area of each planning unit as cost.
	pu_sf <- sa_grid %>%
		mutate(area_h = 500)
	pu_in <- as(pu_sf, "Spatial")
	
	pa_bound <- st_buffer(st_boundary(pa_sf), 1000)
	pa_grid_df_tmp <- as.data.frame(st_intersects(pu_sf, pa_bound, sparse = FALSE))
	
	pa_grid_df_tmp <- pa_grid_df_tmp*1
	pa_grid_df <- pa_grid_df_tmp %>%
		mutate(IN = rowSums(.[1:ncol(pa_grid_df_tmp)])) %>%
		dplyr::select(IN)
	pu_start <- pu_sf %>%
		bind_cols(pa_grid_df) %>%
		dplyr::filter(IN == 1)
	pu_start_in <- as(pu_start, "Spatial")
	
	
	
# =============================================================================
#  Parameters for loop to solve prioritized connectivity iteratively.
# =============================================================================

	# ----------------------
	#  Average planning unit area
	pu_area_h <- 500
	
	# ----------------------
	#  Total area allocated
	tot_h <- 400000
	
	# ----------------------
	#  Number of iterations to run when you still haven't reached a 
	#   different PA.
	#niter_max_corr <- 55
	niter_max_corr <- 60
	
	# ----------------------
	#  Maximum number of iterations to run in total
	niter_max_tot <- tot_h/pu_area_h
	
	# ----------------------
	#  Number of neighbors for each iteration
	n_neigh <- 2
	
	# ----------------------
	#  Number of neighbors for initial problem.
	n_neigh_init <- 1
	
	# ----------------------
	#  Adjustment area for adding wiggle room to cover area increase.
	wiggle_room <- 100
	
	# ----------------------
	#  Area of coverage for initial problems.
	init_area <- pu_area_h * (n_neigh_init + 3) + wiggle_room
	#init_area <- pu_area_h  + wiggle_room
	
	# ----------------------
	#  Area of increase for each iteration: must put this formula in the problem 
	# cur_area_h + (pu_area_h * (n_neigh + 3)) + wiggle_room

	# ----------------------
	#  Targets for for initial problems
	#init_target <- 0.01


	
# =============================================================================
#  Initiate lists to hold loop outputs
# =============================================================================

	# ----------------------
	#  Loop output holders
	whole_solution_out <- list()
	n_sol_steps <- vector("numeric", 50L)
	start_pa_id <- list()
	end_pa_id <- list()

	
	
# =============================================================================
#  Initiate problem and solve first iteration
# =============================================================================

		# ----------------------
		#  Set up initial problem.
		p_init <- problem(x = pu_start_in, features = moves_sm, cost_column = "area_h") %>%
			add_max_utility_objective(505) %>%
			add_binary_decisions() %>%
			add_gurobi_solver()

		# ----------------------
		#  Solve initial problem.
		s_1 <- solve(p_init)
		s_1_sf <- st_as_sf(s_1) %>%
			dplyr::filter(solution_1 == 1) %>%
			mutate(locked_in = 1) %>%
			mutate(select_order = 1)
		cur_area_h <- as.integer(sum(st_area(s_1_sf)) * 0.0001)
		
		# ----------------------
		#  Make solution an sp object.
		s_1_sp <- as(s_1_sf, "Spatial")
		
		# ----------------------
		#  Buffer around solution to make a boundary of "local" area for next iteration of 
		#   optimization.
		bound <- st_buffer(st_union(s_1_sf), 5000) 
		bound_sp <- as(bound, "Spatial")
		
		# ----------------------
		#  Constrain input raster for the next iteration to the same general area as the boundary
		#   generated above.
		local_r <- raster::mask(moves_sm, bound_sp)
		conn_in_upd_1_01 <- local_r
		
		# ----------------------
		#  Constrain planning units for the next iteration to the same general area as the features
		#   generated above.
		pu_upd_1_tmp <- raster::crop(pu_in, bound_sp)
		pu_upd_1_sf <- st_as_sf(pu_upd_1_tmp) %>%
			mutate(size = as.integer(st_area(.)*0.0001)) %>%
			dplyr::filter(size > 490)
		pu_upd_1 <- as(pu_upd_1_sf, "Spatial")
		
		# ----------------------
		#  Check if the solution resulted in a full corridor from one planning unit to another.
		start_pu_sp <- s_1_sp
		start_pa_tmp <- raster::extract(pa_r, start_pu_sp) 
		start_pa <- unlist(lapply(start_pa_tmp, function(x) if (!is.null(x)) mean(x, na.rm=TRUE) else NA ))
		start_pa_id[[1]] <- ifelse(start_pa == "NaN", NA, start_pa)
		#start_pa_id[[1]] <- start_pa_tmp[[1]][1]
		
		end_pu_sp <- s_1_sp
		end_pa_tmp <- raster::extract(pa_r, end_pu_sp) 
		end_pa <- unlist(lapply(end_pa_tmp, function(x) if (!is.null(x)) mean(x, na.rm=TRUE) else NA ))
		end_pa_id[[1]] <- end_pa
		
		diff_pa_tst <- ifelse(isEmpty(setdiff(end_pa_id[[1]], start_pa_id[[1]])), 0, 
			setdiff(end_pa_id[[1]], start_pa_id[[1]])) 
		diff_pa_tst[is.nan(diff_pa_tst)] <- 0
		
		if(any(diff_pa_tst > 0) == FALSE){
		
			for(i in 2:niter_max_corr){
			#for(i in 2:20){
				# ----------------------
				#  Specify that planning unit selected in initiation or previous iteration is locked in.
				nam1 <- paste("s", i-1, "sf", sep = "_")
				tmp1 <- get(nam1)
				locked_in <- as(tmp1, "Spatial")

				nam_upd_feat <- paste("conn_in_upd", i-1, "01", sep = "_")
				upd_feat_tmp <- get(nam_upd_feat)
				
				nam_upd_pu <- paste("pu_upd", i-1, sep = "_")
				upd_pu_tmp <- get(nam_upd_pu)
				
				# ----------------------
				#  Set up problem to add planning units to first planning unit selected.
				p <- problem(x = upd_pu_tmp, features = upd_feat_tmp, cost_column = "area_h") %>%
					add_max_utility_objective(505*i) %>%
					#add_boundary_penalties(0.00002, 0.5) %>%
					add_binary_decisions() %>%
					add_gurobi_solver()
				
				# ----------------------
				#  Solve iteration problem.
				tmp3 <- solve(p)
				tmp4 <- st_as_sf(tmp3) %>%
					dplyr::filter(solution_1 == 1) %>%
					mutate(locked_in = 1) %>%
					mutate(select_order = i)
				
				# ----------------------
				#  Combine with previous solution and remove duplicates ---
				#   used to obtain correct selection order number.
				tmp1_df <- as.data.frame(tmp1)
				tmp1_coords <- as.data.frame(st_coordinates(st_cast(tmp1, "MULTIPOLYGON"))) %>%
					group_by_at(ncol(.)) %>%
					slice(1) %>%
					as.data.frame(.) %>%
					mutate(X_coord = round(X), Y_coord = round(Y)) 
				tmp1_dup <- cbind(tmp1_df, tmp1_coords) %>%
					dplyr::select(area_h, solution_1, locked_in, select_order,
					geometry, X_coord, Y_coord) 
				
				tmp4_df <- as.data.frame(tmp4)
				tmp4_coords <- as.data.frame(st_coordinates(st_cast(tmp4, "MULTIPOLYGON"))) %>%
					group_by_at(ncol(.)) %>%
					slice(1) %>%
					as.data.frame(.) %>%
					mutate(X_coord = round(X), Y_coord = round(Y)) 
				tmp4_dup <- cbind(tmp4_df, tmp4_coords) %>%
					dplyr::select(area_h, solution_1, locked_in, select_order,
					geometry, X_coord, Y_coord) 
				
				tmp14 <- rbind(tmp1_dup, tmp4_dup) %>%
					arrange(X_coord, desc(Y_coord), select_order) %>%
					group_by(X_coord, Y_coord) %>%
					slice(1) %>%
					as.data.frame() %>%
					select(-X_coord, -Y_coord)
				tmp14_sf <- st_sf(tmp14)
				
				# ----------------------
				#  Assign new solution name.
				nam4 <- paste("s", i, "sf", sep = "_")
				assign(nam4, tmp14_sf)
				cur_area_h <- as.integer(sum(st_area(tmp14_sf)) * 0.0001)
			
				# ----------------------
				#  Make solution an sp object.
				tmp5 <- as(tmp14_sf, "Spatial")
				nam5 <- paste("s", i, "sp", sep = "_")
				assign(nam5, tmp5)
	
				# ----------------------
				#  Buffer around solution to make a boundary of "local" area for next iteration of 
				#   optimization.
				bound <- st_buffer(st_union(tmp14_sf), 7500)
				bound_sp <- as(bound, "Spatial")
				
				# ----------------------
				#  Constrain input raster for the next iteration to the same general area as the boundary
				#   generated above.
				local_r <- raster::mask(moves_sm, bound_sp)
				nam7 <- paste("conn_in_upd", i, "01", sep = "_")
				tmp7 <- local_r
				assign(nam7, tmp7)
				
				# ----------------------
				#  Constrain planning units for the next iteration to the same general area as the features
				#   generated above.
				nam8 <- paste("pu_upd", i, sep = "_")
				tmp8a <- raster::crop(pu_in, bound_sp)
				tmp8b <-  st_as_sf(tmp8a) %>%
					mutate(size = as.integer(st_area(.)*0.0001)) %>%
					dplyr::filter(size > 490)
				tmp8 <- as(tmp8b, "Spatial")
				assign(nam8, tmp8)
	#}
	
	
				# ----------------------
				#  Check if the solution resulted in a full corridor from one planning unit to another.
				#   Using start PA from first iteration of this scenario.
				start_pa_id[[i]] <-start_pa_id[[i-1]]
								
				end_pu_sp <- tmp5
				end_pa_tmp <- raster::extract(pa_r, end_pu_sp) 
				end_pa <- unlist(lapply(end_pa_tmp, function(x) if (!is.null(x)) mean(x, na.rm=TRUE) else NA ))
				end_pa_id[[i]] <- end_pa
								
				diff_pa_tst <- ifelse(isEmpty(setdiff(end_pa_id[[i]], start_pa_id[[i]])), 0, 
				setdiff(end_pa_id[[i]], start_pa_id[[i]])) 
				diff_pa_tst[is.nan(diff_pa_tst)] <- 0
				
				if(any(diff_pa_tst > 0) == TRUE| nrow(tmp4) == nrow(tmp1)) {break}
				}
				
			}
			
		j = 1
		n_sol_steps[j] <- i		
		
		last_sol <- get(paste("s", i, "sf", sep = "_"))
		last_sol_sp <- as(last_sol, "Spatial")
		cur_area_h <- as.integer(sum(st_area(last_sol)) * 0.0001)
		last_sol_df <- as.data.frame(last_sol) %>%
			dplyr::mutate(scenario = j)
		
		whole_solution_out[[j]] <- last_sol_df
	
	all_sol <- do.call("rbind", whole_solution_out)
	all_sol_sf <- st_as_sf(all_sol)
	tot_area_h <-  as.integer(sum(st_area(all_sol_sf))*0.0001)
	
	#while(tot_area_h < 50000){
	#j = j + 1
	
	for(j in 2:20){
		
		# ----------------------
		#  Re-initiate loop output holders
		start_pa_id_loop <- list()
		end_pa_id_loop <- list()
		
		# ----------------------
		#  Specify that planning units selected in initiation problem solution that should become 
		#  "locked out" and not available for selection.
		loop_lock_out_tmp <- all_sol_sf %>% 
			mutate(locked_out = 1)
		locked_out <- as(loop_lock_out_tmp, "Spatial")
			
		# ----------------------
		#  Set up looped initial problem.
		p_init_loop <- problem(x = pu_start_in, features = moves_sm, cost_column = "area_h") %>%
			add_max_utility_objective(505) %>%
			add_locked_out_constraints(locked_out) %>%
			add_binary_decisions() %>%
			add_gurobi_solver()
		
			# ----------------------
		#  Solve initial problem.
		s_1_loop <- solve(p_init_loop)
		s_1_loop_sf <- st_as_sf(s_1_loop) %>%
			dplyr::filter(solution_1 == 1) %>%
			mutate(locked_in = 1) %>%
			mutate(select_order = 1)
		cur_area_h <- as.integer(sum(st_area(s_1_loop_sf)) * 0.0001)
			
		# ----------------------
		#  Make solution an sp object and then a raster.
		s_1_loop_sp <- as(s_1_loop_sf, "Spatial")
		
		# ----------------------
		#  Buffer around solution to make a boundary of "local" area for next iteration of 
		#   optimization.
		bound <- st_buffer(st_union(s_1_loop_sf), 5000)
		bound_sp <- as(bound, "Spatial")
		
		# ----------------------
		#  Constrain input raster for the next iteration to the same general area as the boundary
		#   generated above.
		local_r_loop <- raster::mask(moves_sm, bound_sp)
		conn_in_upd_loop_1_01 <- local_r_loop
		
		# ----------------------
		#  Constrain planning units for the next iteration to the same general area as the features
		#   generated above.
		pu_upd_1_loop_tmp <- raster::crop(pu_in, bound_sp)
		pu_upd_1_loop_sf <- st_as_sf(pu_upd_1_loop_tmp) %>%
			mutate(size = as.integer(st_area(.)*0.0001)) %>%
			dplyr::filter(size > 490)
		pu_upd_1_loop <- as(pu_upd_1_loop_sf, "Spatial")
		
		# ----------------------
		#  Check if the solution resulted in a full corridor from one planning unit to another.
		start_pu_sp <- as(s_1_loop_sf, "Spatial")
		start_pa_tmp <- raster::extract(pa_r, start_pu_sp) 
		start_pa <- unlist(lapply(start_pa_tmp, function(x) if (!is.null(x)) mean(x, na.rm=TRUE) else NA ))
		start_pa_id_loop[[1]] <- start_pa
		
		end_pu_sp <- as(s_1_loop_sf, "Spatial")
		end_pa_tmp <- raster::extract(pa_r, end_pu_sp) 
		end_pa <- unlist(lapply(end_pa_tmp, function(x) if (!is.null(x)) mean(x, na.rm=TRUE) else NA ))
		end_pa_id_loop[[1]] <- end_pa
		
		diff_pa_tst_loop <- ifelse(isEmpty(setdiff(end_pa_id_loop[[1]], start_pa_id_loop[[1]])), 0, 
			setdiff(end_pa_id_loop[[1]], start_pa_id_loop[[1]])) 
		diff_pa_tst_loop[is.nan(diff_pa_tst_loop)] <- 0
			
		if(any(diff_pa_tst_loop > 0) == FALSE){			
		
			for(k in 2:niter_max_corr){

				# ----------------------
				#  Specify that planning units selected in initiation problem solution that should become 
				#  "locked out" and not available for selection.
				locked_out_loop_tmp <- all_sol_sf %>% 
					mutate(locked_out = 1)
				locked_out <- as(locked_out_loop_tmp, "Spatial")
		
				# ----------------------
				#  Specify that planning unit selected in initiation or previous iteration is locked in.
				nam1_loop <- paste("s", k-1, "loop", "sf", sep = "_")
				tmp1_loop <- get(nam1_loop)
				locked_in <- as(tmp1_loop, "Spatial")
			
				# ----------------------
				#  Designate planning unit and connectivity features inputs.
				nam_upd_feat <- paste("conn_in_upd_loop", k-1, "01", sep = "_")
				upd_feat_tmp <- get(nam_upd_feat)
				
				nam_upd_pu <- paste("pu_upd", k-1, "loop", sep = "_")
				upd_pu_tmp2 <- get(nam_upd_pu)
				
				# ----------------------
				#  Remove the locked out planning units (those already selected in previous loop runs)
				#   from the planning unit set.
				if(nrow(inter <- st_intersection(all_sol_sf, st_as_sf(upd_pu_tmp2))) < 1){
					upd_pu_tmp <- upd_pu_tmp2} else {
					inter <- st_intersection(all_sol_sf, st_as_sf(upd_pu_tmp2))
					upd_pu_tmp_sf <- st_erase(st_as_sf(upd_pu_tmp2), st_as_sf(locked_out))
					upd_pu_tmp <- as(upd_pu_tmp_sf, "Spatial")
					}
					
				# ----------------------
				#  Set up problem to add planning units to first planning unit selected.
				p_loop <- problem(x = upd_pu_tmp, features = upd_feat_tmp, cost_column = "area_h") %>%
					add_max_utility_objective(505*k) %>%
					#add_boundary_penalties(0.0000002, 0.5) %>%
					add_locked_in_constraints(locked_in) %>%
					add_binary_decisions() %>%
					add_gurobi_solver()
		
				# ----------------------
				#  Solve iteration problem.
				tmp3_loop <- solve(p_loop)
				tmp4_loop <- st_as_sf(tmp3_loop) %>%
					dplyr::filter(solution_1 == 1) %>%
					mutate(locked_in = 1) %>%
					mutate(select_order = k)
				
				# ----------------------
				#  Combine with previous solution and remove duplicates ---
				#   used to obtain correct selection order number.
				tmp1_loop_df <- as.data.frame(tmp1_loop)
				tmp1_loop_coords <- as.data.frame(st_coordinates(st_cast(tmp1_loop, "MULTIPOLYGON"))) %>%
					#st_geometry_type(tmp1_loop)[1]))) %>%
					group_by_at(ncol(.)) %>%
					slice(1) %>%
					as.data.frame(.) %>%
					mutate(X_coord = round(X), Y_coord = round(Y)) 
				tmp1_loop_dup <- cbind(tmp1_loop_df, tmp1_loop_coords) %>%
					dplyr::select(area_h, solution_1, locked_in, select_order,
					geometry, X_coord, Y_coord) 
				
				tmp4_loop_df <- as.data.frame(tmp4_loop)
				tmp4_loop_coords <- as.data.frame(st_coordinates(st_cast(tmp4_loop, "MULTIPOLYGON"))) %>%
					group_by_at(ncol(.)) %>%
					slice(1) %>%
					as.data.frame(.) %>%
					mutate(X_coord = round(X), Y_coord = round(Y))  
				tmp4_loop_dup <- cbind(tmp4_loop_df, tmp4_loop_coords) %>%
					dplyr::select(area_h, solution_1, locked_in, select_order,
					geometry, X_coord, Y_coord) 
				
				tmp14_loop <- rbind(tmp1_loop_dup, tmp4_loop_dup) %>%
					arrange(X_coord, desc(Y_coord), select_order) %>%
					group_by(X_coord, Y_coord) %>%
					slice(1) %>%
					as.data.frame() %>%
					select(-X_coord, -Y_coord)
				tmp14_loop_sf <- st_sf(tmp14_loop)
				
				# ----------------------
				#  Assign new solution name.
				nam4_loop <- paste("s", k, "loop_sf", sep = "_")
				assign(nam4_loop, tmp14_loop_sf)
				cur_area_h <- as.integer(sum(st_area(tmp14_loop_sf)) * 0.0001)
			
				# ----------------------
				#  Make solution an sp object.
				tmp5_loop <- as(tmp14_loop_sf, "Spatial")
				nam5_loop <- paste("s", k, "loop_sp", sep = "_")
				assign(nam5_loop, tmp5_loop)
	
				# ----------------------
				#  Buffer around solution to make a boundary of "local" area for next iteration of 
				#   optimization.
				bound <- st_buffer(st_union(tmp14_loop_sf), 5000)
				bound_sp <- as(bound, "Spatial")
				
				# ----------------------
				#  Constrain input raster for the next iteration to the same general area as the boundary
				#   generated above.
				local_r_loop <- raster::mask(moves_sm, bound_sp)
				nam7_loop <- paste("conn_in_upd_loop", k, "01", sep = "_")
				tmp7_loop <- local_r_loop
				assign(nam7_loop, tmp7_loop)
				
				# ----------------------
				#  Constrain planning units for the next iteration to the same general area as the features
				#   generated above.
				nam8_loop <- paste("pu_upd", k, "loop", sep = "_")
				tmp8a_loop <- raster::crop(pu_in, bound_sp)
				tmp8b_loop <-  st_as_sf(tmp8a_loop) %>%
					mutate(size = as.integer(st_area(.)*0.0001)) %>%
					dplyr::filter(size > 490)
				tmp8_loop <- as(tmp8b_loop, "Spatial")
				assign(nam8_loop, tmp8_loop)

				# ----------------------
				#  Check if the solution resulted in a full corridor from one planning unit to another.
				#   Using start PA from first iteration of this scenario.
				start_pa_id_loop[[k]] <-start_pa_id_loop[[k-1]]
			
				end_pu_sp <- tmp5_loop
				end_pa_tmp <- raster::extract(pa_r, end_pu_sp) 
				end_pa <- unlist(lapply(end_pa_tmp, function(x) if (!is.null(x)) mean(x, na.rm=TRUE) else NA ))
				end_pa_id_loop[[k]] <- end_pa
				
				diff_pa_tst_loop <- ifelse(isEmpty(setdiff(end_pa_id_loop[[k]], start_pa_id_loop[[k]])), 0, 
					setdiff(end_pa_id_loop[[k]], start_pa_id_loop[[k]])) 
				diff_pa_tst_loop[is.nan(diff_pa_tst_loop)] <- 0

				if(any(diff_pa_tst_loop > 0) == TRUE | nrow(tmp4_loop) == nrow(tmp1_loop)) {break}
				}
			
		n_sol_steps[j] <- k			
		last_sol <- get(paste("s", k, "loop_sf", sep = "_"))
		last_sol_df <- as.data.frame(last_sol) %>%
			dplyr::mutate(scenario = j)
		whole_solution_out[[j]] <- last_sol_df
	
		all_sol <- do.call("rbind", whole_solution_out)
		all_sol_sf <- st_as_sf(all_sol)
		tot_area_h <-  as.integer(sum(st_area(all_sol_sf))*0.0001)
		}
		}
		
		

	
	

	
	# ----------------------
	#  Save plot to see planning units selected over time.
	 solution_p <- ggplot() +
		geom_sf(data = main_sabah_sf, colour = "grey50", fill = "transparent", alpha = 0.7) +
		#geom_sf(data = border_sarawak_sf, colour = "grey50", fill = "grey80") +
		#geom_sf(data = border_kali_sf, colour = "grey50", fill = "grey80") +
		#geom_sf(data = acd_agg_sf_for, colour = "#1B792F", fill = "#1B792F", alpha = 0.7) +
		#geom_sf(data = sa_grid_sf, fill = "transparent", colour = "grey20") +
		geom_sf(data = all_sol_sf,  aes(fill = factor(scenario))) +
		#geom_sf(data = s1_sf) +
		#geom_sf(data = corr_sw,  aes(fill = factor(scenario))) +
			#scale_fill_brewer(type = "seq", palette = "Blues", direction = -1) +
		#geom_sf(data = all_sol_sf_clust,  aes(fill = factor(group))) +
		#geom_sf(data = all_sol_sf_arr,  aes(fill = tot_order)) +
		#scale_fill_distiller(type = "seq", palette = "Reds", direction = 1) +	
		geom_sf(data = pa_sf, fill = "darkorange3", colour = "transparent") +
		coord_sf(crs = st_crs(32650)) +
		ylab("Latitude") +
		xlab("Longitude") +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		theme_bw()
	solution_p

	
	
	save(whole_solution_out, file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/whole_solution_out.Rdata")
	save(all_sol_sf, file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/all_sol_sf.Rdata")

	
	
	
	
	
	
	
	
	load(file = "conn_move_sim/conn_optim_out/Final/all_sol_sf.Rdata")
	 
	 
	# ----------------------
	#  Adjust solution so that only complete corridors are selected, not same park
	#   connections. And add variable for converting to raster.
	corr_scen <- all_sol_sf %>%	
		split(.$scenario) %>%
		lapply(st_union) %>%
		do.call(c,  .) %>%
		st_cast("POLYGON")
	corr_scen_sp <- as(corr_scen, "Spatial")
	corr_scen_sf <- st_as_sf(corr_scen_sp) %>% 
		mutate(size = as.integer(st_area(.))*0.0001) %>%
		dplyr::filter(size > 20000) %>%
		mutate(scenario = 1:n()) %>%
		dplyr::filter(scenario < 20) %>%
		dplyr::filter(scenario != 1, scenario != 5, scenario != 6) %>%
		mutate(scenario_new = 1:n()) %>%
		mutate(scenario_r = (n()+1) + desc(scenario_new))
		
	# ----------------------
	#  Plot selected corridors.
	corr_p <- ggplot() +
		geom_sf(data = border_sabah_sf, colour = "grey50", fill = "grey50", alpha = 0.7) +
		geom_sf(data = border_sarawak_sf, colour = "grey50", fill = "grey80") +
		geom_sf(data = border_kali_sf, colour = "grey50", fill = "grey80") +
		geom_sf(data = acd_agg_sf_for, colour = "#1B792F", fill = "#1B792F", alpha = 0.3) +
		#geom_sf(data = sa_grid_sf, colour = "transparent", fill = "grey90") +
		geom_sf(data = corr_scen_sf,  aes(fill = factor(scenario_new))) +
		scale_fill_brewer(type = "seq", palette = "Reds", direction = -1) +
		geom_sf(data = pa_in_b, fill = "darkorange3", colour = "transparent") +
		coord_sf(crs = st_crs(32650)) +
		xlab("Latitude") +
		ylab("Longitude") +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		theme_bw()
	corr_p
	
	# ----------------------
	#  Covert to raster for feature input.
	tmp2 <- smooth(corr_scen_sf, method = "ksmooth", smoothness = 3)
	area_thresh <- units::set_units(100, ha)
	tmp3 <- fill_holes(tmp2, area_thresh)	%>%
		st_buffer(2000)
	tmp3_sp <- as(tmp3, "Spatial")
	corr_r <- rasterize(tmp3_sp, r_template, field = tmp3_sp$scenario)
	x_min <- corr_r@data@min
	x_max <- corr_r@data@max
	corr_out_r <- raster::calc(corr_r, range01)
	
	corr_r <- rasterize(corr_scen_sf, r_template, field = corr_scen_sf$scenario)
	corr_r <- rasterize(corr_scen_sf, r_template, field = corr_scen_sf$scenario)
	x_min <- corr_r@data@min
	x_max <- corr_r@data@max
	corr_out_r <- raster::calc(corr_r, range01)
	
	
	corr_r[corr_r == 6] <- 0
	corr_r[corr_r == 5] <- 0
	#corr_sm <- focal(corr_r,  w = matrix(1, 21, 21))	
	corr_sm <- corr_r
	#corr_sm[corr_sm == 0] <- NA
	x_min <- corr_sm@data@min
	x_max <- corr_sm@data@max
	corr_out_r <- raster::calc(corr_sm, range01)
	
	#writeRaster(corr_out_r, "C:/Users/saraw/Desktop/tmp/corr_out_r_small.grd")
	
	# ----------------------	
	#  Plot raster.
	corr_r_p <- levelplot(corr_out, par.settings = r_theme, 
		#main = "Probability of use during species movement \n",
		xlab= "Longitude (UTM)",
		ylab="Latitude (UTM)",
		margin = FALSE) +
		layer(sp.polygons(main_ssk_sp, lwd = 0.8, col = 'grey40')) +
		#layer(sp.polygons(for_sp, lwd = 0.8, col = "darkgreen")) + #, fill = "darkgreen", alpha = 0.2)) +
		layer(sp.polygons(pa_clust_sp, lwd = 0.8, col = "darkorange", fill = "darkorange"))
	corr_r_p

	

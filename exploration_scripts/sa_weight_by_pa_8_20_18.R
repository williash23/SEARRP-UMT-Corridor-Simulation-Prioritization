###############################################################################
#  Connectivity weighting based on PA area and distance
#  January 23, 2018; last updated February 28, 2018
#  Script to generate scaled weight for each study area grid cell that is based on the area and 
#   distance from the 2 PAs it is closest to. 
#  Sara Williams
###############################################################################

#install.packages(c("reshape", "tidyr"), lib = "C:/Users/saraw/Documents/R/win-library/3.5")

library(tidyr, lib.loc = "C:/Users/saraw/Documents/R/win-library/3.5")
library(reshape, lib.loc = "C:/Users/saraw/Documents/R/win-library/3.5")

# =============================================================================
#  Load packages.
# =============================================================================
library(sf)
library(sp)
library(raster)
library(rgdal)
library(rgeos)
library(reshape)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(smoothr)




# =============================================================================
#  Load data.
# =============================================================================		
	
	# ----------------------
	# Border of around Sabah.
	load("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/border_sabah_sf.Rdata")
	border_sabah_sp <- as(border_sabah_sf, "Spatial")

	# ----------------------
	#  Load connectivity information from animal movement scenarios for template.
	moves_in <- raster("C:/Users/saraw/Documents/SEARRP_Analyses/conn_move_sim/moves_in/1000sim_3000steps_init_loc_wt_pa_size_over5000h_turn_conc80_lg_area/moves_in.grd")
	
	# ----------------------
	#  Load study area grid that holds planning units for connectivity assessment.
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/sa_grid_sf.Rdata")

	# ----------------------
	#  "Locked in" already existing protected areas.
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/ssk_pa_near_sf_clust_1k.Rdata")
	
	
	
# =============================================================================
#  Specify parameters for weighting.
# =============================================================================
	
	# ----------------------
	#  Weighting factors 
	#   Distance: 1
	#   Size of the smaller PA: 2
	#   Size of the larger PA: 5
	dist_wt <- 1
	smaller_pa_wt <- 2
	larger_pa_wt <- 5
	
	
# =============================================================================
#  Make sure PAs have sequential IDs
# =============================================================================	

	# ----------------------
	#  Generate new sequential id for PAs sf
	pas_sf_tmp <- ssk_pa_near_sf_clust_1k %>%
		st_intersection(st_set_crs(st_as_sf(as(raster::extent(231556, 751756, 400000, 783428), 
		"SpatialPolygons")), st_crs(ssk_pa_near_sf_clust_1k))) %>%
		mutate(size = as.integer(st_area(.)*0.0001)) %>%
		dplyr::filter(size > 5000) %>%
		dplyr::select(-group)
	
	area_thresh <- units::set_units(100, km^2)	
	pas_sf <- pas_sf_tmp %>%
		drop_crumbs(area_thresh) %>%
		mutate(new_id = row_number()) %>%
		st_buffer(10)

	pas_sp  <- as(pas_sf, "Spatial")
	
	
	
# =============================================================================
#  Set up planning unit later and obtain centroids of each planning unit.
# =============================================================================		
	
	# ----------------------
	#  If using planning units: input which planning unit size is desired.
	pu_sf <- sa_grid_sf

	# ----------------------
	#  Convert to sp object.
	pu <- as(pu_sf, "Spatial")	
	
	# ----------------------
	#  Get centroids
	# pts <- as.data.frame(coordinates(r_template))
	# coordinates(pts) <- ~x+y
	# #pts <- gCentroid(pu, byid = TRUE)
	# pts_sf <- st_as_sf(pts)



# =============================================================================
#  Get distance (in meters) and ID of the closest PA to each planning unit.
# =============================================================================	

	# ----------------------
	#  Get the distances (m) between each centroid and each PA.
	pu_to_pa_dist <- as.data.frame(gDistance(pas_sp, pu,  byid = TRUE))
	
	# ----------------------
	#  Get the ID of the closest PA to each planning unit centroid.
	closest_pa_to_pu_id <- integer(length(pu))
	for(i in 1:length(pu)){
		pu_i_sf <- pu_sf[i, 1]
		pu_i <- as(pu_i_sf, "Spatial")
		closest_pa_to_pu_id[i] <- pas_sp$new_id[which.min(gDistance(pas_sp, pu_i, byid = TRUE))]
		}

	# ----------------------
	#  Get the distance (m) of the closest PA to each planning unit centroid.
	closest_pa_to_pu_dist <- as.data.frame(apply(gDistance(pu, pas_sp,  byid = TRUE), 2, min)) 
	colnames(closest_pa_to_pu_dist ) <- "dist_to_1"
	closest_pa_to_pu_dist_v <- as.vector(closest_pa_to_pu_dist[,1])


	
# =============================================================================
#  Get distance (in meters) and ID of the SECOND closest PA to each planning unit.
# =============================================================================	
	
	# ----------------------
	#  Temporarily increase the distance of the closest PA to each planning unit (determined above) 
	#   in order to select second closest PA. Increase to a distance larger that any distances between
	#   planning units and PAs.
	fake_dist <- max(pu_to_pa_dist) + 10000

	second_pa_to_pu_dist_m <- matrix(nrow = nrow(pu_to_pa_dist), ncol = ncol(pu_to_pa_dist))
	for(i in 1:nrow(pu_to_pa_dist)){
		tmp1 <- as.numeric(pu_to_pa_dist[i,])
		tmp2 <- closest_pa_to_pu_dist_v[i]
		find_same <- function(x){which(x == tmp2)}
		tmp3 <- find_same(tmp1)
		tmp4 <- tmp1
		tmp4[tmp3] <- fake_dist
		second_pa_to_pu_dist_m[i,] <- tmp4
		}
	
	# ----------------------
	#  Get the distance (m) of the second closest PA to each planning unit centroid.
	second_pa_to_pu_dist_tmp <- matrix(nrow = nrow(pu_to_pa_dist), ncol = 1)
	for(i in 1:nrow(second_pa_to_pu_dist_m)){
		tmp1 <- as.numeric(second_pa_to_pu_dist_m[i,])
		tmp2 <- min(tmp1)
		second_pa_to_pu_dist_tmp[i,1] <- tmp2
		}
		
	# ----------------------
	#  Convert to data frame
	second_pa_to_pu_dist <- as.data.frame(second_pa_to_pu_dist_tmp)
	colnames(second_pa_to_pu_dist) <- "dist_to_2"
	
	# ----------------------
	#  Get the ID of the second closest PA to each planning unit centroid.
	second_pa_to_pu_id <- integer(length(pu))
	for(i in 1:length(pu)){
		tmp1 <- as.vector(pu_to_pa_dist[i,])
		tmp2 <- second_pa_to_pu_dist[i,]
		tmp3 <- which(tmp1 == tmp2)[1]
		second_pa_to_pu_id[i] <- tmp3
		}

		

# =============================================================================
#  Join info for closest and second closest PA to each planning unit to the planning unit data frame.
# =============================================================================	
	
	# ----------------------
	#  Bind ID's and distances for closest PA and second closest PA.
	closest_pa_to_pu_df <- as.data.frame(cbind(closest_pa_to_pu_id, closest_pa_to_pu_dist))
	second_pa_to_pu_df <- as.data.frame(cbind(second_pa_to_pu_id, second_pa_to_pu_dist))
	
	# ----------------------
	#  Bind each data frame generated above to planning unit sf object.
	pu_sf_closest_pa <- bind_cols(pu_sf, closest_pa_to_pu_df) 
	pu_sf_second_pa <- bind_cols(pu_sf, second_pa_to_pu_df) 
	
	# ----------------------
	#  Select only the area, ID, and name columns from the PA sf object and scale the area of all PAs.
	pas_df <- pas_sf# %>%
		#mutate(pa_area_h = as.numeric(st_area(.) / 10000))
	st_geometry(pas_df) <- NULL

	# ----------------------
	#  Join each planing unit sf object, which now has the closest and second closest PA info attached
	#   to it, to the PA sf objects by PA ID.
	pu_sf_1 <- left_join(pu_sf_closest_pa, pas_df, 
			by = c("closest_pa_to_pu_id" = "new_id")) %>%
		dplyr::rename(pa_1_area_h = size)
	pu_sf_2 <- left_join(pu_sf_second_pa, pas_df, 
			by = c("second_pa_to_pu_id" = "new_id")) %>%
		dplyr::rename(pa_2_area_h = size) 



# =============================================================================
#  Create a weighted value for the importance of each planning unit for connectivity.
# =============================================================================

	# ----------------------
	#  Create new column with calculate weight value.
	pu_wt_sf_tmp <- bind_cols(pu_sf_1, pu_sf_2) %>%
		mutate(dist_bt_near_pas_km = (dist_to_1 + dist_to_2)/1000) %>%
		mutate(inv_dist_bt_near_pas_log = 1/log10(dist_bt_near_pas_km)) %>%
		rowwise(.) %>%
		mutate(larger_pa_h = max(pa_1_area_h,  pa_2_area_h)) %>%
		mutate(larger_pa_log = log10(larger_pa_h)) %>%
		mutate(smaller_pa_h = min(pa_1_area_h,  pa_2_area_h)) %>%
		mutate(smaller_pa_log = log10(smaller_pa_h)) %>%
		as.data.frame(.) %>%
		dplyr::select(-geometry1)
	
	# ----------------------
	#  Recast as sf object.
	pu_wt_sf_tmp2 <- pu_wt_sf_tmp %>%
		st_as_sf(sf_column_name = "geometry")
			
	# ----------------------
	#  Multiply be determined weights.
	pu_wt_sf_tmp3 <- pu_wt_sf_tmp2 %>%
		mutate(pu_wt =  (dist_wt*inv_dist_bt_near_pas_log) + 
			(smaller_pa_wt*smaller_pa_log) + 
			(larger_pa_wt*larger_pa_log))
	
	#ssk_pa_buff <- st_buffer(ssk_pa_near_sf_clust_1k, 300)
	
	#pu_wt_sf <- pu_wt_sf_tmp3 %>%
		#st_difference(ssk_pa_buff)
	
	# ----------------------	
	#  Make a smaller data frame for saving and plotting and convert to sp object.
	#pu_wt_sf_small <- pu_wt_sf %>%
	pu_wt_sf_small <- pu_wt_sf_tmp3 %>%
		dplyr::select(pu_wt)
	pu_wt_sp <- as(pu_wt_sf_small, "Spatial")
	pu_wt_sf_plot <- st_as_sf(pu_wt_sp)
	
	# ----------------------	
	#  Save sf and sp objects.
	# save(pu_wt_sf, file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/pu_wt_sf.Rdata")
	# save(pu_wt_sf_plot, file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/pu_wt_sf_plot.Rdata")
	# save(pu_wt_sp, file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/pu_wt_sp.Rdata")
	
	
	
# =============================================================================
#  Rasterize
# =============================================================================
	
	# ----------------------
	#  Load input from above if needed.
	#load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/pu_wt_sp.Rdata")
	
	# ----------------------
	#  Create raster following template of study area (cell values are empty)
	r_mat <- matrix(NA, nrow(moves_in), ncol(moves_in))
	r_template <- raster(r_mat)
	extent(r_template) <- extent(moves_in)
	projection(r_template) <- CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0") 

	# ----------------------
	#  Reformat as raster
	pt_wt_r <- raster::rasterize(pu_wt_sp, r_template, field = pu_wt_sp$pu_wt, fun = modal, na.rm = TRUE)

				# ----------------------	
				#  Smooth using a focal window
				wt_smooth <- focal(pt_wt_r,  w = matrix(1, 33, 33), na.rm = TRUE, pad = TRUE)
		
		
	# ----------------------
	#  Rescale to 0 to 1 scale, if wanted for use as a conservation feature input.
	x_min <- pt_wt_r@data@min
	x_max <- pt_wt_r@data@max
	range01 <- function(x){(x-x_min)/(x_max - x_min)}
	pt_wt_r_01 <-raster::calc(pt_wt_r, range01)
	
	

	
	pa_r <- rasterize(pas_sp, r_template)
	pa_r[pa_r > 1] <- 1
	pa_r[is.na(pa_r)] <- 0
	pa_rs <- raster::resample(pa_r, pt_wt_r_01)
	pa_rs[pa_rs > 0] <- 1
	stk <- raster::stack(pt_wt_r_01, pa_rs)
	r <- sum(stk, na.rm = TRUE)
	r[r > 1] <- 1
	r[r == 0] <- 0.01
	
	tr <- transition(r, transitionFunction = mean, directions = 4)
	tr1 <- geoCorrection(tr, type = "c", multpl = FALSE)
	pas_tmp <- raster::crop(pas_sp, r)
	pa_cents <- gCentroid(pas_tmp, byid = TRUE)
	

	knn <- knearneigh(pa_cents, k=2)$nn
 	
	paths1_ls <- list()
	for(i in 1:15){
			tmp <- shortestPath(tr1, pa_cents[i], pa_cents[knn[i,1]], output="SpatialLines")
			paths1_ls[[i]] <- tmp
			}
			
	paths1 <- do.call(rbind, paths1_ls)
	paths1_sf <- st_as_sf(paths1)
	
	paths2_ls <- list()
	for(i in 1:15){
			tmp <- shortestPath(tr1, pa_cents[i], pa_cents[knn[i,2]], output="SpatialLines")
			paths2_ls[[i]] <- tmp
			}
			
	paths2 <- do.call(rbind, paths2_ls)
	paths2_sf <- st_as_sf(paths2)
	
	paths3_ls <- list()
	for(i in 1:15){
			tmp <- shortestPath(tr1, pa_cents[i], pa_cents[knn[i,3]], output="SpatialLines")
			paths3_ls[[i]] <- tmp
			}
			
	paths3 <- do.call(rbind, paths3_ls)
	paths3_sf <- st_as_sf(paths3)
	
	paths_sf <- rbind(paths1_sf, paths2_sf, paths3_sf)
	paths_sp <- as(paths_sf, "Spatial")
		
		
		
	# ----------------------
	#  Rescale to 0 to 1 scale, if wanted for use as a conservation feature input.
	x_min <- out@data@min
	x_max <- out@data@max
	range01 <- function(x){(x-x_min)/(x_max - x_min)}
	wt <-raster::calc(out, range01)
	
	
	
	# ----------------------
	#  Save rasters.
	# writeRaster(pt_wt_r, "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/pt_wt_r.grd")
	# writeRaster(pt_wt_r_01, "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/pt_wt_r_01.grd")
	
	

# =============================================================================
#  Plot study area with each planning unit weighted by connectivity value 
#   of distance and area of PAs.
# =============================================================================

	# ----------------------
	#  Reload weighted PA data if needed.
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/pu_wt_sf.Rdata")
	
	# ----------------------
	#  All protected areas in Sabah, Sarawak, and Kalimantan.
	load("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/sabah_pa_sf.Rdata")
	load("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/sarawak_pa_sf.Rdata")
	load("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/kali_pa_sf.Rdata")
	
	# ----------------------
	#  Region/country borders.
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/border_sabah_sf.Rdata")
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/border_sarawak_sf.Rdata")
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/border_kali_sf.Rdata")
	
	border_sabah_sp <- as(border_sabah_sf, "Spatial")
	border_sarawak_sp <- as(border_sarawak_sf, "Spatial")
	border_kali_sp <- as(border_kali_sf, "Spatial")
	
	sabah_pa_sp <- as(sabah_pa_sf, "Spatial")
	sarawak_pa_sp <- as(sarawak_pa_sf, "Spatial")
	kali_pa_sp <- as(kali_pa_sf, "Spatial")
	
	
	# ----------------------	
	#  Prep theme for plotting.
	r_theme <- rasterTheme(brewer.pal(10, "Reds"))
	
	# ----------------------	
	#  Plot as raster.
	wt_p <- levelplot(pt_wt_r_01, par.settings = r_theme, 
		xlab= "Longitude (UTM)",
		ylab="Latitude (UTM)",
		margin = FALSE) +
		layer(sp.polygons(border_sabah_sp, lwd = 0.8, col = 'grey40')) + 
		layer(sp.polygons(border_sarawak_sp, lwd = 0.8, col = 'grey60')) + 
		layer(sp.polygons(border_kali_sp, lwd = 0.8, col = 'grey60')) + 
		layer(sp.polygons(sabah_pa_sp, lwd = 0.8, col = 'goldenrod1', fill = 'goldenrod1')) +
		layer(sp.polygons(sarawak_pa_sp, lwd = 0.8, col = 'lightgoldenrod1', fill = 'lightgoldenrod1')) +
		layer(sp.polygons(kali_pa_sp, lwd = 0.8, col = 'lightgoldenrod1', fill = 'lightgoldenrod1')) 
	wt_p
	
	# ----------------------
	#  Plot as ggplot2 sf objects. Fortify PA poygons for plotting ID number on map.
	pas_sp_df <- fortify(pas_sp) %>%
		group_by(id) %>%
		slice(1) %>%
		as.data.frame() %>%
		dplyr::select(long, lat, id)
	
	# ----------------------
	#  Scale weight planning units from 0 to 1.
	x_min <-min(pu_wt_sf$pu_wt)
	x_max <-max(pu_wt_sf$pu_wt)
	pu_wt_sf_plot_sc <- pu_wt_sf %>%
		mutate(wt_sc = (pu_wt - x_min)/(x_max - x_min))
		
	# ----------------------
	#  Plot.
	pu_wt_p <- ggplot() +
		geom_sf(data = border_sabah_sf, colour = "grey50", fill = "grey50", alpha = 0.7) +
		geom_sf(data = border_sarawak_sf, colour = "grey50", fill = "grey80") +
		geom_sf(data = border_kali_sf, colour = "grey50", fill = "grey80") +
		#geom_sf(data = pu_wt_sf_plot, aes(colour = pu_wt), size = 2, alpha = 0.4) +
		geom_sf(data = pu_wt_sf_plot_sc, aes(colour = wt_sc), alpha = 0.4) +
		scale_colour_distiller(type = "seq", palette = "Reds", direction = 1) +		
		geom_sf(data = sabah_pa_sf, colour = "darkseagreen3", fill = "transparent") +
		geom_sf(data =sarawak_pa_sf, colour = "darkseagreen4", fill = "transparent") +
		geom_sf(data = kali_pa_sf, colour = "darkseagreen4", fill = "transparent") +
		#geom_text(data = pas_sp_df, aes(x = long, y = lat, label = id), size = 2) +
		coord_sf(crs = st_crs(32650)) +
		xlab("Latitude") +
		ylab("Longitude") +
		xlim(315000, 755000) +
		ylim(455000, 815000) +
		theme_bw()
	pu_wt_p
	

	
# =============================================================================	
###############################################################################




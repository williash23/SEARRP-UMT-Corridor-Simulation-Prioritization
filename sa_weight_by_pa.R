###############################################################################
#  Connectivity weighting based on PA area and distance
#  January 23, 2018; last updated February 28, 2018
#  Script to generate scaled weight for each study area grid cell that is based on the area and 
#   distance from the 2 PAs it is closest to. 
#  Sara Williams
###############################################################################



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

# =============================================================================
#  Load data.
# =============================================================================		
	
	# ----------------------
	# Border of around Sabah.
	load("C:/Users/saraw/Documents/SEARRP_Analyses/optimization/border_sabah_sf.Rdata")
	border_sabah_sp <- as(border_sabah_sf, "Spatial")

	# ----------------------
	#  Load connectivity information from animal movement scenarios for template.
	spp_moves <- raster("C:/Users/saraw/Documents/SEARRP_Analyses/move_sim/moves_in.grd")
	
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
	pas_sf <- ssk_pa_near_sf_clust_1k %>%
		mutate(new_id = row_number()) %>%
		dplyr::select(-group)
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
	pts <- gCentroid(pu, byid = TRUE)
	#pts_sf <- st_as_sf(pts)


	tmp_pts <- pts[1:5]
	tmp_pa <- pas_sp[1:10,]
	
# =============================================================================
#  Get distance (in meters) and ID of the closest PA to each planning unit.
# =============================================================================	

	# ----------------------
	#  Get the distances (m) between each centroid and each PA.
	all_pts_to_pa_dist <- as.data.frame(gDistance(pas_sp, pts,  byid = TRUE))
		tst_all_pts_to_pa_dist <- as.data.frame(gDistance(tmp_pa, tmp_pts,  byid = TRUE))

	# ----------------------
	#  Get the ID of the closest PA to each planning unit centroid.
	near_1_pa_to_pts_id <- integer(length(pts))
	for(i in 1:length(pts)){
		#near_1_pa_to_pts_id[i] <- pas_sp$new_id[which.min(gDistance(pas_sp, pts[i,], byid = TRUE))]
			t_near_1_pa_to_pts_id[i] <- tmp_pa$new_id[which.min(gDistance(tmp_pa, tmp_pts[i,], byid = TRUE))]
		}

	# ----------------------
	#  Get the distance (m) of the closest PA to each planning unit centroid.
	near_1_pa_to_pts_dist <- as.data.frame(apply(gDistance(pts, pas_sp,  byid = TRUE), 2, min)) 
	colnames(near_1_pa_to_pts_dist) <- "dist_to_1"
	near_1_pa_to_pts_dist_v <- as.vector(near_1_pa_to_pts_dist[,1])

	t_near_1_pa_to_pts_dist <- as.data.frame(apply(gDistance(tmp_pts,tmp_pa,  byid = TRUE), 2, min)) 
	colnames(t_near_1_pa_to_pts_dist) <- "dist_to_1"
	t_near_1_pa_to_pts_dist_v <- as.vector(t_near_1_pa_to_pts_dist[,1])

	
# =============================================================================
#  Get distance (in meters) and ID of the SECOND closest PA to each planning unit.
# =============================================================================	
	
	# ----------------------
	#  Temporarily increase the distance of the closest PA to each planning unit (determined above) 
	#   in order to select second closest PA. Increase to a distance larger that any distances between
	#   planning units and PAs.
	fake_dist <- max(all_pts_to_pa_dist) + 10000
	t_near_2_pa_to_pts_dist_m <- matrix(nrow = nrow(tst_all_pts_to_pa_dist), ncol = ncol(tst_all_pts_to_pa_dist))
	for(i in 1:nrow(tst_all_pts_to_pa_dist)){
		tmp1 <- as.numeric(tst_all_pts_to_pa_dist[i,])
		tmp2 <- t_near_1_pa_to_pts_dist_v[i]
		find_same <- function(x){which(x == tmp2)}
		tmp3 <- find_same(tmp1)
		tmp4 <- tmp1
		tmp4[tmp3] <- fake_dist
		t_near_2_pa_to_pts_dist_m[i,] <- tmp4
		}
	
	near_2_pa_to_pts_dist_m <- matrix(nrow = nrow(all_pts_to_pa_dist), ncol = ncol(all_pts_to_pa_dist))
	for(i in 1:nrow(all_pts_to_pa_dist)){
		tmp1 <- as.numeric(all_pts_to_pa_dist[i,])
		tmp2 <- near_1_pa_to_pts_dist_v[i]
		find_same <- function(x){which(x == tmp2)}
		tmp3 <- find_same(tmp1)
		tmp4 <- tmp1
		tmp4[tmp3] <- fake_dist
		near_2_pa_to_pts_dist_m[i,] <- tmp4
		}
	
	# ----------------------
	#  Get the distance (m) of the second closest PA to each planning unit centroid.

	t_near_2_pa_to_pts_dist_tmp <- matrix(nrow = nrow(tst_all_pts_to_pa_dist), ncol = 1)
	for(i in 1:nrow(t_near_2_pa_to_pts_dist_m)){
		tmp1 <- as.numeric(t_near_2_pa_to_pts_dist_m[i,])
		tmp2 <- min(tmp1)
		t_near_2_pa_to_pts_dist_tmp[i,1] <- tmp2
		}
		
		
	near_2_pa_to_pts_dist_tmp <- matrix(nrow = nrow(all_pts_to_pa_dist), ncol = 1)
	for(i in 1:nrow(near_2_pa_to_pts_dist_m)){
		tmp1 <- as.numeric(near_2_pa_to_pts_dist_m[i,])
		tmp2 <- min(tmp1)
		near_2_pa_to_pts_dist_tmp[i,1] <- tmp2
		}
		
	# ----------------------
	#  Convert to data frame
	near_2_pa_to_pts_dist <- as.data.frame(near_2_pa_to_pts_dist_tmp)
	colnames(near_2_pa_to_pts_dist) <- "dist_to_2"
	
	# ----------------------
	#  Get the ID of the second closest PA to each planning unit centroid.
	near_2_pa_to_pts_id <- integer(length(pts))
	for(i in 1:length(pts)){
		tmp1 <- as.vector(all_pts_to_pa_dist[i,])
		tmp2 <- near_2_pa_to_pts_dist[i,]
		tmp3 <- which(tmp1 == tmp2)[1]
		near_2_pa_to_pts_id[i] <- tmp3
		}

		

# =============================================================================
#  Join info for closest and second closest PA to each planning unit to the planning unit data frame.
# =============================================================================	
	
	# ----------------------
	#  Bind ID's and distances for closest PA and second closest PA.
	near_1_pa_to_pts_df <- as.data.frame(cbind(near_1_pa_to_pts_id, near_1_pa_to_pts_dist))
	near_2_pa_to_pts_df <- as.data.frame(cbind(near_2_pa_to_pts_id, near_2_pa_to_pts_dist))
	
	# ----------------------
	#  Bind each data frame generated above to planning unit sf object.
	pts_sf_near_1_pa <- bind_cols(pts_sf, near_1_pa_to_pts_df) 
	pts_sf_near_2_pa <- bind_cols(pts_sf, near_2_pa_to_pts_df) 
	
	# ----------------------
	#  Select only the area, ID, and name columns from the PA sf object and scale the area of all PAs.
	pas_df <- pas_sf %>%
		mutate(pa_area_tmp = st_area(.) / 10000) %>%
		separate(pa_area_tmp, sep = " ", c("pa_area_h"), drop = TRUE) 
	pas_df$pa_area_h <- as.numeric(pas_df$pa_area_h)
	st_geometry(pas_df) <- NULL

	# ----------------------
	#  Join each planing unit sf object, which now has the closest and second closest PA info attached
	#   to it, to the PA sf objects by PA ID.
	pts_sf_1 <- left_join(pts_sf_near_1_pa, pas_df, 
			by = c("near_1_pa_to_pts_id" = "new_id")) %>%
		dplyr::rename(pa_1_area_h = pa_area_h)
	pts_sf_2 <- left_join(pts_sf_near_2_pa, pas_df, 
			by = c("near_2_pa_to_pts_id" = "new_id")) %>%
		dplyr::rename(pa_2_area_h = pa_area_h) 



# =============================================================================
#  Create a weighted value for the importance of each planning unit for connectivity.
# =============================================================================

	# ----------------------
	#  Create new column with calculate weight value.
	pts_wt_sf_tmp <- bind_cols(pts_sf_1, pts_sf_2) %>%
		mutate(dist_bt_near_pas_m = dist_to_1 + dist_to_2) %>%
		mutate(inv_dist_bt_near_pas_log = 1/log(dist_bt_near_pas_m)) %>%
		rowwise(.) %>%
		mutate(larger_pa_h = max(pa_1_area_h,  pa_2_area_h)) %>%
		mutate(larger_pa_log = log(larger_pa_h)) %>%
		mutate(smaller_pa_h = min(pa_1_area_h,  pa_2_area_h)) %>%
		mutate(smaller_pa_log = log(smaller_pa_h)) %>%
		as.data.frame(.) %>%
		dplyr::select(-geometry1)
	
	# ----------------------
	#  Recast as sf object.
	pts_wt_sf_tmp2 <- pts_wt_sf_tmp %>%
		st_as_sf(sf_column_name = "geometry")
			
	# ----------------------
	#  Multiply be determined weights.
	pts_wt_sf_tmp3 <- pts_wt_sf_tmp2 %>%
		mutate(pu_wt =  (dist_wt*inv_dist_bt_near_pas_log) + 
			(smaller_pa_wt*smaller_pa_log) + 
			(larger_pa_wt*larger_pa_log))
	
	ssk_pa_buff <- st_buffer(ssk_pa_near_sf_clust, 300)
	
	pts_wt_sf <- pts_wt_sf_tmp3 %>%
		st_difference(ssk_pa_buff)
	
	# ----------------------	
	#  Make a smaller data frame for saving and plotting and convert to sp object.
	pts_wt_sf_small <- pts_wt_sf %>%
		dplyr::select(pu_wt)
	pts_wt_sp <- as(pts_wt_sf_small, "Spatial")
	pts_wt_sf_plot <- st_as_sf(pts_wt_sp)
	
	# ----------------------	
	#  Save sf and sp objects.
	# save(pts_wt_sf, file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/pts_wt_sf.Rdata")
	# save(pts_wt_sf_plot, file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/pts_wt_sf_plot.Rdata")
	# save(pts_wt_sp, file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/pts_wt_sp.Rdata")
	
	
	
# =============================================================================
#  Rasterize
# =============================================================================
	
	# ----------------------
	#  Load input from above if needed.
	#load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/pts_wt_sp.Rdata")
	
	# ----------------------
	#  Create raster following template of study area (cell values are empty)
	r_mat <- matrix(NA, nrow(moves_in), ncol(moves_in))
	r_template <- raster(r_mat)
	extent(r_template) <- extent(moves_in)
	projection(r_template) <- CRS("+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0") 

	# ----------------------
	#  Reformat as raster
	pt_wt_r <- raster::rasterize(pts_wt_sp, r_template, 
		field = pts_wt_sp$pu_wt, fun = modal)

				# ----------------------	
				#  Smooth using a focal window
				# wt_smooth <- focal(pt_wt_r,  w = matrix(1, 3, 3), mean, na.rm = TRUE)
		
		
	# ----------------------
	#  Rescale to 0 to 1 scale, if wanted for use as a conservation feature input.
	x_min <- pt_wt_r@data@min
	x_max <- pt_wt_r@data@max
	range01 <- function(x){(x-x_min)/(x_max - x_min)}
	pt_wt_r_01 <-raster::calc(pt_wt_r, range01)
	
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
	load(file = "C:/Users/saraw/Documents/SEARRP_Analyses/optimization/pts_wt_sf.Rdata")
	
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
	x_min <-min(pts_wt_sf$pu_wt)
	x_max <-max(pts_wt_sf$pu_wt)
	pts_wt_sf_plot_sc <- pts_wt_sf %>%
		mutate(wt_sc = (pu_wt - x_min)/(x_max - x_min))
		
	# ----------------------
	#  Plot.
	pu_wt_p <- ggplot() +
		geom_sf(data = border_sabah_sf, colour = "grey50", fill = "grey50", alpha = 0.7) +
		geom_sf(data = border_sarawak_sf, colour = "grey50", fill = "grey80") +
		geom_sf(data = border_kali_sf, colour = "grey50", fill = "grey80") +
		#geom_sf(data = pts_wt_sf_plot, aes(colour = pu_wt), size = 2, alpha = 0.4) +
		geom_sf(data = pts_wt_sf_plot_sc, aes(colour = wt_sc), alpha = 0.4) +
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




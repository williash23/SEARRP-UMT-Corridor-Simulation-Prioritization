###############################################################################
#  Connectivity through movement simulation for bookend scenarios.
#  August 25, 2018; udpated March 4, 2019
#  Script to prepare existing protected areas input for connectivity movement simulations.
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
library(dplyr)
library(smoothr)



	
# =============================================================================
#  Load base layer input data.
# =============================================================================

	
	# ----------------------
	#  Exiting protected areas, clustered (within 5m of each other become single PA)
	load(file = "C:/Users/saraw/Documents/Prioritization/study_area_boundaries_and_pas/ssk_pa_near_sf_clust.Rdata")
	load(file = "C:/Users/saraw/Documents/Prioritization/study_area_boundaries_and_pas/ssk_pa_sf.Rdata")
	load(file = "C:/Users/saraw/Documents/Prioritization/study_area_boundaries_and_pas/main_sabah_sf.Rdata")
		
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
	
	
	
# =============================================================================
#  Process.
# =============================================================================

	# ----------------------
	#  Select only PA cluster with area of 5000 ha.
	pa_tmp1 <- st_intersects(ssk_pa_near_sf_clust, 
		st_set_crs(st_as_sf(as(raster::extent(231556, 751756, 400000, 783428), 
		"SpatialPolygons")), st_crs(ssk_pa_near_sf_clust)))
	
	pa_tmp2 <- ssk_pa_sf %>% #ssk_pa_near_sf_clust %>% #pa_tmp1 %>%
		mutate(size = as.integer(st_area(.)*0.0001)) %>%
		dplyr::select(size)
	
		thresh <- units::set_units(1, km)
		pa_tmp3 <- clusterSF(pa_tmp2, thresh)
		
		
		area_thresh <- units::set_units(1000, km^2)
		pa_tmp4 <- pa_tmp3 %>%
			mutate(size = as.integer(st_area(.)*0.0001)) %>%
			dplyr::filter(size > 1000)
			
		md_tmp1 <- pa_tmp4 %>%
			dplyr::filter(group == 18) %>%
			st_buffer(1500) %>%
			fill_holes(threshold = area_thresh)  %>%
			st_buffer(-1500) 
		md_tmp2 <- pa_tmp3 %>%
			filter(group == 18) %>%
			st_union(md_tmp1) %>%
			dplyr::select(-group.1)
		
		pa_tmp5 <- pa_tmp4 %>%
			dplyr::filter(group != 18) %>%
			rbind(md_tmp2)
		
		close_bound <- st_buffer(main_sabah_sf, 5000)
		close_pas <- as.numeric(which(st_intersects(pa_tmp5, close_bound, sparse = FALSE) == TRUE))
			
		pa_corr_sim <- pa_tmp5[close_pas,]
		save(pa_corr_sim, file = "C:/Users/saraw/Documents/Prioritization/pa_corr_sim.Rdata")
	
	
	#save(pa_sp_move_sim, file = "C:/Users/saraw/Documents/SEARRP_Analyses/conn_move_sim/pa_sp_move_sim.Rdata")
	#save(pa_sf_move_sim, file = "C:/Users/saraw/Documents/SEARRP_Analyses/conn_move_sim/pa_sf_move_sim.Rdata")
	
	
	
	
	pa_sf_move_sim <- rbind(pa_tmp6, md2_fill) %>%
		dplyr::filter(size != 70754) %>%
		dplyr::filter(size != 7559) %>%
		dplyr::filter(size != 17693) 
	pa_sf_move_sim$size[pa_sf$size > 300000] <- max(pa_tmp4$size)
	pa_sp_move_sim <- as(pa_sf, "Spatial")
	
	#save(pa_sp_move_sim, file = "C:/Users/saraw/Documents/SEARRP_Analyses/conn_move_sim/pa_sp_move_sim.Rdata")
	#save(pa_sf_move_sim, file = "C:/Users/saraw/Documents/SEARRP_Analyses/conn_move_sim/pa_sf_move_sim.Rdata")
	
	pa_sf_smooth <- pa_sf %>%
		st_buffer(1000) %>%
		smooth(method = "chaikin") %>%
		fill_holes(threshold = area_thresh) %>%
		st_buffer(-2500)
	pa_sp_smooth_move_sim <- as(pa_sf_smooth, "Spatial")
	
	#save(pa_sp_smooth_move_sim, file = "C:/Users/saraw/Documents/SEARRP_Analyses/conn_move_sim/pa_sp_smooth_move_sim.Rdata")
	
	
	
	
	
	
	
	
	
	
	

	
# =============================================================================
#  Get distance (in meters) and ID of the closest PA to each planning unit.
# =============================================================================	

	
	pa_pts_sf <- st_cast(pa_sf_move_sim, "MULTIPOINT") %>%
		mutate(id = 1:n())
	pa_sf <- pa_sf_move_sim %>%
		mutate(id = 1:n())

	start_t <- Sys.time()

		#for(i in 1:34){
			
			tmp1 <- as(st_cast(pa_pts_sf[i,], "POINT"), "Spatial")
			
			for(j in (i+1):35){
			for(j in 24:35){
			#if(j ==8) next

			i = 33
			j = 34
			
				tmp1 <- as(st_cast(pa_pts_sf[i,], "POINT"), "Spatial")
				tmp2 <- as(st_cast(pa_pts_sf[j,], "POINT"), "Spatial")
				
				pts_dist <- as.data.frame(gDistance(tmp1, tmp2,  byid = TRUE))
				#ncol = nrow(tmp1)
				#nrow = nrow(tmp2)
				col_min <- which.min(apply(pts_dist,MARGIN=2,min))
				row_min <- which.min(apply(pts_dist,MARGIN=1,min))
				
				pt1 <- st_as_sf(tmp1[col_min,]) %>%
					mutate(dist_id = paste(i, j, sep = "_")) %>%
					mutate(pt_id = i)
				pt2 <- st_as_sf(tmp2[row_min,]) %>%
					mutate(dist_id = paste(i, j, sep = "_")) %>%
					mutate(pt_id = j)
				
			path_pts <- rbind(path_pts, pt1, pt2)
			
			
			}
			
		#}
	
	end_t <- Sys.time()
	span_t <- end_t - start_t
	span_t
	

	short_path <- path_pts %>% 
		group_by(dist_id) %>% 
		#mutate(pt_dist = st_distance(., lead(.)) %>%
		summarize(m = mean(size)) %>% 
		st_cast("LINESTRING")
paths_fin_tmp1 <- short_path[-16,]
paths_fin_tmp2 <- paths_fin_tmp1[-47,]
paths_fin_tmp3 <- paths_fin_tmp2[-47,]
paths_fin <- paths_fin_tmp3%>%
	separate(dist_id, c("pt1", "pt2")) %>%
	mutate(line_id = paste(pt1, pt2, sep = "_"))
paths_fin$pt1 <- as.numeric(paths_fin$pt1)
paths_fin$pt2 <- as.numeric(paths_fin$pt2)
paths_fin_pt1 <- left_join(paths_fin, pa_df, by = c("pt1" = "id")) %>%
	select(pt1, pt2,line_id, pt1_size = size, geometry)	 %>%
	select(-geometry) 
	
paths_fin_pt2 <- left_join(paths_fin_pt1, pa_df, by = c("pt2" = "id")) %>%
	select(pt1, pt2,line_id, pt1_size, pt2_size = size, geometry.x)	%>%
	mutate(dist_ln = as.numeric(st_length(.))/1000) %>%
	st_buffer(100) %>%
	as.data.frame(.)  
corr_wts_tmp <- paths_fin_pt2 %>%
	rowwise(.) %>%
	mutate(larger_pa = max(pt1_size, pt2_size)) %>%
	mutate(smaller_pa = min(pt1_size, pt2_size)) %>%
	mutate(corr_wt = (1/(log(dist_ln)) + (2*log(smaller_pa)) + (5*log(larger_pa)))) %>%
	as.data.frame()
	
corr_wts <- corr_wts_tmp %>%
	st_as_sf(sf_column_name = "geometry.x")
corr_wts_sp <- as(corr_wts, "Spatial")
corr_r <- rasterize(corr_wts_sp, r_template, field = corr_wts_sp@data$corr_wt)

m <- matrix(1, ncol=30, nrow=30)
corr_sm <- focal(corr_r, m, fun="mean", na.rm=TRUE, NAonly=TRUE, pad=TRUE) 

corr_out <- corr_sm - 60

writeRaster(corr_out, "C:/Users/saraw/Documents/SEARRP_Analyses/conn_move_sim/corr_out.grd")

#writeRaster(corr_sm, "C:/Users/saraw/Documents/SEARRP_Analyses/conn_move_sim/corr_sm.grd")
#writeRaster(corr_r, "C:/Users/saraw/Documents/SEARRP_Analyses/conn_move_sim/corr_r.grd")
	


	
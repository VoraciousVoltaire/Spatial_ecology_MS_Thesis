# Was doing hypothesis 2 again
# There was a problem with the Iberian Coastal months- 4,5,6,7 ; and with Kara Sea- 6,7,8,9,10; and 
# Barent's sea_9

# Checking out Barent's sea 9----
# Couldn't figure out its deal

# for loop for computing input hr rasters for multiplication
  sub <- df %>% filter(LME_NAME == "Barents Sea" & month == 9)
  
    tracks_wgs <- sub
    nrow(tracks_wgs) # 20,381
   
    sub <- subset(sub, lat < 90)  
    png(filename = paste0(dir_kernels, "/distributions/", "Barents Sea_",month_number,".png"))  
    plot(lat~lon, data = sub, type= "n", asp = 1, 
         xlim=c(min(all_data$lon) + 0.2,max(all_data$lon) - 0.2), 
         ylim=c(min(all_data$lat) + 0.5, max(all_data$lat) - 0.5), 
         main="", frame = T, xlab="", ylab="")
    plot(land, col='lightgrey', add= T)
    points(lat~lon, data = sub, pch = 16, cex = 0.5, col="blue")
    points(lat~lon, data = tracks_wgs, pch = 16, cex = 0.5, col="red")
    dev.off()
    
    # Kernel density estimation
    
    if(nrow(tracks_wgs) > 4){
      
      if(min(tracks_wgs$lon) <= -179 ){ lon_min <- -180
      } else {lon_min <- floor(min(tracks_wgs$lon))-1 }
      
      if(max(tracks_wgs$lon) >= 179){ lon_max <- 180
      } else { lon_max <- ceiling(max(tracks_wgs$lon))+1 }
      
      if(min(tracks_wgs$lat) <= -89 ){ lat_min <- -90 
      } else { lat_min <- floor(min(tracks_wgs$lat))-1 }
      
      if(max(tracks_wgs$lat) >= 89){ lat_max <- 90
      } else { lat_max <- ceiling(max(tracks_wgs$lat))+1 }
      
      so.grid <- expand.grid( LON = seq(lon_min, lon_max, by=1), 
                              LAT = seq(lat_min, lat_max, by=1))
      
      sp::coordinates(so.grid) <- ~LON+LAT
      crs(so.grid) <- proj_wgs84
      
      # Setting a colony-centric crs
      mean_loc <- geosphere::geomean(cbind(tracks_wgs$lon,tracks_wgs$lat))
      DgProj <- sp::CRS(paste0("+proj=laea +lon_0=",mean_loc[1],
                               " +lat_0=",mean_loc[2])) 
      
      so.grid.proj <- sp::spTransform(so.grid, CRS = DgProj)
      coords <- so.grid.proj@coords
      
      c <- min(coords[,1])-1000000   ## to check my min lon
      d <- max(coords[,1])+1000000   ## to check my max lon
      
      e <- min(coords[,2])-1000000   ## to check my min lat
      f <- max(coords[,2])+1000000   ## to check my max lat
      
      a <- seq(c, d, by=10000)
      b <- seq(e, f, by=10000)
      null.grid <- expand.grid(x = a,y = b)
      sp::coordinates(null.grid) <- ~x+y
      sp::gridded(null.grid) <- TRUE
      
      # Converting tracks_wgs into a spatial points data frame
      sp::coordinates(tracks_wgs) <- ~lon+lat
      crs(tracks_wgs) <- proj_wgs84
      
      tracks <- sp::spTransform(tracks_wgs, CRS = DgProj)
      tracks$month <- factor(tracks@data$month)
      
      kudl <- adehabitatHR::kernelUD(tracks[,"month"], 
                                     grid = null.grid, h = 200000)  ## smoothing factor equals 200 km for GLS data
      
      vud <- adehabitatHR::getvolumeUD(kudl)
      
      fud <- vud[[1]]
      hr95 <- as.data.frame(fud)[,1]
      hr95 <- as.numeric(hr95 <= 95)
      hr95 <- data.frame(hr95)
      coordinates(hr95) <- coordinates(fud)
      sp::gridded(hr95) <- TRUE
      
      kde_spixdf <- adehabitatHR::estUDm2spixdf(kudl)
      kern95 <- kde_spixdf
      
      stk_100 <- raster::stack(kern95)
      stk_95 <- raster::stack(hr95)
      
      sum_all_100 <- stk_100[[1]]
      sum_all_95 <- stk_95[[1]]
      
      sum_all_raw <- sum_all_100 * sum_all_95
      
      rast <- sum_all_raw/sum(raster::getValues(sum_all_raw))
      rast[rast == 0] <- NA
      
      x.matrix <- is.na(as.matrix(rast))
      colNotNA <- which(colSums(x.matrix) != nrow(rast))
      rowNotNA <- which(rowSums(x.matrix) != ncol(rast))
      
      croppedExtent <- raster::extent(rast, 
                                      r1 = rowNotNA[1]-2, 
                                      r2 = rowNotNA[length(rowNotNA)]+2,
                                      c1 = colNotNA[1]-2, 
                                      c2 = colNotNA[length(colNotNA)]+2)
      
      cropped <- raster::crop(rast, croppedExtent)
      cropped[is.na(cropped)] <- 0
      plot(cropped)
      
      # Changing land's projection for masking
      mask_proj <- sp::spTransform(land, DgProj)
      mask_proj_pol <- as(mask_proj, "SpatialPolygons")
      
      plot(cropped)
      extent(cropped)
      plot(mask_proj_pol, add = T)
      mask_proj_pol_2 <- crop(mask_proj_pol, extent(cropped))
       rast_mask_na <- raster::mask(cropped, mask_proj_pol_2, inverse = T) ### Doesn't work for some reason
      
      rast_mask_na <- raster::mask(cropped, mask_proj_pol, inverse = T)
      rast_mask <- rast_mask_na
      plot(rast_mask)
      
      rast_mask[is.na(rast_mask)] <- 0
      rast_mask_sum1 <- rast_mask/sum(raster::getValues(rast_mask))
      
      rast_mask_final <- raster::mask(rast_mask_sum1, mask_proj_pol, inverse = TRUE)
      rast_mask_final2 <- rast_mask_final 
      
      raster::writeRaster(rast_mask_final, filename = paste0(datadir,"Second_hyp_tifs/",KDE_ref,".tif"), 
                          format = "GTiff", overwrite = TRUE)
      
      mask_wgs84 <- projectRaster(rast_mask_final2, crs = proj_wgs84, over = F)
      
      KDE_ref <- paste0(i, "_", months[month_number])
      
      raster::writeRaster(mask_wgs84, filename = paste0(dir_kernels,"/unique_tifs/",KDE_ref,".tif"), 
                          format = "GTiff", overwrite = TRUE)
      
      png(filename = paste0(dir_kernels,"/unique_distributions/",KDE_ref,".png"))
      plot(mask_wgs84, main = paste0("",i))
      plot(land, add = T, col = "#66000000")
      dev.off()
      
      
    } # First if loop ends
  



    
# Checking out Greenland sea 2---- 
  # Couldn't figure out its deal 
    
    # for loop for computing input hr rasters for multiplication
    sub <- df %>% filter(LME_NAME == "Greenland Sea" & month == 2)
    
    tracks_wgs <- sub
    nrow(tracks_wgs) # 4050
    
    sub <- subset(sub, lat < 90)  
    png(filename = paste0(dir_kernels, "/distributions/", "Barents Sea_",month_number,".png"))  
    plot(lat~lon, data = sub, type= "n", asp = 1, 
         xlim=c(min(all_data$lon) + 0.2,max(all_data$lon) - 0.2), 
         ylim=c(min(all_data$lat) + 0.5, max(all_data$lat) - 0.5), 
         main="", frame = T, xlab="", ylab="")
    plot(land, col='lightgrey', add= T)
    points(lat~lon, data = sub, pch = 16, cex = 0.5, col="blue")
    points(lat~lon, data = tracks_wgs, pch = 16, cex = 0.5, col="red")
    dev.off()
    
    # Kernel density estimation
    
    if(nrow(tracks_wgs) > 4){
      
      if(min(tracks_wgs$lon) <= -179 ){ lon_min <- -180
      } else {lon_min <- floor(min(tracks_wgs$lon))-1 }
      
      if(max(tracks_wgs$lon) >= 179){ lon_max <- 180
      } else { lon_max <- ceiling(max(tracks_wgs$lon))+1 }
      
      if(min(tracks_wgs$lat) <= -89 ){ lat_min <- -90 
      } else { lat_min <- floor(min(tracks_wgs$lat))-1 }
      
      if(max(tracks_wgs$lat) >= 89){ lat_max <- 90
      } else { lat_max <- ceiling(max(tracks_wgs$lat))+1 }
      
      so.grid <- expand.grid( LON = seq(lon_min, lon_max, by=1), 
                              LAT = seq(lat_min, lat_max, by=1))
      
      sp::coordinates(so.grid) <- ~LON+LAT
      crs(so.grid) <- proj_wgs84
      
      # Setting a colony-centric crs
      mean_loc <- geosphere::geomean(cbind(tracks_wgs$lon,tracks_wgs$lat))
      DgProj <- sp::CRS(paste0("+proj=laea +lon_0=",mean_loc[1],
                               " +lat_0=",mean_loc[2])) 
      
      so.grid.proj <- sp::spTransform(so.grid, CRS = DgProj)
      coords <- so.grid.proj@coords
      
      c <- min(coords[,1])-1000000   ## to check my min lon
      d <- max(coords[,1])+1000000   ## to check my max lon
      
      e <- min(coords[,2])-1000000   ## to check my min lat
      f <- max(coords[,2])+1000000   ## to check my max lat
      
      a <- seq(c, d, by=10000)
      b <- seq(e, f, by=10000)
      null.grid <- expand.grid(x = a,y = b)
      sp::coordinates(null.grid) <- ~x+y
      sp::gridded(null.grid) <- TRUE
      
      # Converting tracks_wgs into a spatial points data frame
      sp::coordinates(tracks_wgs) <- ~lon+lat
      crs(tracks_wgs) <- proj_wgs84
      
      tracks <- sp::spTransform(tracks_wgs, CRS = DgProj)
      tracks$month <- factor(tracks@data$month)
      
      kudl <- adehabitatHR::kernelUD(tracks[,"month"], 
                                     grid = null.grid, h = 200000)  ## smoothing factor equals 200 km for GLS data
      
      vud <- adehabitatHR::getvolumeUD(kudl)
      
      fud <- vud[[1]]
      hr95 <- as.data.frame(fud)[,1]
      hr95 <- as.numeric(hr95 <= 95)
      hr95 <- data.frame(hr95)
      coordinates(hr95) <- coordinates(fud)
      sp::gridded(hr95) <- TRUE
      
      kde_spixdf <- adehabitatHR::estUDm2spixdf(kudl)
      kern95 <- kde_spixdf
      
      stk_100 <- raster::stack(kern95)
      stk_95 <- raster::stack(hr95)
      
      sum_all_100 <- stk_100[[1]]
      sum_all_95 <- stk_95[[1]]
      
      sum_all_raw <- sum_all_100 * sum_all_95
      
      rast <- sum_all_raw/sum(raster::getValues(sum_all_raw))
      rast[rast == 0] <- NA
      
      x.matrix <- is.na(as.matrix(rast))
      colNotNA <- which(colSums(x.matrix) != nrow(rast))
      rowNotNA <- which(rowSums(x.matrix) != ncol(rast))
      
      croppedExtent <- raster::extent(rast, 
                                      r1 = rowNotNA[1]-2, 
                                      r2 = rowNotNA[length(rowNotNA)]+2,
                                      c1 = colNotNA[1]-2, 
                                      c2 = colNotNA[length(colNotNA)]+2)
      
      cropped <- raster::crop(rast, croppedExtent)
      cropped[is.na(cropped)] <- 0
      plot(cropped)
      
      # Changing land's projection for masking
      mask_proj <- sp::spTransform(land, DgProj)
      mask_proj_pol <- as(mask_proj, "SpatialPolygons")
      
      plot(cropped)
      extent(cropped)
      plot(mask_proj_pol, add = T)
      mask_proj_pol_2 <- crop(mask_proj_pol, extent(cropped))
      rast_mask_na <- raster::mask(cropped, mask_proj_pol_2, inverse = T) ### Doesn't work for some reason
      
      rast_mask_na <- raster::mask(cropped, mask_proj_pol, inverse = T)
      rast_mask <- rast_mask_na
      plot(rast_mask) # blank
      
      rast_mask[is.na(rast_mask)] <- 0
      rast_mask_sum1 <- rast_mask/sum(raster::getValues(rast_mask))
      
      rast_mask_final <- raster::mask(rast_mask_sum1, mask_proj_pol, inverse = TRUE)
      rast_mask_final2 <- rast_mask_final 
      
      raster::writeRaster(rast_mask_final, filename = paste0(datadir,"Second_hyp_tifs/",KDE_ref,".tif"), 
                          format = "GTiff", overwrite = TRUE)
      
      mask_wgs84 <- projectRaster(rast_mask_final2, crs = proj_wgs84, over = F)
      
      KDE_ref <- paste0(i, "_", months[month_number])
      
      raster::writeRaster(mask_wgs84, filename = paste0(dir_kernels,"/unique_tifs/",KDE_ref,".tif"), 
                          format = "GTiff", overwrite = TRUE)
      
      png(filename = paste0(dir_kernels,"/unique_distributions/",KDE_ref,".png"))
      plot(mask_wgs84, main = paste0("",i))
      plot(land, add = T, col = "#66000000")
      dev.off()
      
      
    } # First if loop ends
    
    
# Checking out Iberian Coastal months 4,5,6 and 7----
# Couldn't figure out its deal but sample size is too small (<50) anyways so no reliable KDE can be estimated 
# for them, hence didn't consider at all. Sample sizes for 4,5,6 and 7 month: 36, 34, 31, 30 respectively
    
    # for loop for computing input hr rasters for multiplication
    sub <- df %>% filter(LME_NAME == "Iberian Coastal" & month == 7)
    
    tracks_wgs <- sub
    nrow(tracks_wgs)
    
    
    sub <- subset(sub, lat < 90)  
    png(filename = paste0(dir_kernels, "/distributions/", "Barents Sea_",month_number,".png"))  
    plot(lat~lon, data = sub, type= "n", asp = 1, 
         xlim=c(min(all_data$lon) + 0.2,max(all_data$lon) - 0.2), 
         ylim=c(min(all_data$lat) + 0.5, max(all_data$lat) - 0.5), 
         main="", frame = T, xlab="", ylab="")
    plot(land, col='lightgrey', add= T)
    points(lat~lon, data = sub, pch = 16, cex = 0.5, col="blue")
    points(lat~lon, data = tracks_wgs, pch = 16, cex = 0.5, col="red")
    dev.off()
    
    # Kernel density estimation
    
    if(nrow(tracks_wgs) > 4){
      
      if(min(tracks_wgs$lon) <= -179 ){ lon_min <- -180
      } else {lon_min <- floor(min(tracks_wgs$lon))-1 }
      
      if(max(tracks_wgs$lon) >= 179){ lon_max <- 180
      } else { lon_max <- ceiling(max(tracks_wgs$lon))+1 }
      
      if(min(tracks_wgs$lat) <= -89 ){ lat_min <- -90 
      } else { lat_min <- floor(min(tracks_wgs$lat))-1 }
      
      if(max(tracks_wgs$lat) >= 89){ lat_max <- 90
      } else { lat_max <- ceiling(max(tracks_wgs$lat))+1 }
      
      so.grid <- expand.grid( LON = seq(lon_min, lon_max, by=1), 
                              LAT = seq(lat_min, lat_max, by=1))
      
      sp::coordinates(so.grid) <- ~LON+LAT
      crs(so.grid) <- proj_wgs84
      
      # Setting a colony-centric crs
      mean_loc <- geosphere::geomean(cbind(tracks_wgs$lon,tracks_wgs$lat))
      DgProj <- sp::CRS(paste0("+proj=laea +lon_0=",mean_loc[1],
                               " +lat_0=",mean_loc[2])) 
      
      so.grid.proj <- sp::spTransform(so.grid, CRS = DgProj)
      coords <- so.grid.proj@coords
      
      c <- min(coords[,1])-1000000   ## to check my min lon
      d <- max(coords[,1])+1000000   ## to check my max lon
      
      e <- min(coords[,2])-1000000   ## to check my min lat
      f <- max(coords[,2])+1000000   ## to check my max lat
      
      a <- seq(c, d, by=10000)
      b <- seq(e, f, by=10000)
      null.grid <- expand.grid(x = a,y = b)
      sp::coordinates(null.grid) <- ~x+y
      sp::gridded(null.grid) <- TRUE
      
      # Converting tracks_wgs into a spatial points data frame
      sp::coordinates(tracks_wgs) <- ~lon+lat
      crs(tracks_wgs) <- proj_wgs84
      
      tracks <- sp::spTransform(tracks_wgs, CRS = DgProj)
      tracks$month <- factor(tracks@data$month)
      
      kudl <- adehabitatHR::kernelUD(tracks[,"month"], 
                                     grid = null.grid, h = 200000)  ## smoothing factor equals 200 km for GLS data
      
      vud <- adehabitatHR::getvolumeUD(kudl)
      
      fud <- vud[[1]]
      hr95 <- as.data.frame(fud)[,1]
      hr95 <- as.numeric(hr95 <= 95)
      hr95 <- data.frame(hr95)
      coordinates(hr95) <- coordinates(fud)
      sp::gridded(hr95) <- TRUE
      
      kde_spixdf <- adehabitatHR::estUDm2spixdf(kudl)
      kern95 <- kde_spixdf
      
      stk_100 <- raster::stack(kern95)
      stk_95 <- raster::stack(hr95)
      
      sum_all_100 <- stk_100[[1]]
      sum_all_95 <- stk_95[[1]]
      
      sum_all_raw <- sum_all_100 * sum_all_95
      
      rast <- sum_all_raw/sum(raster::getValues(sum_all_raw))
      rast[rast == 0] <- NA
      
      x.matrix <- is.na(as.matrix(rast))
      colNotNA <- which(colSums(x.matrix) != nrow(rast))
      rowNotNA <- which(rowSums(x.matrix) != ncol(rast))
      
      croppedExtent <- raster::extent(rast, 
                                      r1 = rowNotNA[1]-2, 
                                      r2 = rowNotNA[length(rowNotNA)]+2,
                                      c1 = colNotNA[1]-2, 
                                      c2 = colNotNA[length(colNotNA)]+2)
      
      cropped <- raster::crop(rast, croppedExtent)
      cropped[is.na(cropped)] <- 0
      plot(cropped)
      
      # Changing land's projection for masking
      mask_proj <- sp::spTransform(land, DgProj)
      mask_proj_pol <- as(mask_proj, "SpatialPolygons")
      
      plot(cropped)
      extent(cropped)
      plot(mask_proj_pol, add = T)
      mask_proj_pol_2 <- crop(mask_proj_pol, extent(cropped))
      rast_mask_na <- raster::mask(cropped, mask_proj_pol_2, inverse = T) ### Doesn't work for some reason
      
      rast_mask_na <- raster::mask(cropped, mask_proj_pol, inverse = T)
      rast_mask <- rast_mask_na
      plot(rast_mask)
      
      rast_mask[is.na(rast_mask)] <- 0
      rast_mask_sum1 <- rast_mask/sum(raster::getValues(rast_mask))
      
      rast_mask_final <- raster::mask(rast_mask_sum1, mask_proj_pol, inverse = TRUE)
      rast_mask_final2 <- rast_mask_final 
      
      raster::writeRaster(rast_mask_final, filename = paste0(datadir,"Second_hyp_tifs/",KDE_ref,".tif"), 
                          format = "GTiff", overwrite = TRUE)
      
      mask_wgs84 <- projectRaster(rast_mask_final2, crs = proj_wgs84, over = F)
      
      KDE_ref <- paste0(i, "_", months[month_number])
      
      raster::writeRaster(mask_wgs84, filename = paste0(dir_kernels,"/unique_tifs/",KDE_ref,".tif"), 
                          format = "GTiff", overwrite = TRUE)
      
      png(filename = paste0(dir_kernels,"/unique_distributions/",KDE_ref,".png"))
      plot(mask_wgs84, main = paste0("",i))
      plot(land, add = T, col = "#66000000")
      dev.off()
      
      
    } # First if loop ends
    

    
# Checking out the deal with Kara Sea months 6,7,8,9 and 10----
# Sample sizes for 6,7,8,9 and 10: 45, 53, 219, 1437, 779 respectively  
      
    # for loop for computing input hr rasters for multiplication
    sub <- df %>% filter(LME_NAME == "Kara Sea" & month == 10)
    
    tracks_wgs <- sub
    nrow(tracks_wgs)
    
    
    sub <- subset(sub, lat < 90)  
    png(filename = paste0(dir_kernels, "/distributions/", "Barents Sea_",month_number,".png"))  
    plot(lat~lon, data = sub, type= "n", asp = 1, 
         xlim=c(min(all_data$lon) + 0.2,max(all_data$lon) - 0.2), 
         ylim=c(min(all_data$lat) + 0.5, max(all_data$lat) - 0.5), 
         main="", frame = T, xlab="", ylab="")
    plot(land, col='lightgrey', add= T)
    points(lat~lon, data = sub, pch = 16, cex = 0.5, col="blue")
    points(lat~lon, data = tracks_wgs, pch = 16, cex = 0.5, col="red")
    dev.off()
    
    # Kernel density estimation
    
    if(nrow(tracks_wgs) > 4){
      
      if(min(tracks_wgs$lon) <= -179 ){ lon_min <- -180
      } else {lon_min <- floor(min(tracks_wgs$lon))-1 }
      
      if(max(tracks_wgs$lon) >= 179){ lon_max <- 180
      } else { lon_max <- ceiling(max(tracks_wgs$lon))+1 }
      
      if(min(tracks_wgs$lat) <= -89 ){ lat_min <- -90 
      } else { lat_min <- floor(min(tracks_wgs$lat))-1 }
      
      if(max(tracks_wgs$lat) >= 89){ lat_max <- 90
      } else { lat_max <- ceiling(max(tracks_wgs$lat))+1 }
      
      so.grid <- expand.grid( LON = seq(lon_min, lon_max, by=1), 
                              LAT = seq(lat_min, lat_max, by=1))
      
      sp::coordinates(so.grid) <- ~LON+LAT
      crs(so.grid) <- proj_wgs84
      
      # Setting a colony-centric crs
      mean_loc <- geosphere::geomean(cbind(tracks_wgs$lon,tracks_wgs$lat))
      DgProj <- sp::CRS(paste0("+proj=laea +lon_0=",mean_loc[1],
                               " +lat_0=",mean_loc[2])) 
      
      so.grid.proj <- sp::spTransform(so.grid, CRS = DgProj)
      coords <- so.grid.proj@coords
      
      c <- min(coords[,1])-1000000   ## to check my min lon
      d <- max(coords[,1])+1000000   ## to check my max lon
      
      e <- min(coords[,2])-1000000   ## to check my min lat
      f <- max(coords[,2])+1000000   ## to check my max lat
      
      a <- seq(c, d, by=10000)
      b <- seq(e, f, by=10000)
      null.grid <- expand.grid(x = a,y = b)
      sp::coordinates(null.grid) <- ~x+y
      sp::gridded(null.grid) <- TRUE
      
      # Converting tracks_wgs into a spatial points data frame
      sp::coordinates(tracks_wgs) <- ~lon+lat
      crs(tracks_wgs) <- proj_wgs84
      
      tracks <- sp::spTransform(tracks_wgs, CRS = DgProj)
      tracks$month <- factor(tracks@data$month)
      
      kudl <- adehabitatHR::kernelUD(tracks[,"month"], 
                                     grid = null.grid, h = 200000)  ## smoothing factor equals 200 km for GLS data
      
      vud <- adehabitatHR::getvolumeUD(kudl)
      
      fud <- vud[[1]]
      hr95 <- as.data.frame(fud)[,1]
      hr95 <- as.numeric(hr95 <= 95)
      hr95 <- data.frame(hr95)
      coordinates(hr95) <- coordinates(fud)
      sp::gridded(hr95) <- TRUE
      
      kde_spixdf <- adehabitatHR::estUDm2spixdf(kudl)
      kern95 <- kde_spixdf
      
      stk_100 <- raster::stack(kern95)
      stk_95 <- raster::stack(hr95)
      
      sum_all_100 <- stk_100[[1]]
      sum_all_95 <- stk_95[[1]]
      
      sum_all_raw <- sum_all_100 * sum_all_95
      
      rast <- sum_all_raw/sum(raster::getValues(sum_all_raw))
      rast[rast == 0] <- NA
      
      x.matrix <- is.na(as.matrix(rast))
      colNotNA <- which(colSums(x.matrix) != nrow(rast))
      rowNotNA <- which(rowSums(x.matrix) != ncol(rast))
      
      croppedExtent <- raster::extent(rast, 
                                      r1 = rowNotNA[1]-2, 
                                      r2 = rowNotNA[length(rowNotNA)]+2,
                                      c1 = colNotNA[1]-2, 
                                      c2 = colNotNA[length(colNotNA)]+2)
      
      cropped <- raster::crop(rast, croppedExtent)
      cropped[is.na(cropped)] <- 0
      plot(cropped)
      
      # Changing land's projection for masking
      mask_proj <- sp::spTransform(land, DgProj)
      mask_proj_pol <- as(mask_proj, "SpatialPolygons")
      
      plot(cropped)
      extent(cropped)
      plot(mask_proj_pol, add = T)
      extent(mask_proj_pol)
      summary(cropped)
      mask_proj_pol_2 <- crop(mask_proj_pol, extent(cropped)) # tried matching extents
      rast_mask_na <- raster::mask(cropped, mask_proj_pol_2, inverse = T) ### Doesn't work for some reason
      
      rast_mask_na <- raster::mask(cropped, mask_proj_pol, inverse = T)
      rast_mask <- rast_mask_na
      plot(rast_mask) # blank
      
      rast_mask[is.na(rast_mask)] <- 0
      rast_mask_sum1 <- rast_mask/sum(raster::getValues(rast_mask))
      
      rast_mask_final <- raster::mask(rast_mask_sum1, mask_proj_pol, inverse = TRUE)
      rast_mask_final2 <- rast_mask_final 
      
      raster::writeRaster(rast_mask_final, filename = paste0(datadir,"Second_hyp_tifs/",KDE_ref,".tif"), 
                          format = "GTiff", overwrite = TRUE)
      
      mask_wgs84 <- projectRaster(rast_mask_final2, crs = proj_wgs84, over = F)
      
      KDE_ref <- paste0(i, "_", months[month_number])
      
      raster::writeRaster(mask_wgs84, filename = paste0(dir_kernels,"/unique_tifs/",KDE_ref,".tif"), 
                          format = "GTiff", overwrite = TRUE)
      
      png(filename = paste0(dir_kernels,"/unique_distributions/",KDE_ref,".png"))
      plot(mask_wgs84, main = paste0("",i))
      plot(land, add = T, col = "#66000000")
      dev.off()
      
      
    } # First if loop ends
    
    
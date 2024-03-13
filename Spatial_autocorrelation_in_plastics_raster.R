# Checking for spatial autocorrelation in plastics raster

library(spData)
library(spdep)
library(raster)
library(sp)





# Replace missing values with 0
plastics_zero <- setValues(plastics, ifelse(is.na(values(plastics)), 0, values(plastics)))

# Convert raster layer to SpatialPointsDataFrame
plastics_sp <- rasterToPoints(plastics_zero)

# Create a spatial weights matrix using queen contiguity
nb <- cell2nb(nrow(plastics_zero), ncol(plastics_zero), type = "queen")
listw <- nb2listw(nb, style = "B", zero.policy = TRUE)

# Check for missing values in the raster object
plastics_values <- getValues(plastics_zero)
if (any(is.na(plastics_values))) {
  stop("Missing values found in the raster data. Please handle missing values.")
}

# Calculate Moran's I
moran_result <- moran.test(plastics_values, listw = listw, alternative = "two.sided")

# Print results
moran_result




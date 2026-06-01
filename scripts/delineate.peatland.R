rm(list = ls())

library(terra)

A <- rast("./data/shapefiles/Dargie_Congo_Peat_LikeliestClassAndProbability.tif")

# 0 - no data
# 1xxx - Water
# 2xxx - Savanna
# 3xxx - Terra firme forest
# 4xxx - Palm-dominated swamp
# 5xxx - Hardwood swamp
B <- round(A/1000)

C <- B
C[B %in% c(0:3)] <- 0
C[B > 3] <- 1

plot(C)

coord <- expand.grid(lon = seq(-179.75,179.75,0.05),
                     lat = seq(-30.25,30.25,0.05)) %>%
  mutate(value = 1)
craster <- rast(raster::rasterFromXYZ(coord))

C.agg <- resample(C,craster)
C.agg[C.agg >= 0.25] <- 1
C.agg[C.agg < 0.25] <- 0
plot(C.agg)

polys <- as.polygons(C.agg, dissolve = TRUE)

# Step 4: Simplify geometry for rough boundary
polys_simplified <- simplifyGeom(polys,
                                 tolerance = 0.01)  # adjust tolerance as needed

# Step 5: Convert to sf and save as shapefile
polys_sf <- st_as_sf(polys_simplified)
names(polys_sf)

# Suppose the column is named "layer" (common with terra::as.polygons)
polys_1 <- polys_sf %>%
  filter(Dargie_Congo_Peat_LikeliestClassAndProbability == 1)
st_crs(polys_1) <- 4326  # EPSG code for WGS84
st_write(polys_1, "./data/shapefiles/Peatland.shp",append = FALSE)
plot(polys_1)



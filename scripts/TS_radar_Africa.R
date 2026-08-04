rm(list = ls())

library(R.matlab)
library(terra)
library(dplyr)
library(ggplot2)
library(lubridate)

# ============================================================
# Load .mat file
# ============================================================

africa <- readMat(
  "~/Downloads/Cband_radar_africa/Merged_LHSCAT_2025_all_afmask_v6.mat"
)

lonlat <- africa$Af.cell.centers
radar  <- africa$Final.time.series.FullRadar

# Fix orientation if needed
if (nrow(radar) == nrow(lonlat)) {
  radar <- t(radar)
}

stopifnot(ncol(radar) == nrow(lonlat))

# ============================================================
# Dates
# ============================================================

ntime <- nrow(radar)

dates <- seq(
  as.Date("1992-01-01"),
  by = "1 month",
  length.out = ntime
)

# ============================================================
# Build raster only over available coordinates
# ============================================================

res <- 0.25

lat <- lonlat[, 1]
lon <- lonlat[, 2]

# Optional: add small buffer around available points
lon_min <- floor(min(lon, na.rm = TRUE) / res) * res
lon_max <- ceiling(max(lon, na.rm = TRUE) / res) * res
lat_min <- floor(min(lat, na.rm = TRUE) / res) * res
lat_max <- ceiling(max(lat, na.rm = TRUE) / res) * res

template <- rast(
  xmin = lon_min,
  xmax = lon_max,
  ymin = lat_min,
  ymax = lat_max,
  resolution = res,
  crs = "EPSG:4326"
)

df <- data.frame(
  lon = lon,
  lat = lat
)

df[paste0("m", seq_len(ntime))] <- as.data.frame(t(radar))

pts <- vect(
  df,
  geom = c("lon", "lat"),
  crs = "EPSG:4326"
)

crast <- rasterize(
  pts,
  template,
  field = paste0("m", seq_len(ntime)),
  fun = mean,
  na.rm = TRUE
)

names(crast) <- format(dates, "%Y_%m")
time(crast) <- dates

# ============================================================
# Quick checks
# ============================================================

print(crast)
print(range(lon, na.rm = TRUE))
print(range(lat, na.rm = TRUE))
print(dates[1])
print(dates[ntime])

plot(crast[[1]])
plot(crast[[ntime]])

# ============================================================
# Climatological mean
# ============================================================

clim <- mean(crast, na.rm = TRUE)

plot(clim)

# ggplot climatology
df_clim <- as.data.frame(clim, xy = TRUE, na.rm = FALSE)
names(df_clim) <- c("lon", "lat", "value")

ggplot(df_clim, aes(lon, lat, fill = value)) +
  geom_raster() +
  coord_equal(expand = FALSE) +
  theme_bw() +
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "value"
  )

# ============================================================
# Regional monthly mean
# ============================================================

monthly_mean <- global(crast, mean, na.rm = TRUE)[, 1]

df_ts <- data.frame(
  date = dates,
  mean = monthly_mean
)

ggplot(df_ts, aes(date, mean)) +
  geom_line() +
  stat_smooth(method = "lm") +
  theme_bw() +
  labs(
    x = "Date",
    y = "Monthly mean",
    title = "Regional monthly mean"
  )

# ============================================================
# Example monthly map
# ============================================================

target_date <- as.Date("2024-01-01")
target_layer <- which(dates == target_date)

plot(crast[[target_layer]])

df_map <- as.data.frame(crast[[target_layer]], xy = TRUE, na.rm = FALSE)
names(df_map) <- c("lon", "lat", "value")

ggplot(df_map, aes(lon, lat, fill = value)) +
  geom_raster() +
  coord_equal(expand = FALSE) +
  theme_bw() +
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "value",
    title = paste("Radar monthly value —", target_date)
  )

# ============================================================
# Save raster
# ============================================================

dir.create("./outputs", showWarnings = FALSE)

writeRaster(
  crast,
  "./outputs/Radar_all.years.tif",
  overwrite = TRUE,
  gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "TFW=YES")
)

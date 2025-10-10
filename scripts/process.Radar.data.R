rm(list = ls())

library(R.matlab)
library(ggplot2)
library(terra)
library(raster)
library(reshape2)
library(dplyr)
library(lubridate)

# ---- Load .mat ----
africa <- readMat('~/Downloads/Africa_Final_full_radar_YYT.mat')  # 1992–Sep 2024

# Input data (cell centers & time series)
lonlat.mask <- africa$Af.cell.centers          # [ncell x 2] (lon, lat) in degrees
radar.mask  <- africa$Final.time.series.FullRadar  # [ntime x ncell], ntime = 393 (Jan 1992–Sep 2024)

# ---- Target grid: lon -25..65, lat -25..25 at 0.25° ----
# Grid meta
res      <- 0.25
lon_min  <- -25
lon_max  <-  65
lat_min  <- -25
lat_max  <-  25

# number of cols/rows
ncol <- as.integer((lon_max - lon_min) / res)  # 90 / 0.25 = 360
nrow <- as.integer((lat_max - lat_min) / res)  # 50 / 0.25 = 200
ntime <- nrow(radar.mask)                      # should be 393

# cell-center origins used by your indexing scheme
lon0_center <- lon_min + res/2   # -24.875
lat0_top    <- lat_max - res/2   # +24.875 (top row center)

# Allocate arrays
radar.africa.data  <- array(0,  dim = c(ncol, nrow, ntime))
radar.africa.count <- array(0,  dim = c(ncol, nrow, ntime))

# ---- Bin each station time series into the 0.25° grid ----
num.pixels <- nrow(lonlat.mask)
for (rr in 1:num.pixels) {
  this_lon <- lonlat.mask[rr, 1]
  this_lat <- lonlat.mask[rr, 2]

  # column index (1..ncol) using cell centers
  lonlon <- floor((this_lon - lon0_center) / res) + 1L

  # row index (1..nrow), top-to-bottom (row 1 is highest latitude),
  # consistent with your original lat indexing
  latlat <- floor((lat0_top - this_lat) / res) + 1L

  if (lonlon > 0 && lonlon <= ncol && latlat > 0 && latlat <= nrow) {
    # add the full time series for this cell
    radar.africa.data[lonlon, latlat, ]  <- radar.africa.data[lonlon, latlat, ]  + radar.mask[, rr]
    radar.africa.count[lonlon, latlat, ] <- radar.africa.count[lonlon, latlat, ] + as.numeric(!is.na(radar.mask[, rr]))
  }
}

# Avoid 0/0; convert zeros to NA then compute grid means
radar.africa.data[radar.africa.data == 0]   <- NA
radar.africa.count[radar.africa.count == 0] <- NA
radar.africa.grid <- radar.africa.data / radar.africa.count   # [ncol x nrow x ntime]

# ---- Reshape to monthly [year, month] cubes as you had ----
# Years 1992..2023 (32 years), plus Jan–Sep 2024 (9 months) = 393
africa.radar.monthly <- array(NA_real_, dim = c(ncol, nrow, 33, 12))

# Flip latitude for storage (so 1..nrow becomes south->north later when plotting)
for (year in 1992:2023) {
  idx <- ((year - 1992) * 12 + 1):((year - 1992) * 12 + 12)
  month.radar.grid <- radar.africa.grid[, nrow:1, idx]
  africa.radar.monthly[, , year - 1991, ] <- month.radar.grid
}

# 2024 (first 9 months)
year <- 2024
africa.radar.monthly[, , 33, 1:9] <- radar.africa.grid[, nrow:1, ((year - 1992) * 12 + 1):ntime]

# Annual means per cell per year, then climatological mean
africa.radar.annual <- apply(africa.radar.monthly, c(1, 2, 3), mean, na.rm = TRUE)
africa.mean         <- apply(africa.radar.annual,  c(1, 2),     mean, na.rm = TRUE)  # [ncol x nrow]

# ---- Build lon/lat for plotting (use true coordinates) ----
lon_vec <- seq(lon_min + res/2, lon_max - res/2, by = res)     # length = ncol
lat_vec <- seq(lat_min + res/2, lat_max - res/2, by = res)     # length = nrow

# Note: africa.mean is [lon, lat_flipped]; we flipped rows earlier (nrow:1),
# so align with lat_vec accordingly:
df2plot <- melt(africa.mean, varnames = c("ix", "iy"), value.name = "value")
df2plot$lon <- lon_vec[df2plot$ix]
df2plot$lat <- lat_vec[df2plot$iy]

# ---- Plot ----
ggplot(df2plot,
       aes(x = lon, y = lat, fill = value)) +
  geom_raster() +
  coord_fixed(expand = FALSE) +
  labs(x = "Longitude (°E)", y = "Latitude (°)", fill = "Value") +
  theme_bw()

monthly_mean <- apply(radar.africa.grid, 3, mean, na.rm = TRUE)  # length 393

# 2) Build the date index (Jan 1992 → Sep 2024 = 393 months)
dates <- seq(as.Date("1992-01-01"), by = "1 month", length.out = length(monthly_mean))

df_ts <- data.frame(date = dates, mean = as.numeric(monthly_mean))

ggplot(df_ts, aes(x = date, y = mean)) +
  geom_line() +
  labs(x = "Date", y = "Monthly mean (region-wide)", title = "Regional Monthly Mean") +
  theme_bw()


df_monthly <- reshape2::melt(
  africa.radar.monthly,
  varnames = c("lon_idx","lat_idx","year_idx","month"),
  value.name = "value")

# Map indices -> coordinates & calendar
df_monthly$lon   <- lon_vec[df_monthly$lon_idx]
df_monthly$lat   <- lat_vec[df_monthly$lat_idx]
df_monthly$year  <- 1992 + df_monthly$year_idx - 1L
df_monthly$month <- as.integer(df_monthly$month)

# Example: map for June 2024
ggplot(df_monthly %>%
         filter(year == 2013, month == 12),
       aes(lon, lat, fill = value)) +
  geom_raster() +
  coord_equal(expand = FALSE) +
  labs(x = "Longitude (°E)", y = "Latitude (°)", fill = "Value",
       title = "Radar monthly value — June 2024") +
  theme_bw()


times <- df_monthly %>%
  dplyr::select(year,month) %>%
  distinct() %>%
  arrange(year,month)

dates <- as.Date(sprintf("%04d-%02d-01", times$year, times$month))
all.rast <- list()

for (itime in seq(1,nrow(times))){

  print(itime)
  cyear <- times$year[itime]
  cmonth <- times$month[itime]
  cdf <- df_monthly %>%
    filter(year == cyear, month == cmonth)
  all.rast[[itime]] <-
    rast(rasterFromXYZ(cdf %>%
                         dplyr::select(lon,lat,value)))

}

crast <- rast(all.rast)
time(crast) <- dates

writeRaster(crast,
            paste0("./outputs/","Radar_all.years.tif"),
            overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

plot(crast[[392]])

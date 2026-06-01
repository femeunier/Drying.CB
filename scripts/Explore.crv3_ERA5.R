rm(list = ls())

library(ncdf4)
library(lubridate)
library(dplyr)
library(ggplot2)

nc <- nc_open("~/Downloads/20crv3-era5_obsclim_tas_global_daily_2022_2024.nc")

lon  <- ncvar_get(nc, "lon")
lat  <- ncvar_get(nc, "lat")
time <- ncvar_get(nc, "time")

# Target coordinate
target_lon <- 20
target_lat <- -2

# Find closest grid cell
ilon <- which.min(abs(lon - target_lon))
ilat <- which.min(abs(lat - target_lat))

lon_sel <- lon[ilon]
lat_sel <- lat[ilat]

# Extract tas time series at that cell
tas <- ncvar_get(
  nc,
  "tas",
  start = c(ilon, ilat, 1),
  count = c(1, 1, -1)
)

nc_close(nc)

# Convert time
dates <- as.Date(time, origin = "1860-01-01")

# Kelvin to °C
tas_c <- tas - 273.15

df <- data.frame(
  date = dates,
  tas_c = as.numeric(tas_c)
)

ggplot(df, aes(date, tas_c)) +
  geom_line() +
  labs(
    title = paste0("Daily near-surface air temperature at closest grid cell"),
    subtitle = paste0("Target: 0.5°N, 25.5°E | Grid cell: ", lat_sel, "°N, ", lon_sel, "°E"),
    x = NULL,
    y = "Temperature (°C)"
  ) +
  theme_minimal()

################################################################################

nc <- nc_open("~/Downloads/20crv3-era5_obsclim_pr_global_daily_2022_2024.nc")

lon  <- ncvar_get(nc, "lon")
lat  <- ncvar_get(nc, "lat")
time <- ncvar_get(nc, "time")

ilon <- which.min(abs(lon - target_lon))
ilat <- which.min(abs(lat - target_lat))

lon_sel <- lon[ilon]
lat_sel <- lat[ilat]

pr <- ncvar_get(
  nc,
  "pr",
  start = c(ilon, ilat, 1),
  count = c(1, 1, -1)
)

nc_close(nc)

dates <- as.Date(time, origin = "1860-01-01")

df_daily <- data.frame(
  date = dates,
  pr_mm_day = as.numeric(pr) * 86400
)

df_monthly <- df_daily %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>%
  summarise(
    pr_mm_month = sum(pr_mm_day, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(df_monthly, aes(month, pr_mm_month)) +
  geom_col() +
  labs(
    title = "Monthly precipitation",
    subtitle = paste0(
      "Target: 0.5°N, 25.5°E | Grid cell: ",
      lat_sel, "°N, ", lon_sel, "°E"
    ),
    x = NULL,
    y = "Precipitation (mm/month)"
  ) +
  theme_minimal()

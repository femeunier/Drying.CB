rm(list = ls())

library(ncdf4)
library(dplyr)
library(terra)
library(lubridate)

ncfiles <- list.files("/data/gent/vo/000/gvo00074/ED_common_data/met/CB/ERA5/",
                      pattern = "*pressure*",full.names = TRUE)

df.ts <- data.frame()
monthly.maps <- list()

for (ncfile in ncfiles){

  print(basename(ncfile))
  nc <- nc_open(ncfile)

  data <- ncvar_get(nc,"w")
  times <- as.Date(ncvar_get(nc,"valid_time")/86400,
                   origin = "1970-01-01")
  months = month(times)
  years = year(times)

  lats <- ncvar_get(nc,"latitude")
  lons <- ncvar_get(nc,"longitude")

  monthly.m <- array(
    unlist(lapply(1:12, function(m) {
      apply(data[, , months == m, drop = FALSE], c(1, 2), mean, na.rm = TRUE)
    })),
    dim = c(length(lons), length(lats), 12)
  )

  if (length(lats) > 1 && lats[1] > lats[length(lats)]) {
    monthly.m <- monthly.m[, length(lats):1, , drop = FALSE]
    lats <- rev(lats)
  }

  r_monthly <- rast(
    monthly.m,
    ext = ext(min(lons), max(lons), min(lats), max(lats)),
    crs = "EPSG:4326"
  )

  time(r_monthly) <- as.Date(paste0(unique(years),"-",1:12,"-01"))
  monthly.maps[[basename(ncfile)]] <- r_monthly


  timeseries <- data[which.min(abs(lons - 25.5)),
                     which.min(abs(lats - 0.5)),]

  df.ts <- bind_rows(df.ts,
                     data.frame(time = times,
                                year = years,
                                month = months,
                                value = timeseries))
  nc_close(nc)
}


saveRDS(rast(monthly.maps),"./outputs/Monthly.maps.Congo.wp.RDS")
saveRDS(df.ts,"./outputs/ts.Congo.wp.RDS")


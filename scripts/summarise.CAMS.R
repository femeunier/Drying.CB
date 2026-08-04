rm(list = ls())

# https://psl.noaa.gov/data/gridded/data.ghcncams.html

library(reshape2)
library(dplyr)
library(lubridate)

ncfile1 <- "/data/gent/vo/000/gvo00074/ED_common_data/met/CAMS/air.mon.mean.nc"

nc1 <- nc_open(ncfile1)

climatology <- ncvar_get(nc1,"air")

times <- as.Date(ncvar_get(nc1,"time")/24,
                 origin = "1800-1-1 00:00:00")
lats <- ncvar_get(nc1,"lat")
lons <- ncvar_get(nc1,"lon")
lons[lons > 180] <- lons[lons > 180]-360

nc_close(nc1)

df <- melt(climatology) %>%
  filter(!is.na(value)) %>%
  rename(lon = Var1,
         lat = Var2,
         time = Var3) %>%
  mutate(lon = lons[lon],
         lat = lats[lat],
         time = times[time]) %>%
  mutate(year = year(time),
         month = month(time)) %>%
  dplyr::select(-time)

df.select <- df %>%
  filter(abs(lat) <= 30) %>%
  # filter(year >= 1980) %>%
  ungroup()

saveRDS(df.select %>%
          rename(tas = value) %>%
          filter(!is.na(tas)),
        "/data/gent/vo/000/gvo00074/felicien/R/outputs/df.CAMS.Tropics.climate.RDS")

# r <- rasterFromXYZ(df %>%
#                      filter(year == year[1],
#                             month == 1) %>%
#                      dplyr::select(lon,lat,value))

# scp /Users/felicien/Documents/projects/Drying.CB/scripts/summarise.CAMS.R hpc:/data/gent/vo/000/gvo00074/felicien/R


# saveRDS(df.select %>%
#           rename(tas = value) %>%
#           filter(!is.na(tas)),
#         "~/Downloads/df.CAMS.Tropics.climate.RDS")
#
# system2("rsync",
#         c("-avz",
#           "~/Downloads/df.CAMS.Tropics.climate.*",
#           "hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/"))

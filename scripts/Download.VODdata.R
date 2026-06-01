rm(list = ls())

library(ncdf4)
library(terra)
library(raster)
library(reshape2)

file <- "/home/femeunier/Downloads/SM_SCIE_MIR_L4AGBC_20110101T000000_20241231T235959_100_001_8.nc"

nc <- nc_open(file)
lats <- ncvar_get(nc,"lat")
lons <- ncvar_get(nc,"lon")
times <- as.Date(paste0(ncvar_get(nc,"year"),"/01/01"))

AGBs <- ncvar_get(nc,"AGB")

df <- melt(AGBs) %>%
  rename(lon = Var1,
         lat = Var2,
         time = Var3) %>%
  mutate(lon = lons[lon],
         lat = lats[lat],
         time = times[time]) %>%
  na.omit()

all.times <- sort(unique(df$time))

all.R <- list()

rx <- ry <- 1

for (itime in seq(1,length(all.times))){

  print(itime)

  ctime <- all.times[itime]
  cdf <- df %>%
    filter(time == ctime)

  cdf <- subset(cdf,
                is.finite(lon) & is.finite(lat) & is.finite(value))
  cdf <- aggregate(value ~ lon + lat, cdf, mean)  # collapse duplicates

  extnt <- ext(min(cdf$lon), max(cdf$lon), min(cdf$lat), max(cdf$lat))
  r0 <- rast(extnt, resolution = c(rx, ry), crs = "EPSG:4326")

  xy <- cbind(cdf$lon, cdf$lat)
  r  <- crop(rasterize(xy, r0, cdf$value, fun = mean, background = NA),
             extent(-180,180,-25,25))

  all.R[[itime]] <- r

}

all.Rs <- rast(all.R)
time(all.Rs) <- all.times

Mask <- read_sf("./data/shapefiles/Rainforests.shp")
cdata.mask <- crop(mask(all.Rs,
                        Mask),
                   ext(-25,65,-25,25))

N <- global(all.Rs,function(i)
  length(which(!is.na(as.vector(i)))))

plot(N$global)

ts <- global(all.Rs,mean,na.rm = TRUE)
plot(ts$mean)

ts.mask <- global(cdata.mask,mean,na.rm = TRUE)
plot(all.times,
     ts.mask$mean,type = "p")

plot((cdata.mask[[14]] - cdata.mask[[1]])/
       cdata.mask[[1]])

nc_close(nc)

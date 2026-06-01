rm(list = ls())

library(terra)
library(dplyr)

vpd.files <- list.files("/data/gent/vo/000/gvo00074/felicien/R/outputs/all.climate",
                        pattern = ".*vpd_all.years.tif",full.names = TRUE)

clon = 24.82 ; clat = 0.82
Ygb <- vect(data.frame(
  lon = c(clon),
  lat = c(clat)
), geom = c("lon", "lat"), crs = "EPSG:4326")


all.ts <- data.frame()
for (cfile in vpd.files){

  cR <- rast(cfile)
  cts <- as.vector(as.numeric(extract(cR,Ygb)))

  cdf.ts <- data.frame(time = time(cR),
                       vpd = (cts)[2:length(cts)])

  all.ts <- bind_rows(all.ts,
                      cdf.ts %>%
                        mutate(source = strsplit(basename(cfile),"\\_")[[1]][1]))


}

saveRDS(all.ts,
        "./outputs/VPD.YGB.RDS")

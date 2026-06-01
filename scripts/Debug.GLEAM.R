rm(list = ls())

library(terra)
library(sf)

Mask <- read_sf("./data/shapefiles/Rainforests.shp")

all.data <- rast("~/Downloads/SMrz_all_years.tif")
all.data.msk <- crop(mask(all.data,
                     Mask),ext(-25,65,-25,25))
GLEAM <- crop(all.data,
              ext(-180,180,-25,25))

ts <- global(GLEAM,mean,na.rm = TRUE)

plot(ts$mean[517:540],type = "l")

GLEAM.2 <- GLEAM[[which(names(GLEAM) == "2024_1")]]
GLEAM.1 <- GLEAM[[which(names(GLEAM) == "2023_12")]]
diff <- GLEAM.2 - GLEAM.1
plot(diff)
hist(as.vector(diff))

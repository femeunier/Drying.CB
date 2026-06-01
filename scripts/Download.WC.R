rm(list = ls())

library(terra)
library(geodata)

setwd("/Users/felicien/Documents/projects/Drying.CB")

worldclim_global("tmin", 0.5,
                 "./data/", version="2.1")
worldclim_global("tmax", 0.5,
                 "./data/", version="2.1")
worldclim_global("tavg", 0.5,
                 "./data/", version="2.1")
worldclim_global("prec", 0.5,
                 "./data/", version="2.1")


A <- rast("./data/Climate/wc2.1_5m/wc2.1_5m_tmin_01.tif")
plot(crop(A,
          ext(7,33,-7,7)))


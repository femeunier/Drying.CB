




system2("rsync",
        c("-avz","hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/all.climate/GPCC_pre_all.years.tif",
          "./outputs/"))

Mask <- read_sf("~/Documents/projects/Congo.vs.Amazon/data/Rainforests.shp")
A <- crop(mask(rast("./outputs/GPCC_pre_all.years.tif"),
          Mask),ext(Mask))

times = terra::time(A)
plot(A[[times == as.Date("2024-08-01")]])

B <- crop(mask(rast("./outputs/Maps/GPCC_pre_Zanomalies_CA.tif"),
               Mask),ext(Mask))
plot(B[[times == as.Date("2024-12-01")]])

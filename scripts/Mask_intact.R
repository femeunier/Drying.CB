rm(list = ls())

library(terra)

E <- ext(-25, 65, -25, 15)

# TMF
files <- list.files("/data/gent/vo/000/gvo00074/felicien/TMF/",
                    pattern = "\\.tif$",
                    full.names = TRUE)
undisturbed <- crop(vrt(files),E)
is.intact_30m <- (undisturbed == 1)

# Radar grid
Radar.data <- rast("/data/gent/vo/000/gvo00074/felicien/R/data/Radar/Radar_all.years.tif")[[12]]
Radar.data <- crop(Radar.data, E)
radar_template <- Radar.data
values(radar_template) <- NA

# CCI
LC <- rast("/data/gent/vo/000/gvo00074/felicien/ELIE/C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1.nc")[[1]]
LC <- crop(LC, E)
LC.evergreen <- ifel(LC %in% c(50, 160, 170), 1, 0)

writeRaster(LC.evergreen,
            paste0("/data/gent/vo/000/gvo00074/felicien/R/data/LC.evergreen.tif"),
            overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

# Fractions within radar pixels
evergreen_prop <- project(LC.evergreen, radar_template, method = "average")

# intact_prop    <- project(is.intact_30m, radar_template, method = "average")
intact_prop <- resample(is.intact_30m, radar_template, method = "average")


writeRaster(evergreen_prop,
            paste0("/data/gent/vo/000/gvo00074/felicien/R/data/evergreen_prop.tif"),
            overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

writeRaster(intact_prop,
            paste0("/data/gent/vo/000/gvo00074/felicien/R/data/intact_prop.tif"),
            overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

# Final mask
is.evergreen_radar <- (evergreen_prop > 0.5)
is.intact_radar    <- (intact_prop >= 0.95)
mask_final         <- (is.evergreen_radar & is.intact_radar)

writeRaster(is.evergreen_radar,
            paste0("/data/gent/vo/000/gvo00074/felicien/R/data/is.evergreen_radar.tif"),
            overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

writeRaster(is.intact_radar,
            paste0("/data/gent/vo/000/gvo00074/felicien/R/data/is.intact_radar.tif"),
            overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

# Apply to radar data
Radar.intact <- ifel(mask_final, Radar.data, NA)

plot(mask_final)
plot(Radar.intact)


writeRaster(Radar.intact,
            paste0("./data/Intact_Radar_mask.tif"),
            overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

writeRaster(mask_final,
            paste0("./data/Intact_mask.tif"),
            overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

###############################################################################

system2("scp",
        c("hpc:/data/gent/vo/000/gvo00074/felicien/R/data/evergreen_prop.tif",
          "./data/shapefiles/"))

system2("scp",
        c("hpc:/data/gent/vo/000/gvo00074/felicien/R/data/intact_prop.tif",
          "./data/shapefiles/"))

system2("scp",
        c("hpc:/data/gent/vo/000/gvo00074/felicien/R/data/is.evergreen_radar.tif",
          "./data/shapefiles/"))

system2("scp",
        c("hpc:/data/gent/vo/000/gvo00074/felicien/R/data/is.intact_radar.tif",
          "./data/shapefiles/"))

system2("scp",
        c("hpc:/user/gent/425/vsc42558/R/data/Intact_Radar_mask.tif",
          "./data/shapefiles/"))

system2("scp",
        c("hpc:/user/gent/425/vsc42558/R/data/Intact_mask.tif",
          "./data/shapefiles/"))

evergreen_prop <- rast("./data/shapefiles/evergreen_prop.tif")
plot(evergreen_prop)
intact_prop <- rast("./data/shapefiles/intact_prop.tif")
plot(intact_prop)

is.evergreen_radar <- rast("./data/shapefiles/is.evergreen_radar.tif")
plot(is.evergreen_radar)
is.intact_radar <- rast("./data/shapefiles/is.intact_radar.tif")
plot(is.intact_radar)

Intact_Radar_mask <- rast("./data/shapefiles/Intact_Radar_mask.tif")
plot(Intact_Radar_mask)
Intact_mask <- rast("./data/shapefiles/Intact_mask.tif")
plot(Intact_mask)


is.intact_radar    <- (intact_prop >= 0.9)
mask_final_mod         <- (is.evergreen_radar & is.intact_radar)
plot(mask_final_mod)



writeRaster(mask_final_mod,
            paste0("./data/shapefiles/Intact_mask_mod.tif"),
            overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

# scp /Users/felicien/Documents/projects/Drying.CB/scripts/Mask_intact.R hpc:/data/gent/vo/000/gvo00074/felicien/R/
# scp /Users/felicien/Documents/projects/Drying.CB/scripts/Mask_intact.R hydra:/user/gent/425/vsc42558/R



rm(list = ls())

library(terra)

E <- ext(-180, 180, -25, 25)

# CCI
LC <- rast("~/Downloads//C3S-LC-L4-LCCS-Map-300m-P1Y-2022-v2.1.1.nc")[[3]]
LC <- crop(LC, E)
LC.evergreen <- ifel(LC %in% c(50, 160, 170), 1, 0)[[1]]

template <- rast(
  ext = align(ext(E), 0.1),
  resolution = 0.1,
  crs = "EPSG:4326"
)

# LC.evergreen.01 <- project(LC.evergreen,
#                            template,
#                            method = "average")

fact <- round(0.1 / res(LC))

LC.evergreen.01 <- aggregate(
  LC.evergreen,
  fact = fact,
  fun = "mean",
  na.rm = TRUE
)

LC.evergreen.01.mask <- (LC.evergreen.01 >= 0.5)

writeRaster(
  LC.evergreen.01.mask,
  filename = "./outputs/LC_evergreen_01deg_mask.tif",
  datatype = "INT1U",
  NAflag = 255,
  overwrite = TRUE,
  gdal = c("COMPRESS=DEFLATE", "TILED=YES")
)

plot(LC.evergreen.01.mask)

system2("scp",
        c("./outputs/LC_evergreen_01deg_mask.tif",
          "hpc:/kyukon/data/gent/vo/000/gvo00074/felicien/R/data/shapefiles/"))


# scp /Users/felicien/Documents/projects/Drying.CB/scripts/Mask_intact.R hpc:/data/gent/vo/000/gvo00074/felicien/R/
# scp /Users/felicien/Documents/projects/Drying.CB/scripts/Mask_intact.R hydra:/user/gent/425/vsc42558/R



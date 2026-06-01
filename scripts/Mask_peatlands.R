rm(list = ls())

library(terra)

A <- rast("~/Downloads/Crezee_2022_Smoothed_7x7_Classification_Most_likely_class.tif")
A[A<4] <- 0
A[A>=4] <- 1

plot(A)

writeRaster(A,
            paste0("./data/shapefiles/Peatlands_mask.tif"),
            overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))


A_agg <- aggregate(A, fact = 250,
                   fun = "modal", na.rm = TRUE)
A_agg[A_agg == 0] <- NA
plot(A_agg)
v1 <- as.polygons(A_agg, dissolve = TRUE, na.rm = TRUE)
plot(v1)
writeVector(v1, "./data/shapefiles/Peatlands_mask.shp", overwrite = TRUE)

system2("rsync",
        c("-avz",
          "./data/shapefiles/Peatlands_mask*",
          "hpc:/kyukon/data/gent/vo/000/gvo00074/felicien/R/data/"))

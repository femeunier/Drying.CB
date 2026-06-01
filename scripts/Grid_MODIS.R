rm(list = ls())

library(terra)

A <- rast("~/Downloads/MOD13C2.A2004306.061.2020211202721.hdf")
plot(A[[1]])

B <- crop(A[[1]],ext(-25,65,-20,15))
plot(B)

writeRaster(B,
            paste0("~/Downloads/Grid.tif"),
            overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))
plot(B)

target <- rast(
  ext(B),
  resolution = 0.01/2,
  crs = crs(B)
)

# Resample
B_001 <- resample(B, target, method = "bilinear")**0
B_001 <- !is.na(B_001)   # TRUE/FALSE grid
B_001 <- as.int(B_001)   # 0/1


writeRaster(
  B_001,
  "~/Downloads/Grid_fine.tif",
  overwrite = TRUE,
  datatype = "INT1U",
  gdal = c(
    "COMPRESS=DEFLATE",
    "PREDICTOR=2",
    "ZLEVEL=9",
    "TILED=YES",
    "BIGTIFF=YES"
  )
)
plot(B_001)

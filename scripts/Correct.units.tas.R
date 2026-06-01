rm(list = ls())

library(raster)
library(terra)

dir <- "/data/gent/vo/000/gvo00074/felicien/R/outputs/all.climate"
files <- list.files(dir,
                    pattern = "*.tif",
                    full.names = TRUE)

files.split <- strsplit(basename(files),"\\_")
models <- (sapply(files.split,"[[",1))
models.uni <- sort(unique(models))

vars <- sapply(files.split,"[[",2)
vars.uni <- c("tas","tasmin","tasmax")

for (cvar in vars.uni){

  for (cmodel in models.uni){

    print(paste0(cmodel,"-",cvar))

    cfiles <- files[which((vars == cvar) & (models == cmodel))]

    if (length(cfiles) == 0){ next()}

    cR <- rast(cfiles)

    if (global(cR[[1]],mean,na.rm = TRUE)[["mean"]] > 200){
      cR <- cR -273.15

      print("Correcting")

      writeRaster(cR,
                  paste0("/data/gent/vo/000/gvo00074/felicien/R/outputs/all.climate/",
                         cmodel,"_",cvar,"_all.years.tif"),
                  overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

    }
  }
}

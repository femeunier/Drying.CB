rm(list = ls())

library(raster)
library(terra)

dir <- "/data/gent/vo/000/gvo00074/ED_common_data/met/Precip.Tropics/CMIP6/"
files <- list.files(dir,
                    pattern = "*.tif",
                    full.names = TRUE)
files <- files[!grepl("MEM",files)]

files.split <- strsplit(basename(files),"\\.")
models <- sapply(files.split,"[[",3)
scenarios <- sapply(files.split,"[[",4)
scenarios.uni <- sort(unique(scenarios))

vars <- gsub("_CA","",sapply(files.split,"[[",5))
vars.uni <- sort(unique(vars))

for (cscenario in scenarios.uni){
  for (cvar in vars.uni){

    cfiles <- files[which(vars == cvar & scenarios == cscenario)]
    cR <- rast(cfiles)

    all.times <- time(cR)
    ctimes <- sort(unique(all.times))

    MEM <- list()
    for (itime in seq(1,length(ctimes))){

      ctime <- ctimes[itime]
      print(paste0(cscenario,"-",cvar,"-",ctime))
      pos <- which(all.times == ctime)
      MEM[[itime]] <- mean(cR[[pos]])

    }

    all.MEM <- rast(MEM)
    time(all.MEM) <- ctimes

    writeRaster(all.MEM,
                paste0("/data/gent/vo/000/gvo00074/ED_common_data/met/Precip.Tropics/CMIP6/",
                       "Timeseries.MCWD.MEM.",cscenario,".",cvar,"_CA.tif"),
                overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

  }
}


# scp /home/femeunier/Documents/projects/Drying.CB/scripts/Compute.MEM.CMIP6.R hpc:/data/gent/vo/000/gvo00074/felicien/R/

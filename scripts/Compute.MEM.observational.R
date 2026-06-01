rm(list = ls())

library(raster)
library(terra)

dir <- "/data/gent/vo/000/gvo00074/felicien/R/outputs/all.climate"
files <- list.files(dir,
                    pattern = "*.tif",
                    full.names = TRUE)
files <- files[!grepl("MEM",files)]

files.split <- strsplit(basename(files),"\\_")
models <- sapply(files.split,"[[",1)
vars <- sapply(files.split,"[[",2)

vars.uni <- sort(unique(vars))
vars.uni <- c("pre",'tasmin',"vpd")

for (cvar in vars.uni){

  cfiles <- files[which(vars == cvar)]
  #cR <- rast(cfiles)

  rlist_raw <- lapply(cfiles, rast)
  common_ext <- Reduce(intersect, lapply(rlist_raw, ext))
  template <- crop(rlist_raw[[2]], common_ext)
  rlist <- lapply(rlist_raw, function(r) {
    r_crop <- crop(r, common_ext)
    resample(r_crop, template, method = "bilinear")
  })

  cR <- rast(rlist)

  all.times <- time(cR)
  ctimes <- sort(unique(all.times))

  MEM <- list()
  for (itime in seq(1,length(ctimes))){

    ctime <- ctimes[itime]
    print(paste0(cvar,"-",ctime))
    pos <- which(all.times == ctime)
    MEM[[itime]] <- mean(cR[[pos]],na.rm = TRUE)

  }

  all.MEM <- rast(MEM)
  time(all.MEM) <- ctimes

  writeRaster(all.MEM,
              paste0("/data/gent/vo/000/gvo00074/felicien/R/outputs/all.climate/",
                     "MEM_",cvar,"_all.years.tif"),
              overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

}

# scp /Users/felicien/Documents/projects/Drying.CB/scripts/Compute.MEM.observational.R hpc:/data/gent/vo/000/gvo00074/felicien/R/


rm(list = ls())

library(raster)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(lubridate)
library(TrENDY.analyses)
library(tidyr)
library(raster)
library(terra)

files <- list.files("/data/gent/vo/000/gvo00074/ED_common_data/met/Precip.Tropics/CMIP6/",
                    pattern = "*.RDS",
                    full.names = TRUE)

overwrite = TRUE

for (ifile in seq(1,length(files))){

  cfile <- files[ifile]
  cmodel <- strsplit(basename(cfile),
                     "\\.")[[1]][3]
  cscenario <-strsplit(basename(cfile),
                       "\\.")[[1]][4]
  print(paste0(cmodel," - ",cscenario))



  cdf <- readRDS(cfile) %>%
    ungroup()
  cdf.selected <- cdf %>%
    filter(abs(lat) <= 25,
           lon > -25, lon <= 65)

  cn <- colnames(cdf.selected)
  cn2keep <- cn[!(cn %in% c("year","month","lon","lat",
                            "scenario","variant","model",
                            "period"))]

  times <- cdf.selected %>%
    dplyr::select(year,month) %>%
    distinct()

  for (cvar in cn2keep){

    crast <- list()

    OP.file <- paste0(
      dirname(cfile),"/",
      basename(tools::file_path_sans_ext(cfile)),".",
      cvar,
      "_CA.tif")

    if (file.exists(OP.file) & !overwrite){
      next()
    }


    for (itime in seq(1,nrow(times))){

      print(paste0(" -",cmodel,"=",cscenario,"|",cvar,": ",itime,"/",nrow(times)))

      cyear <- times$year[itime]
      cmonth <- times$month[itime]

      ccdf <- cdf.selected %>%
        filter(year == cyear,
               month == cmonth)

      crast[[itime]] <- rast(rasterFromXYZ(ccdf %>%
                                             dplyr::select(lon,lat,!!cvar)))



    }

    all.rast <- rast(crast)
    time(all.rast) <- times %>%
      mutate(date = as.Date(paste0(year,"/",month,"/01"))) %>%
      pull(date)

    writeRaster(all.rast,
                OP.file,
                overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

  }

}

# scp /home/femeunier/Documents/projects/Drying.CB/scripts/convert.CMIP6.to.rast.R hpc:/data/gent/vo/000/gvo00074/felicien/R

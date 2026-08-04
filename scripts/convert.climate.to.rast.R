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


files <- list.files("./outputs/",
                    pattern = "*.climate.RDS",
                    full.names = TRUE)
files <- (files[grepl("3IMERG|Berk|CAMS|chirps|chirpsv3|CRU|ERA5|GLDAS|GPCC|MSWEP|NCEP|GSMaP|CPC|JRA3Q",
                     files)])
# files <- files[grepl("Berk|CAMS|chirps|chirpsv3|CRU|ERA5|GLDAS|GPCC|MSWEP|NCEP|GSMaP|CPC|JRA3Q",
#                      files)]
files <- files[!grepl("CHELSA|CRUJRA|MERRA2|W5E5|GSWP",
                     files)]

overwrite = TRUE

for (ifile in seq(1,length(files))){

  cfile <- files[ifile]
  cdataset <- strsplit(basename(cfile),
                       "\\.")[[1]][2]
  print(paste0("-",cdataset))

  OP.file <- paste0(tools::file_path_sans_ext(cfile),".tif")

  if (file.exists(OP.file) & !overwrite){
    next()
  }

  ctemp <- readRDS(files[ifile]) %>%
    ungroup()

  if ("MAP" %in% colnames(ctemp)) {
    cdf <- ctemp %>%
      rename(pre = MAP)
  } else {
    cdf <- ctemp
  }

  cn <- colnames(cdf)

  cn2keep <- cn[!(cn %in% c("year","month","lon","lat"))]
  print(cn2keep)

  times <- cdf %>%
    ungroup() %>%
    dplyr::select(year,month) %>%
    distinct()

  for (cvar in cn2keep){

    crast <- list()
    for (itime in seq(1,nrow(times))){

      print(paste0(" -",cdataset,"|",cvar,": ",itime,"/",nrow(times)))

      cyear <- times$year[itime]
      cmonth <- times$month[itime]

      ccdf <- cdf %>%
        filter(year == cyear,
               month == cmonth) %>%
        mutate(lat = round(lat,digits = 2),
               lon = round(lon,digits = 2))

      crast[[itime]] <- rast(rasterFromXYZ(ccdf %>%
                                        dplyr::select(lon,lat,!!cvar)))



    }

    all.rast <- rast(crast)
    time(all.rast) <- times %>%
      mutate(date = as.Date(paste0(year,"/",month,"/01"))) %>%
      pull(date)

    writeRaster(all.rast,
                paste0("./outputs/all.climate/",cdataset,"_",cvar,"_all.years.tif"),
                overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

  }

}

# scp /Users/felicien/Documents/projects/Drying.CB/scripts/convert.climate.to.rast.R hpc:/data/gent/vo/000/gvo00074/felicien/R

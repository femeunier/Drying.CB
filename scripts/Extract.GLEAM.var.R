rm(list = ls())

library(dplyr)
library(lubridate)
library(terra)
library(sf)
library(Drying.CB)

files <- list.files("/data/gent/vo/000/gvo00074/felicien/GLEAM/",
                    pattern = "*.tif$",
                    full.names = TRUE)

Mask <- read_sf("./data/Rainforests.shp")

df.all <- data.frame()

baseline_start <- as.Date("1961-01-01")
baseline_end   <- as.Date("2014-12-31")

for (ifile in seq(1,length(files))){

  cfile <- files[ifile]
  file.split <- strsplit(basename(cfile),"_")[[1]]
  cproduct <- "GLEAM"
  cvar <- file.split[1]

  if (!(cvar %in% c("SMrz","SMs","Ep","S"))) next()

  print(paste0(cproduct,"-",cvar))

  cdata <- rast(cfile)

  dates <- as.Date(paste0(substr(names(cdata),1,4),"/",
                          substr(names(cdata),6,7),"/01"))

  cdata.msk <- crop(mask(cdata,Mask),
                    ext(-25,65,-25,25))

  ts <- global(cdata.msk,mean,na.rm = TRUE)

  cdf2save <- data.frame(time = dates,
                         mean = ts$mean) %>%
    mutate(year = year(time),
           month = month(time)) %>%
    rename(value = mean) %>%
    mutate(var = cvar,
           product = cproduct)

  if ((cvar %in% c("tas","tasmin","tasmax")) &
      (mean(cdf2save$value,na.rm = TRUE) > 200)){
    cdf2save <- cdf2save %>%
      mutate(value = value - 273.15)
  }

  df.all <- bind_rows(df.all,
                      cdf2save)

  ##########################################################
  # Anomalies
  anomalies <- anomalies_spatraster(input = cdata.msk,
                                    baseline_start = baseline_start,
                                    baseline_end   = baseline_end,
                                    detrend = TRUE)

  writeRaster(anomalies$trend,
              paste0("./outputs/CMIP6.CA/",
                     "GLEAM_",cvar,"_trends.tif"),
              overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

  time(anomalies$anom) <- as.Date(dates)
  writeRaster(anomalies$anom,
              paste0("./outputs/CMIP6.CA/",
                     "GLEAM_",cvar,"_anomalies.tif"),
              overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))



}

saveRDS(df.all,
        "./outputs/All.GLEAM.CA.RDS")

# scp /home/femeunier/Documents/projects/Drying.CB/scripts/Extract.GLEAM.var.R hpc:/data/gent/vo/000/gvo00074/felicien/R/


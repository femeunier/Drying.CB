rm(list = ls())

library(dplyr)
library(lubridate)
library(terra)
library(sf)
library(Drying.CB)

vars <- c("NDVI",
          "EVI",
          "Pixel.Reliability")

files <- list.files("/kyukon/data/gent/vo/000/gvo00074/felicien/R/outputs/",
                    pattern = paste0(paste0(vars,".*all.years.tif$"),collapse = "|"),
                    full.names = TRUE)

Mask <- read_sf("./data/Rainforests.shp")

df.all <- data.frame()

baseline_start <- as.Date("1961-01-01")
baseline_end   <- as.Date("2014-12-31")

for (ifile in seq(1,length(files))){

  cfile <- files[ifile]
  file.split <- strsplit(basename(cfile),"_")[[1]]
  cproduct <- "MODIS"
  cvar <- file.split[1]

  print(paste0(cproduct,"-",cvar))

  cdata <- rast(cfile)

  dates <- time(cdata)

  cdata.msk <- crop(mask(cdata,Mask),
                    ext(-25,65,-25,25))
  cdata.msk <- project(cdata.msk,
                       "EPSG:4326")

  ts <- global(cdata.msk,mean,na.rm = TRUE)

  cdf2save <- data.frame(time = dates,
                         mean = ts$mean) %>%
    mutate(year = year(time),
           month = month(time)) %>%
    rename(value = mean) %>%
    mutate(var = cvar,
           product = cproduct)

  df.all <- bind_rows(df.all,
                      cdf2save)


  ##########################################################
  # Anomalies
  anomalies <- anomalies_spatraster(input = cdata.msk,
                                    baseline_start = baseline_start,
                                    baseline_end   = baseline_end,
                                    detrend = TRUE)

  writeRaster(anomalies$trend,
              paste0("./outputs/",
                     "MODIS_",cvar,"_trends.tif"),
              overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

  time(anomalies$anom) <- as.Date(dates)
  writeRaster(anomalies$anom,
              paste0("./outputs/",
                     "MODIS_",cvar,"_anomalies.tif"),
              overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))


}

saveRDS(df.all,
        "./outputs/All.MODIS.CA.RDS")

# scp /home/femeunier/Documents/projects/Drying.CB/scripts/Extract.MODIS.var.R hpc:/data/gent/vo/000/gvo00074/felicien/R/


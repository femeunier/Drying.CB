rm(list = ls())

library(dplyr)
library(lubridate)
library(terra)
library(sf)
library(Drying.CB)

files <- list.files("./data/GRACE/",
                    pattern = "*.tif$",
                    full.names = TRUE)

Mask <- read_sf("./data/Rainforests.shp")

baseline_start <- as.Date("2000-01-01")
baseline_end   <- as.Date("2024-12-31")

df.all <- data.frame()

for (ifile in seq(1,length(files))){

  cfile <- files[ifile]
  file.split <- strsplit(basename(cfile),"_")[[1]]
  cproduct <- "GRACE"
  cvar <- file.split[1]

  if (!(cvar %in% c("tws","std","model","leakage"))) next()

  print(paste0(cproduct,"-",cvar))

  cdata <- rast(cfile)

  dates <- time(cdata)

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

  df.all <- bind_rows(df.all,
                      cdf2save)

  ##########################################################
  # Anomalies
  anomalies <- anomalies_spatraster_roll(input = cdata.msk,
                                         baseline_start = baseline_start,
                                         baseline_end   = baseline_end,
                                         detrend = FALSE)

  writeRaster(anomalies$trend,
              paste0("/data/gent/vo/000/gvo00074/felicien/R/outputs/Drying.CB/",
                     "GRACE_",cvar,"_trends.tif"),
              overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))


  time(anomalies$anom) <- as.Date(dates)
  writeRaster(anomalies$anom,
              paste0("/data/gent/vo/000/gvo00074/felicien/R/outputs/Drying.CB/",
                     "GRACE_",cvar,"_anomalies.tif"),
              overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))


  time(anomalies$z_anom) <- as.Date(dates)
  writeRaster(anomalies$z_anom,
              paste0("/data/gent/vo/000/gvo00074/felicien/R/outputs/Drying.CB/",
                     "GRACE_",cvar,"_Zanomalies.tif"),
              overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

  time(anomalies$roll_mean_input) <- as.Date(anomalies$roll_times)
  writeRaster(anomalies$roll_mean_input,
              paste0("/data/gent/vo/000/gvo00074/felicien/R/outputs/Drying.CB/",
                     "GRACE_",cvar,"_Rollmeaninput.tif"),
              overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

  writeRaster(anomalies$trend_z_anom,
              paste0("/data/gent/vo/000/gvo00074/felicien/R/outputs/Drying.CB/",
                     "GRACE_",cvar,"_trendsZanomalies.tif"),
              overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

  writeRaster(anomalies$trend_anom,
              paste0("/data/gent/vo/000/gvo00074/felicien/R/outputs/Drying.CB/",
                     "GRACE_",cvar,"_trendsanomalies.tif"),
              overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))


}

saveRDS(df.all,
        "./outputs/All.GRACE.CA.RDS")

# scp /home/femeunier/Documents/projects/Drying.CB/scripts/Extract.GRACE.var.R hpc:/data/gent/vo/000/gvo00074/felicien/R/


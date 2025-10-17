rm(list = ls())

library(dplyr)
library(lubridate)
library(terra)
library(sf)
library(Drying.CB)

products <- c("GOSIF","FLUXSAT","MODIS")

dirs <- c("GOSIF.GPP","FluxSat","MODIS_GPP") # "Zheng"

main.dir <- '/data/gent/vo/000/gvo00074/felicien/GPP_data'

Mask <- read_sf("./data/Rainforests.shp")

df.all <- data.frame()

baseline_start <- as.Date("2000-01-01")
baseline_end   <- as.Date("2024-12-31")

for (iproduct in seq(1,length(products))){

  cproduct <- products[iproduct]
  cfile <- file.path(main.dir,dirs[iproduct],
                     paste0("gpp.",cproduct,".tif"))

  if (!file.exists(cfile)) next()

  print(cproduct)

  cdata <- rast(cfile)
  dates <- paste0(substr(names(cdata),1,4),"/",
                  substr(names(cdata),6,7),"/01")
  time(cdata) <- as.Date(dates)
  cdata.msk <- crop(mask(cdata,Mask),
                    ext(-25,65,-25,25))*10

  ts <- global(cdata.msk,mean,na.rm = TRUE)

  cdf2save <- data.frame(time = time(cdata),
                         mean = ts$mean) %>%
    mutate(year = year(time),
           month = month(time)) %>%
    rename(value = mean) %>%
    mutate(var = "GPP",
           product = cproduct)

  df.all <- bind_rows(df.all,
                      cdf2save)

  ##########################################################
  # Anomalies
  anomalies <- anomalies_spatraster(input = cdata.msk,
                                    baseline_start = baseline_start,
                                    baseline_end   = baseline_end,
                                    detrend = FALSE)

  writeRaster(anomalies$trend,
              paste0("/data/gent/vo/000/gvo00074/felicien/R/outputs/Drying.CB/",
                     cproduct,"_GPP_trends.tif"),
              overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))


  time(anomalies$anom) <- as.Date(dates)
  writeRaster(anomalies$anom,
              paste0("/data/gent/vo/000/gvo00074/felicien/R/outputs/Drying.CB/",
                     cproduct,"_GPP_anomalies.tif"),
              overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))


  time(anomalies$z_anom) <- as.Date(dates)
  writeRaster(anomalies$z_anom,
              paste0("/data/gent/vo/000/gvo00074/felicien/R/outputs/Drying.CB/",
                     cproduct,"_GPP_Zanomalies.tif"),
              overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

  time(anomalies$roll_mean_input) <- as.Date(anomalies$roll_times)
  writeRaster(anomalies$roll_mean_input,
              paste0("/data/gent/vo/000/gvo00074/felicien/R/outputs/Drying.CB/",
                     cproduct,"_GPP_Rollmeaninput.tif"),
              overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

  writeRaster(anomalies$trend_z_anom,
              paste0("/data/gent/vo/000/gvo00074/felicien/R/outputs/Drying.CB/",
                     cproduct,"_GPP_trendsZanomalies.tif"),
              overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

  writeRaster(anomalies$trend_anom,
              paste0("/data/gent/vo/000/gvo00074/felicien/R/outputs/Drying.CB/",
                     cproduct,"_GPP_trendsanomalies.tif"),
              overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

}

saveRDS(df.all,
        "./outputs/All.GPP.CA.RDS")

# scp /home/femeunier/Documents/projects/Drying.CB/scripts/Extract.GPP.var.R hpc:/kyukon/data/gent/vo/000/gvo00074/felicien/R/



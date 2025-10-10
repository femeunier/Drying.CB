rm(list = ls())

library(dplyr)
library(lubridate)
library(terra)
library(sf)
library(Drying.CB)

files <- list.files("./data/GRACE/",
                    pattern = "*.tif$",
                    full.names = TRUE)[4]

Mask <- read_sf("./data/shapefiles/Rainforests.shp")

baseline_start <- as.Date("1981-01-01")
baseline_end   <- as.Date("2010-12-31")

df.all <- data.frame()

for (ifile in seq(1,length(files))){

  cfile <- files[ifile]
  file.split <- strsplit(basename(cfile),"_")[[1]]
  cproduct <- "GRACE"
  cvar <- file.split[1]

  if (!(cvar %in% c("tws","std","model","leakage"))) next()

  print(paste0(cproduct,"-",cvar))

  cdata <- rast(cfile)

  dates <- t <- time(cdata)

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
  anomalies <- anomalies_spatraster(input,
                                    baseline_start = baseline_start,
                                    baseline_end   = baseline_end,
                                    detrend = FALSE)

  writeRaster(z_anom,
              paste0("./outputs/",
                     cvar,"_zanomalies.tif"),
              overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))


}

plot(z_anom[[which(year(time(z_anom)) == 2024)]])

saveRDS(df.all,
        "./outputs/All.GRACE.CA.RDS")

# scp /home/femeunier/Documents/projects/Drying.CB/scripts/Extract.GRACE.var.R hpc:/data/gent/vo/000/gvo00074/felicien/R/


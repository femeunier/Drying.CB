rm(list = ls())

library(dplyr)
library(lubridate)
library(terra)
library(sf)
library(Drying.CB)

files <- list.files("/data/gent/vo/000/gvo00074/ED_common_data/met/Precip.Tropics/CMIP6/",
                    pattern = "*.tif$",
                    full.names = TRUE)

Mask <- read_sf("./data/Rainforests.shp")

df.all <- data.frame()

for (ifile in seq(1,length(files))){

  cfile <- files[ifile]
  file.split <- strsplit(basename(cfile),"\\.")[[1]]
  cmodel <- file.split[3]
  cscenario <- file.split[4]
  cvar <- strsplit(file.split[5], "_")[[1]][1]

  if (cscenario == "historical"){
    baseline_start <- as.Date("2000-01-01")
    baseline_end   <- as.Date("2014-12-31")
  } else {
    baseline_start <- as.Date("2000-01-01")
    baseline_end   <- as.Date("2024-12-31")
  }

  print(paste0(cmodel,"-",cscenario,"-",cvar,"-",ifile,"/",length(files)))

  cdata <- rast(cfile)

  dates <- time(cdata)

  cdata.msk <- crop(mask(cdata,Mask),
                    ext(-25,65,-25,25))


  if (cvar == "pr"){
    cdata.msk <- cdata.msk*86400*365.25/12
  }

  ts <- global(cdata.msk,mean,na.rm = TRUE)

  cdf2save <- data.frame(time = dates,
                         mean = ts$mean) %>%
    mutate(year = year(time),
           month = month(time)) %>%
    rename(value = mean) %>%
    mutate(var = cvar,
           model = cmodel,
           scenario = cscenario)

  if (cvar == "pr"){

    cdf2save <- cdf2save %>%
      mutate(var = "pre")
    cvar <- "pre"

  }

  df.all <- bind_rows(df.all,
                      cdf2save)

  ##########################################################
  # Anomalies
  anomalies <- anomalies_spatraster(input = cdata.msk,
                                    baseline_start = baseline_start,
                                    baseline_end   = baseline_end,
                                    detrend = FALSE)

  writeRaster(anomalies$trend,
              paste0("/data/gent/vo/000/gvo00074/felicien/R/outputs/Drying.CB/CMIP6/",
                     cmodel,"_",cscenario,"_",cvar,"_trends.tif"),
              overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

  time(anomalies$anom) <- as.Date(dates)
  writeRaster(anomalies$anom,
              paste0("/data/gent/vo/000/gvo00074/felicien/R/outputs/Drying.CB/CMIP6/",
                     cmodel,"_",cscenario,"_",cvar,"_anomalies.tif"),
              overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))


  time(anomalies$z_anom) <- as.Date(dates)
  writeRaster(anomalies$z_anom,
              paste0("/data/gent/vo/000/gvo00074/felicien/R/outputs/Drying.CB/CMIP6/",
                     cmodel,"_",cscenario,"_",cvar,"_Zanomalies.tif"),
              overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))


}

saveRDS(df.all,
        "./outputs/All.CMIP6.vars.CA.RDS")

# scp /home/femeunier/Documents/projects/Drying.CB/scripts/Extract.CMIP6.climate.var.R hpc:/data/gent/vo/000/gvo00074/felicien/R/


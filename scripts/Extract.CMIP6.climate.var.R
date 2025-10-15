rm(list = ls())

library(dplyr)
library(lubridate)
library(terra)
library(sf)
library(Drying.CB)

files <- list.files("/data/gent/vo/000/gvo00074/ED_common_data/met/Precip.Tropics/CMIP6/",
                    pattern = "*.tif$",
                    full.names = TRUE)
files <- files[grepl("E3SM-2-0",files) &
                 grepl("historical",files) &
                 grepl("pr",files)]

Mask <- read_sf("./data/Rainforests.shp")

df.all <- data.frame()

baseline_start <- as.Date("1961-01-01")
baseline_end   <- as.Date("2014-12-31")

for (ifile in seq(1,length(files))){

  cfile <- files[ifile]
  file.split <- strsplit(basename(cfile),"\\.")[[1]]
  cmodel <- file.split[3]
  cscenario <- file.split[4]
  cvar <- strsplit(file.split[5], "_")[[1]][1]

  if (cscenario == "historical"){
    baseline_start <- as.Date("1961-01-01")
    baseline_end   <- as.Date("2014-12-31")
  } else {
    baseline_start <- as.Date("2015-01-01")
    baseline_end   <- as.Date("2099-12-31")
  }

  print(paste0(cmodel,"-",cscenario,"-",cvar))

  cdata <- rast(cfile)

  dates <- time(cdata)

  cdata.msk <- crop(mask(cdata,Mask),
                    ext(-25,65,-25,25))

  pt <- vect(data.frame(x = -13.75, y = 9.75),
             geom = c("x","y"), crs = "EPSG:4326")
  pt_r <- project(pt, crs(cdata.msk))
  vals <- extract(cdata.msk, pt_r, ID = FALSE)*86400*365.25/12

  summary(lm(data = data.frame(time = dates,
                       value = as.numeric(as.vector(vals))) %>%
               filter(time >= baseline_start,
                      time <= baseline_end),
     formula = value ~ time))

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

  }

  df.all <- bind_rows(df.all,
                      cdf2save)


  input = cdata.msk[[1:10]]

  for (i in seq(1,nlyr(input))){
    print(i)
    input[[i]] <- input[[i]]*0 + i
  }

  stop()
  ##########################################################
  # Anomalies
  anomalies <- anomalies_spatraster(input = cdata.msk,
                                    baseline_start = baseline_start,
                                    baseline_end   = baseline_end,
                                    detrend = FALSE)

  writeRaster(anomalies$trend,
              paste0("./outputs/CMIP6.CA/",
                     cmodel,"_",cscenario,"_",cvar,"_trends.tif"),
              overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

  time(anomalies$anom) <- as.Date(dates)
  writeRaster(anomalies$anom,
              paste0("./outputs/CMIP6.CA/",
                     cmodel,"_",cscenario,"_",cvar,"_anomalies.tif"),
              overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

}

saveRDS(df.all,
        "./outputs/All.CMIP6.vars.CA.RDS")

# scp /home/femeunier/Documents/projects/Drying.CB/scripts/Extract.CMIP6.climate.var.R hpc:/data/gent/vo/000/gvo00074/felicien/R/


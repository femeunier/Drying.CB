rm(list = ls())

library(dplyr)
library(lubridate)
library(terra)
library(sf)

files <- list.files("/data/gent/vo/000/gvo00074/felicien/R/outputs/all.climate",
                    pattern = "*.tif$",
                    full.names = TRUE)
files <- files[!grepl("CRUJRA_",files)]

Mask <- read_sf("./data/Rainforests.shp")

df.all <- data.frame()

for (ifile in seq(1,length(files))){

  cfile <- files[ifile]
  file.split <- strsplit(basename(cfile),"_")[[1]]
  cproduct <- file.split[1]
  cvar <- file.split[2]

  if (!(cvar %in% c("tas","pre","tasmin","tasmax"))) next()

  print(paste0(cproduct,"-",cvar))

  cdata <- rast(cfile)
  cdata.msk <- crop(mask(cdata,Mask),
                    ext(-25,65,-25,25))

  ts <- global(cdata.msk,mean,na.rm = TRUE)

  cdf2save <- data.frame(time = time(cdata),
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
}

saveRDS(df.all,
        "./outputs/All.climatevars.CA.RDS")

# scp /home/femeunier/Documents/projects/Drying.CB/scripts/Extract.allclimate.vars.R hpc:/data/gent/vo/000/gvo00074/felicien/R/


rm(list = ls())

library(dplyr)
library(lubridate)
library(terra)
library(sf)

products <- c("GOSIF","FLUXSAT","MODIS")

dirs <- c("GOSIF.GPP","FluxSat","MODIS_GPP") # "Zheng"

main.dir <- '/data/gent/vo/000/gvo00074/felicien/GPP_data'

Mask <- read_sf("./data/Rainforests.shp")

df.all <- data.frame()
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
                    ext(-25,65,-25,25))

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


}

saveRDS(df.all,
        "./outputs/All.GPP.CA.RDS")

# scp /home/femeunier/Documents/projects/Drying.CB/scripts/Extract.GPP.R hpc:/kyukon/data/gent/vo/000/gvo00074/felicien/R/



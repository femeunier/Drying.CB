rm(list = ls())

library(dplyr)
library(ncdf4)
library(reshape2)
library(lubridate)
library(matlab)

init.years <- seq(2000,2024,1)

WD <- "/data/gent/vo/000/gvo00074/ED_common_data/met/GLDAS"

df.CRU <- data.frame()


for (iyear in seq(1,length(init.years))){

  print(iyear/length(init.years))

  for (month in seq(1,12)){

    cfile <- file.path(WD,paste0("GLDAS_NOAH025_M.A",
                                 init.years[iyear],sprintf("%02d",month),
                                 ".021.nc4"))

    if (file.exists(cfile)){


      nc <- nc_open(cfile)

      prate <- ncvar_get(nc,"Rainf_tavg")
      lats <- (as.vector(ncvar_get(nc,"lat")))
      lons <- ncvar_get(nc,"lon")
      times <- as.Date(ncvar_get(nc,"time"),
                      origin = "2000-01-01")

      nc_close(nc)

      N <- lubridate::days_in_month(times)

      prate.df <- melt(prate) %>%
        mutate(lon = (lons)[Var1],
               lat = (lats)[Var2]) %>%
        filter(!is.na(value)) %>%
        filter(lon >= -180 & lon <= 180,
               lat >= -30 & lat <= 30) %>%
        dplyr::select(lat,lon,value) %>%
        mutate(prate = value*N*86400) %>%
        dplyr::select(-value)

      prate.df.time <- prate.df %>%
        mutate(date = times) %>%
        mutate(month = month(date),
               year = year(date))

      cdf <- prate.df.time %>%
        group_by(year,month,lat,lon) %>%
        summarise(MAP = sum(prate),
                  .groups = "keep")

    }

    df.CRU <- bind_rows(list(df.CRU,
                             cdf))

  }
}

saveRDS(df.CRU,"./outputs/df.GLDAS.Tropics.RDS")

# scp /home/femeunier/Documents/projects/YGB/scripts/summarise.precip.GLDAS.R hpc:/kyukon/data/gent/vo/000/gvo00074/felicien/R



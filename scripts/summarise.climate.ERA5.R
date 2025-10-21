rm(list = ls())

library(dplyr)
library(ncdf4)
library(reshape2)
library(lubridate)

init.years <- seq(1940,2024,1)
vars <- c("tp","t2m","d2m")

WD <- "/data/gent/vo/000/gvo00074/ED_common_data/met/global"

df.ERA5 <- data.frame()
for (iyear in seq(1,length(init.years))){

  print(iyear/length(init.years))

  cfile <- file.path(WD,paste0("ERA5_global_",init.years[iyear],".nc"))
  if (file.exists(cfile)){

    nc <- nc_open(cfile)
    lats <- ncvar_get(nc,"latitude")
    lons <- ncvar_get(nc,"longitude")

    cyear <- init.years[iyear]

    if (cyear < 2024){
        times <- as.Date(ncvar_get(nc,"time")/24,
                     origin = "1900-01-01")
        months <- month(times)
    } else {
        times <- as.Date(ncvar_get(nc,"valid_time")/86400,
                     origin = "1970-01-01")
    }

    for (cvar in vars){
      prate <- ncvar_get(nc,cvar)

      prate.df.time <- melt(prate) %>%
        ungroup() %>%
        mutate(lon = (lons)[Var1],
               lat = (lats)[Var2],
               time = times[Var3]) %>%
        dplyr::select(lat,lon,time,value) %>%
        rename(!!cvar := value) %>%
        filter(abs(lat) <= 30) %>%
        mutate(month = month(time),
               year = year(time),
               day = day(time))

      if (cvar == "tp"){
        cdf <- prate.df.time %>%
          group_by(year,month,lat,lon) %>%
          summarise(pre = sum(tp)*1000*3,
                    .groups = "keep")
      } else if (cvar == "d2m") {

        tmp <- ncvar_get(nc,"t2m")

        tmp.df.time <- melt(tmp) %>%
          ungroup() %>%
          mutate(lon = (lons)[Var1],
                 lat = (lats)[Var2],
                 time = times[Var3]) %>%
          dplyr::select(lat,lon,time,value) %>%
          rename(t2m = value) %>%
          filter(abs(lat) <= 30) %>%
          mutate(month = month(time),
                 year = year(time),
                 day = day(time))

        cdf.all <- prate.df.time %>%
          left_join(tmp.df.time,
                    by = c("lon","lat","time","year","month","day")) %>%
          mutate(vpd = max(0,
                           0.6108*exp(17.27*t2m/(t2m + 237.3)) -
                             0.6108*exp(17.27*d2m/(d2m + 237.3))))


        cdf <- cdf %>%
          left_join(cdf.all %>%
                      group_by(year,month,lat,lon) %>%
                      summarise(vpd = mean(vpd,na.rm = TRUE),
                                .groups = "keep"),
                    c("year","month","lon","lat"))



      } else {

        ccdf <- data.frame()

        for (imonth in seq(1,12)){

          print(paste0("- ",imonth))

          ccdf <- bind_rows(ccdf,
                            prate.df.time %>%
                              filter(month == imonth) %>%
            group_by(year,month,day,lat,lon) %>%
            summarise(tmin = min(t2m),
                      tmax = max(t2m),
                      tmean = mean(t2m),
                      .groups = "keep") %>%
            group_by(year,month,lat,lon) %>%
            summarise(tmin = mean(tmin),
                      tmax = mean(tmax),
                      tmp = mean(tmean),
                      .groups = "keep"))
        }

        cdf <- cdf %>%
          left_join(ccdf,
          by = c("year","month","lon","lat"))
      }
    }

    nc_close(nc)

  }

  df.ERA5 <- bind_rows(df.ERA5,
                       cdf)
}

saveRDS(df.ERA5,"./outputs/df.ERA5.Tropics.climate.RDS")

# scp /home/femeunier/Documents/projects/Drying.CB/scripts/summarise.climate.ERA5.R hpc:/kyukon/data/gent/vo/000/gvo00074/felicien/R


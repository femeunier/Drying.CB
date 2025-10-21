rm(list = ls())

library(dplyr)
library(ncdf4)
library(reshape2)
library(plotbiomes)
library(raster)
library(ggplot2)
library(sf)
library(PEcAn.data.atmosphere)
library(lubridate)
library(Drying.CB)

years <- 1901:2024

vars <- c("pre","tmp","tmin","tmax","spfh")
# vars <- c("pre")

# dir <- "/home/femeunier/Documents/projects/TrENDY.analyses/data"
dir <- "/data/gent/vo/000/gvo00074/felicien/TrENDY/inputs/CRUJRA3Q"

days2months <- c(31,28,31,30,31,30,
                 31,31,30,31,30,31)

all.days <- rep(1:12,days2months)
all.days.hours <- sort(rep(all.days,4))

df.all <- df.all.monthly <- data.frame()

prefix <- rev(c("3","3.5","3.02"))


origin = as.Date("1901-01-01 00:00:00")

for (cyear in years){

  print(paste0(cyear))
  cyr.df <- data.frame()

  temp.array <- array(data = NA,
                      dim = c(720,360,1460,length(vars)))

  ivar = 1
  for (cvar in vars){

    print(paste0("- ",cvar))

    i = 1 ; zip.file.exist = FALSE
    while (i <= length(prefix) & !zip.file.exist){

      zip.file <- file.path(dir,
                            paste0("crujra.v",prefix[i],".5d.",cvar,".",cyear,".365d.noc.nc.gz"))

      zip.file.exist <- file.exists(zip.file)

      i = i + 1
    }

    if (!zip.file.exist){
      warning(paste0("Zip file not found:", cvar," - ",cyear))
      next()
    }

    system2("gunzip",
            paste("-k",zip.file))

    nc.file <- file.path(dir,
                         paste0("crujra.v",prefix[i-1],".5d.",cvar,".",cyear,".365d.noc.nc"))

    nc <- nc_open(nc.file)
    lats <- ncvar_get(nc,"lat")
    lons <- ncvar_get(nc,"lon")

    data <- ncvar_get(nc,cvar)

    ivar <- ivar + 1

   days_since_origin <- ncvar_get(nc,"time")

  year <- year(origin) + floor(days_since_origin[1] / 365)
  days_remaining <- days_since_origin %% 365
  months <- month(as.Date(days_remaining,"1901/01/01"))
  days <- day(as.Date(days_remaining,"1901/01/01"))
  dates <- as.Date(paste0(year,"/",months,"/",days))

  # Construct result date by adding years and days manually


   cdf <- melt(data) %>%
      filter(!is.na(value)) %>% ungroup() %>%
      rename(lon = Var1,
             lat = Var2,
             time = Var3) %>%
      mutate(lon = lons[lon],
             lat = lats[lat],
             time = dates[time]) %>%
      mutate(month = month(time),
             day = day(time))


   if (cvar %in% c("tmin")){

     cdf <- cdf %>%
        group_by(lon,lat,month,day) %>%
        summarise(value = min(value,na.rm = TRUE),.groups = "keep")

   } else if (cvar %in% c("tmax")){

     cdf <- cdf %>%
        group_by(lon,lat,month,day) %>%
        summarise(value = max(value,na.rm = TRUE),.groups = "keep")


   }


   if (cvar == "tmp"){
     df2keep <- cdf
   }

   if (cvar == "spfh"){
     cdf <- cdf %>%
       left_join(df2keep %>%
                   rename(spfh = value),
                 by = c("lon","lat","month","day","time")) %>%
       mutate(vpd = vpd_from_T_q(value,spfh,101.325)) %>%
       ungroup() %>%
       mutate(value = vpd)
       dplyr::select(-c("spfh","vpd"))

       cdf.monthly <- cdf %>%
         group_by(lon,lat,month) %>%
         summarise(value = mean(value,na.rm = TRUE),.groups = "keep") %>%
         ungroup() %>%
         rename(vpd = value) %>%
         filter(abs(lat) < 30)

   } else {
     cdf.monthly <- cdf %>%
       group_by(lon,lat,month) %>%
       summarise(value = mean(value,na.rm = TRUE),.groups = "keep") %>%
       ungroup() %>%
       rename(!!cvar := value) %>%
       filter(abs(lat) < 30)
   }


    if (cvar == vars[1]){
      # cyr.df <- cdf
      cmonth.df <- cdf.monthly
    } else {
      # cyr.df <- cyr.df %>%
      #   left_join(cdf,
      #             by = c("lat","lon"))

      cmonth.df <- cmonth.df %>%
        left_join(cdf.monthly,
                  by = c("lat","lon","month"))
    }



    system2("rm",nc.file)

  }


  df.all.monthly <- bind_rows(list(df.all.monthly,
                                   cmonth.df %>%
                                     mutate(year = cyear)))

}

saveRDS(df.all.monthly %>%
            mutate(Ndays = as.numeric(lubridate::days_in_month(paste0(year,"/",month,"/01")))) %>%
            mutate(pre = pre*4*Ndays) %>%
            dplyr::select(-Ndays),
        paste0("./outputs/df.CRUJRA3Q.Tropics.climate.RDS"))

# scp /home/femeunier/Documents/projects/Drying.CB/scripts/summarise.climate.CRUJRA3Q.R hpc:/data/gent/vo/000/gvo00074/felicien/R/

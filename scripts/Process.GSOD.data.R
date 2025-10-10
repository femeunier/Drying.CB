rm(list = ls())

library(readr)
library(dplyr)
library(fs)
library(readr)
library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(dplyr)
library(readr)
library(fs)
library(stringr)
library(GSODTools)
library(GSODR)

selected.stations <- gsodstations %>%
  filter(LAT <= 25, LAT >= -25,
         LON >= -25, LON <= 65)

all.data <- data.frame()
for (irow in seq(1,nrow(selected.stations))){

  print(irow/nrow(selected.stations))
  year_init <- year(selected.stations$BEGIN)[irow]
  year_end <- year(selected.stations$END)[irow]
  station_id <- paste0(selected.stations$USAF[irow],
                       "-",
                       selected.stations$WBAN[irow])

  gsod_data <- tryCatch(get_GSOD(years = year_init:year_end,
                                 station = station_id),
                        error = function(e) return(NULL))

  if (!is.null(gsod_data)){

    if (nrow(gsod_data) == 0) next

    cdata <- gsod_data %>%
      rename(STATION = STNID) %>%
      dplyr::select(any_of(c(
        "STATION","YEAR","MONTH","DAY",
        "LATITUDE","LONGITUDE","PRCP","DEWP","TEMP","MIN","MAX"
      ))) %>%
      # parse DATE robustly (adjust if your format is e.g. "%Y%m%d")
      rename(year = YEAR,
             month = MONTH,
             day = DAY) %>%
      # ensure numeric measurement columns
      mutate(across(any_of(c("PRCP","DEWP","TEMP","MIN","MAX")), ~ suppressWarnings(as.numeric(.)))) %>%
      pivot_longer(cols = any_of(c("PRCP","DEWP","TEMP","MIN","MAX")),
                   names_to = "variable",
                   values_to = "value") %>%
      mutate(value = case_when(value == 9999.9 ~ NA_real_,
                               TRUE ~ value)) %>%
      group_by(STATION, LATITUDE, LONGITUDE, year, month, variable) %>%
      summarise(
        value.m = if_else(variable[1] %in% c("PRCP","EVAP"),
                          sum(value, na.rm = TRUE),   # monthly totals
                          mean(value, na.rm = TRUE)), # monthly means
        N = sum(!is.na(value)),                        # number of daily obs used
        .groups = "drop"
      )

    data <- cdata %>%
      rename(station_id = STATION) %>%
      dplyr::select(-c(LATITUDE,LONGITUDE)) %>%
      mutate(Ndays = days_in_month(as.Date(paste0(year,'/',month,"/01")))) %>%
      mutate(flag = case_when(Ndays == N ~ "full",
                              N == 0 ~ "empty",
                              TRUE ~ "partial")) %>%
      filter(flag != 'empty') %>%
      filter((flag == "full")|
               (flag == "partial" &
                  variable %in% c("TEMP","MIN","MAX","DEWP") &
                  N >= 20))


    all.data <- bind_rows(all.data,
                          data)

    print(paste0("-",nrow(data)))

  }
}


all.data %>%
  filter(variable == "PRCP") %>%
  pull(value.m) %>%
  hist()

metadata <- all.data %>%
  dplyr::select(station_id) %>%
  distinct() %>%
  left_join(selected.stations %>%
              dplyr::select(USAF,WBAN,LAT,LON,`ELEV(M)`) %>%
              mutate(station_id = paste0(USAF,"-",WBAN)) %>%
              rename(lon = LON,
                     lat = LAT,
                     elevation = `ELEV(M)`) %>%
              dplyr::select(station_id,lon,lat,elevation),
            by = c("station_id"))

saveRDS(metadata,
        "./outputs/metadata.GSOD.RDS")


data.wide <- all.data %>%
  ungroup() %>%
  distinct() %>%
  dplyr::select(-c(N,flag,Ndays)) %>%
  pivot_wider(values_from = "value.m",
              names_from = "variable") %>%
  rename(tas = TEMP,
         tasmin = MIN,
         tasmax = MAX,
         dewpoint = DEWP,
         pre = PRCP)

saveRDS(data.wide,
        "./outputs/data.GSOD.RDS")

system2("rsync",
        c("hpc:/kyukon/data/gent/vo/000/gvo00074/felicien/R/outputs/metadata.GSOD.RDS",
          "./data/station_data/GSOD/"))

system2("rsync",
        c("hpc:/kyukon/data/gent/vo/000/gvo00074/felicien/R/outputs/data.GSOD.RDS",
          "./data/station_data/GSOD/"))

# scp /home/femeunier/Documents/projects/Drying.CB/scripts/Process.GSOD.data.R hpc:/data/gent/vo/000/gvo00074/felicien/R/


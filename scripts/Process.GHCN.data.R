rm(list = ls())

library(readr)
library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(rnoaa)

# # Step 1
metadata <- read_fwf(
  # "./data/station_data/GHCN/ghcnd-stations.txt",
  "/data/gent/vo/000/gvo00074/felicien/station/GHCN/ghcnd-stations.txt",
  col_positions = fwf_positions(
    start = c(1, 13, 22, 32, 39, 42, 73, 77, 81),
    end   = c(11, 20, 30, 37, 40, 71, 75, 79, 85),
    col_names = c("id", "lat", "lon", "elevation", "state", "name",
                  "gsn_flag", "hcn_crn_flag", "wmo")
  ),
  col_types = cols(
    id = col_character(),
    lat = col_double(),
    lon = col_double(),
    elevation = col_double(),
    state = col_character(),
    name = col_character(),
    gsn_flag = col_character(),
    hcn_crn_flag = col_character(),
    wmo = col_character()
  ),
  trim_ws = TRUE
) %>%
  mutate(
    name = str_squish(name),
    gsn = if_else(gsn_flag == "GSN", TRUE, FALSE, missing = FALSE),
    wmo = na_if(str_trim(wmo), ""),
    state = na_if(str_trim(state), "")
  ) %>%
  dplyr::select(id, lat, lon, elevation, name, state, gsn, hcn_crn_flag, wmo)

stations_selected <- metadata %>%
  filter(abs(lat) <= 25,
         lon >= -25, lon <= 65)

monitors <- stations_selected$id
all.data <- meteo_pull_monitors(monitors,
                                var = c("PRCP","EVAP","TAVG",
                                        "RHAV","TMAX", "TMIN")) %>%
  mutate(year = year(date),
         month = month(date),
         day = day(date)) %>%
  left_join(metadata %>%
              dplyr::select(id,lon,lat),
            by = "id")

all.data.sum <- all.data %>%
  rename(station_id = id) %>%
  dplyr::select(-any_of(c("lat","lon","date"))) %>%
  dplyr::select(any_of(c(
    "station_id","year","month",
    "prcp","tavg","tmin","tmax","evap","rhav"))) %>%
  pivot_longer(cols = any_of(c("prcp","tavg","tmin","tmax")),
               names_to = "variable",
               values_to = "value") %>%
  group_by(station_id, year, month, variable) %>%
  summarise(
    value.m = case_when(variable[1] %in% c("prcp","evap") ~ sum(value, na.rm = TRUE),
                        variable[1] %in% c("rhav") ~ mean(value, na.rm = TRUE),
                        TRUE ~  mean(value, na.rm = TRUE))/10,
    N = sum(!is.na(value)),
    .groups = "drop"
  )


metadata.selected <- stations_selected %>%
  rename(station_id = id) %>%
  dplyr::select(station_id,lon,lat,elevation) %>%
  distinct()

saveRDS(metadata.selected,
        "./outputs/metadata.GHCN.RDS")

data <- all.data.sum %>%
  mutate(Ndays = days_in_month(as.Date(paste0(year,'/',month,"/01")))) %>%
  mutate(flag = case_when(Ndays == N ~ "full",
                          N == 0 ~ "empty",
                          TRUE ~ "partial")) %>%
  filter(flag != 'empty') %>%
  filter((flag == "full")|
           (flag == "partial" &
              variable %in% c("tavg","tmin","tmax","evap","rhav") &
              N >= 20))

data.wide <- data %>%
  dplyr::select(-c(N,flag,Ndays)) %>%
  pivot_wider(values_from = "value.m",
              names_from = "variable") %>%
  rename(tas = tavg,
         tasmin = tmin,
         tasmax = tmax,
         pre = prcp)

summary(data.wide %>%
          filter(pre>0) %>%
          pull(pre))

saveRDS(data.wide,
        "./outputs/data.GHCN.RDS")


system2("rsync",
        c("hpc:/kyukon/data/gent/vo/000/gvo00074/felicien/R/outputs/metadata.GHCN.RDS",
          "./data/station_data/GHCN/"))

system2("rsync",
        c("hpc:/kyukon/data/gent/vo/000/gvo00074/felicien/R/outputs/data.GHCN.RDS",
          "./data/station_data/GHCN/"))


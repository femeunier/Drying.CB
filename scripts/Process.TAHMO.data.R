rm(list = ls())

library(sf)
library(terra)
library(dplyr)
library(purrr)
library(ggplot2)
library(rnaturalearth)
library(stringr)
library(lubridate)
library(tidyr)
library(stringr)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

stations <- read.csv("./data/station_data/TAHMO/stations_metadata.csv")

ggplot() +
  # scale_fill_gradient(name = "Tree cover (%)") +
  guides(fill = "none") +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  geom_point(data = stations,
             aes(x = longitude, y = latitude), color = "red", size = 1.2) +
  coord_sf(xlim = c(-20, 60), ylim = c(-1, 1)*25, expand = FALSE) +
  theme_bw()

saveRDS(stations %>%
          rename(station_id = station.code,
                 lat = latitude,
                 lon = longitude,
                 elevation = elevation..m.) %>%
          dplyr::select(station_id,lon,lat,elevation),
        "./data/station_data/TAHMO/metadata.TAHMO.RDS")

station.files <- list.files(path = "./data/station_data/TAHMO/",
                            pattern = "^TA.*.csv$",full.names = TRUE)

all.data <- data_frame()

for (ifile in seq(1,length(station.files))){

  cfile <- station.files[ifile]
  cstation <- tools::file_path_sans_ext(basename(cfile))
  cdata <- read.csv(cfile)

  print(paste0(cstation,' - ',ifile,'/',length(station.files)))

  cdata.long <- cdata %>%
    mutate(station = cstation,
           timestamp = as.Date(timestamp)) %>%
    mutate(year = year(timestamp),
           month = month(timestamp)) %>%
    pivot_longer(-c(timestamp,year,month,station),
                 names_to = "variable",
                 values_to = "value") %>%
    mutate(var = word(variable,1,sep = "\\.")) %>%
    mutate(sensor = word(variable,2,sep = "\\."),
           type = word(variable,3,sep = "\\."))

  all.data <- bind_rows(all.data,
                        cdata.long)

}

all.data <- all.data %>%
  mutate(type = case_when(sensor %in% c("AVG","MIN","MAX") ~ sensor,
                          type == "mm" ~ "",
                          TRUE ~ type),
         sensor = case_when(sensor %in% c("AVG","MIN","MAX") ~ "",
                          TRUE ~ sensor))

colnames(all.data)
unique(all.data$var)
unique(all.data$sensor)
unique(all.data$type)

all.data.sum <- all.data %>%
  group_by(station,year,month,variable,var,sensor,type) %>%
  summarise(value.m = case_when(var[1] == "precipitation" ~ sum(value,na.rm = TRUE),
                                var[1] == "temperature" ~ mean(value,na.rm = TRUE),
                                TRUE ~ FALSE),
            .groups = "keep")


data <- all.data.sum %>%
  ungroup() %>%
  mutate(type = case_when(type == "AVG" ~ "",
                          type == "" ~ "",
                          TRUE ~ tolower(type))) %>%
  mutate(var = case_when(var == "precipitation" ~ "pre",
                              var == "temperature" ~ "tas")) %>%
  mutate(variable = paste0(var,type)) %>%
  dplyr::select(-c(type,var,sensor)) %>%
  group_by(station,year,month,variable) %>%
  summarise(value.m = mean(value.m,na.rm = TRUE),
            .groups = "keep") %>%
  ungroup() %>%
  pivot_wider(names_from = variable,
              values_from = value.m)

saveRDS(data %>%
          rename(station_id = station),
        "./data/station_data/TAHMO/data.TAHMO.RDS")

ggplot(data = data %>%
         pivot_longer(cols = c(pre,tas,tasmax,tasmin),
                      names_to = "variable",
                      values_to = "value.m"),
      aes(x = year + (month - 1/2)/12,
          y = value.m,
          color = as.factor(station))) +
  geom_line() +
  facet_wrap(~interaction(variable),
             scales = "free") +
  theme_bw()

all.data.sum.anomaly <- all.data.sum %>%
  group_by(station,var,type,month) %>%
  mutate(m = mean(value.m,na.rm = TRUE)) %>%
  group_by(station,var,type) %>%
  mutate(anomaly = value.m - m)

stations2remove <- all.data.sum.anomaly %>%
  filter(var == "temperature" & anomaly > 2.5) %>%
  pull(station) %>%
  unique()


ggplot(data = all.data.sum.anomaly %>%
         filter(!(station %in% c(stations2remove))),
       aes(x = year + (month - 1/2)/12,
           y = anomaly,
           color = as.factor(station))) +
  geom_line() +
  facet_wrap(~interaction(var,type),
             scales = "free") +
  theme_bw()

all.data.sum.anomaly.m <- all.data.sum.anomaly %>%
  filter(!(station %in% c(stations2remove))) %>%
  group_by(var,type,year,month) %>%
  summarise(anomaly.m = mean(anomaly,na.rm = TRUE),
            .groups = "keep")

ggplot(data = all.data.sum.anomaly.m,
       aes(x = year + (month - 1/2)/12,
           y = anomaly.m)) +
  geom_line() +
  facet_wrap(~interaction(var,type),
             scales = "free") +
  theme_bw()


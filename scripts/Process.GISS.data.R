rm(list = ls())

library(dplyr)
library(ggplot2)
library(tidyr)
library(Drying.CB)


stations <- read_giss_inventory("./data/station_data/GISS/v4.temperature.inv.txt")

stations.selected <- stations %>%
  rename(elevation = elevation_m) %>%
  dplyr::select(station_id,lat,lon,elevation) %>%
  filter(abs(lat) <= 25,
         lon >= -25,lon <= 65)




df <- read_giss_v4("/home/femeunier/Documents/projects/Drying.CB/data/station_data/GISS/v4.mean_GISS_homogenized.txt",
                   scale = 0.01)
df.selected <- df %>%
  filter(station_id %in% (stations.selected %>%
                            pull(station_id)))

unique(df.selected$station_id)
hist(df.selected$latitude)
hist(df.selected$longitude)
hist(df.selected$value)


world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  guides(fill = "none") +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  geom_point(data = stations.selected %>%
               ungroup() %>%
               dplyr::select(lon,lat) %>%
               distinct(),
             aes(x = lon, y = lat),
             color = "red",
             size = 0.6, alpha = 0.5) +
  coord_sf(xlim = c(-30, 65), ylim = c(-25, 25), expand = FALSE) +
  theme_bw()

saveRDS(df.selected %>%
          ungroup() %>%
          dplyr::select(-c(longitude,latitude,name,elevation_m)) %>%
          rename(tas = value),
        "./data/station_data/GISS/data.GISS.RDS")

saveRDS(stations.selected %>%
          ungroup() %>%
          dplyr::select(station_id,lon,lat,elevation) %>%
          distinct(),
        "./data/station_data/GISS/metadata.GISS.RDS")



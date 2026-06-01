rm(list = ls())

library(dplyr)
library(ggplot2)
library(tidyr)
library(Dryind.CB)

files <- c("tmn.1708171135.clean.dtb",
           "tmp.1708041519.clean.dtb",
           "tmx.1708171215.clean.dtb",
           "pre.1704241136.clean.dtb")

vars <- c("tasmin","tas","tasmax","pre")

df.all <- data.frame()

for (ivar in seq(1,length(vars))){

  cvar <- vars[ivar]

  print(cvar)

  df <- read_cru_ts(
    paste0("/home/femeunier/Documents/projects/Drying.CB/data/station_data/CRU/",
           files[ivar]))

  df.selected <- df %>%
    filter(abs(lat) <= 25,
           lon >= -25,lon <= 65) %>%
    filter(!(lon == 0 & lat == 0)) %>%
    dplyr::select(wmo,lon,lat,elev_m,year,month,value) %>%
    rename(station_id = wmo,
           elevation = elev_m) %>%
    mutate(value = value/10) %>%
    na.omit()

  df.all <- bind_rows(df.all,
                      df.selected %>%
                        mutate(variable = cvar))

}

df.all.wide <- df.all %>%
  group_by(station_id,lon,lat,elevation,variable,year,month) %>%
  summarise(value = mean(value,na.rm = TRUE),
            .groups = "keep") %>%
  pivot_wider(names_from = "variable",
              values_from = "value")


hist(df.all.wide$lon)
hist(df.all.wide$lat)
hist(df.all.wide$elevation)

df.all.wide %>%
  filter(lon < -180)

sort(unique(df.all.wide$country))
sort(unique(df$station_name))



hist(df.selected$value)

hist(df.selected$year)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  guides(fill = "none") +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  geom_point(data = df.all.wide %>%
               ungroup() %>%
               dplyr::select(lon,lat) %>%
               distinct(),
             aes(x = lon, y = lat),
             color = "red",
             size = 0.6, alpha = 0.5) +
  coord_sf(xlim = c(-30, 65), ylim = c(-25, 25), expand = FALSE) +
  theme_bw()

saveRDS(df.all.wide %>%
          ungroup() %>%
          dplyr::select(-c(lon,lat)) %>%
          ungroup() ,
        "./data/station_data/CRU/data.CRU.RDS")

metadata2save <- df.all.wide %>%
  ungroup() %>%
  dplyr::select(station_id,lon,lat,elevation) %>%
  distinct()

saveRDS(metadata2save,
        "./data/station_data/CRU/metadata.CRU.RDS")


hist(metadata2save$elevation)

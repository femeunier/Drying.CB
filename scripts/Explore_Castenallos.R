rm(list = ls())

library(dplyr)
library(ggplot2)
library(ggthemes)

md <- readRDS("/Users/felicien/Documents/projects/Drying.CB/data/station_data/Castellanos/metadata.Castellanos.RDS")


long_df <- readRDS("/Users/felicien/Documents/projects/Drying.CB/data/station_data/Castellanos/data.Castellanos.RDS")

cdf <- long_df %>%
  filter(year >= 1985, year <= 2014)

N <- cdf %>%
  group_by(source,station_id) %>%
  summarise(N.years = length(unique(year)),
            years = paste0(min(year),"-",max(year))) %>%
  filter(N.years > 20)

N %>%
  pull(source) %>%
  unique()

stations <- cdf %>%
  filter(station_id %in% c(N %>%
                             pull(station_id))) %>%
  dplyr::select(station_id) %>%
  distinct() %>%
  left_join(md, by = "station_id")

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  geom_sf(data = world,fill = NA,color = "grey") +
  geom_point(data = stations ,
             aes(x = lon, y = lat,
                 color = source), alpha = 0.7,
             shape = 16, size = 1) +
  theme_map() +
  theme(panel.grid.major = element_blank(),
        text = element_text(size = 24)) +
  scale_fill_manual(values = c("white",c("#72a83d"),"darkgreen")) +
  scale_x_continuous(limits = c(-20,60),expand = c(0,0)) +
  scale_y_continuous(limits = c(-25,25),expand = c(0,0))


Climato <- cdf %>%
  filter(station_id %in% c(N %>%
                             pull(station_id))) %>%
  group_by(source,station_id,month) %>%
  summarise(pre.m = mean(pre),
            .groups = "keep") %>%
  group_by(source,station_id) %>%
  summarise(MAP = sum(pre.m),
            N = sum(pre.m < 100),
            .groups = "keep") %>%
  left_join(md, by = "station_id")

ggplot() +
  geom_sf(data = world,fill = NA,color = "grey") +
  geom_point(data = Climato ,
             aes(x = lon, y = lat,
                 color = MAP), alpha = 0.7,
             shape = 16, size = 1) +
  theme_map() +
  theme(panel.grid.major = element_blank(),
        text = element_text(size = 24)) +
  scale_color_gradient2(limits = c(0,2000),
                        midpoint = 1000) +
  scale_x_continuous(limits = c(-20,60),expand = c(0,0)) +
  scale_y_continuous(limits = c(-25,25),expand = c(0,0))

cdf2plot <- cdf %>%
  filter(station_id %in% c(N %>%
                             pull(station_id))) %>%
  group_by(source,station_id,month) %>%
  summarise(pre.m = mean(pre),
            .groups = "keep") %>%
  left_join(md, by = "station_id") %>%
  filter(lat < 0, lat > -5,
         lon > 15, lon < 25) %>%
  group_by(station_id) %>%
  filter(length(unique(month)) == 12)


ggplot(data = cdf2plot) +
  geom_line(aes(x = month,y = pre.m,
                color = as.factor(station_id))) +
  theme_bw() +
  guides(color = "none")

###############################################################
# Trend

cdf <- long_df %>%
  group_by(station_id) %>%
  filter(year >= 1960) %>%
  mutate(N.years = length(unique(year))) %>%
  filter(N.years > 10) %>%
  mutate(time = year + (month -1/2)/12) %>%
  mutate(slope = coef(lm(pre ~ time))[2],
         pval = summary((lm(pre ~ time)))[["coefficients"]][2,4])

hist(cdf$slope)

stations.slopes <- cdf %>%
  left_join(md, by = "station_id") %>%
  dplyr::select(-c(year,month,pre))

stations2plot <- stations.slopes %>%
  dplyr::select(station_id,lon,lat, pval, slope) %>%
  ungroup() %>%
  distinct()

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  geom_sf(data = world,fill = NA,color = "grey") +
  geom_point(data = stations2plot ,
             aes(x = lon, y = lat,
                 color = slope*10), alpha = 0.7,
             shape = 16, size = 1) +
  theme_map() +
  theme(panel.grid.major = element_blank(),
        text = element_text(size = 24)) +
  scale_color_gradient2(limits = c(-1,1)*0.25*10, oob = scales::squish) +
  scale_x_continuous(limits = c(-20,60),expand = c(0,0)) +
  scale_y_continuous(limits = c(-25,25),expand = c(0,0))



ggplot(data )

rm(list = ls())

library(ncdf4)
library(reshape2)
library(lubridate)
library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)

nc <- nc_open("~/Downloads/Congo_Q_2004-2025.nc")

Discharge <- ncvar_get(nc,"discharge")
dim(Discharge)

times <- as.Date(ncvar_get(nc,"time"),
                 origin = "0000-01-01")

nc_close(nc)

df <- melt(Discharge) %>%
  rename(station = Var1,
         time = Var2) %>%
  mutate(time = times[time],
         year = year(time),
         month = month(time)) %>%
  ungroup() %>%
  mutate(value = case_when(value < 0 ~ 0,
                           TRUE ~ value))

stations <- unique(df$station) ; Nstation = 100
Samples <- sample(unique(stations),Nstation)


df_extremes <- df %>%
  group_by(station) %>%
  mutate(
    q95 = quantile(value, 0.95, na.rm = TRUE),
    q05 = quantile(value, 0.05, na.rm = TRUE),
    flood = value > q95,
    drought = value < q05
  ) %>%
  mutate(type = case_when(flood ~ "flood",
                          drought ~"drought",
                          TRUE ~ "else"))

df.month <- df %>%
  group_by(station,year,month) %>%
  summarise(value = mean(value,na.rm = TRUE),
            .groups = "keep")


ggplot(data = df.month %>%
         filter(station %in% Samples)) +
  geom_line(aes(x = year + (month-1/2)/12, y = value,
                group = as.factor(station))) +
  theme_bw() +
  guides(color = "none")


df_extremes_months <- df.month %>%
  group_by(station,month) %>%
  mutate(
    q50 = quantile(value, 0.5, na.rm = TRUE),
    qm = mean(value,na.rm = TRUE)) %>%
  group_by(station,month) %>%
  mutate(
    q95 = quantile(value, 0.975, na.rm = TRUE),
    q05 = quantile(value, 0.025, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(flood = value > q95,
         drought = value < q05,
         type = case_when(flood ~ "flood",
                          drought ~"drought",
                          TRUE ~ "else"))

df_extremes_months_sum <- df_extremes_months %>%
  group_by(year,month) %>%
  summarise(drought = 100*sum(drought)/n(),
            flood = 100*sum(flood)/n(),
            .groups = "keep") %>%
  pivot_longer(cols = c(drought,flood),
               values_to = "frac",
               names_to = "type")

ggplot(data = df_extremes_months_sum) +
  geom_line(aes(x = year + (month - 1/2)/12,
                y = frac, color = type)) +
  theme_bw()

df_extremes_months_sum_large <- df_extremes_months %>%
  ungroup() %>%
  filter(qm > 5000) %>%
  group_by(year,month) %>%
  summarise(drought = 100*sum(drought)/n(),
            flood = 100*sum(flood)/n(),
            .groups = "keep") %>%
  pivot_longer(cols = c(drought,flood),
               values_to = "frac",
               names_to = "type")


AAA <- df_extremes_months %>%
  filter(year == 2024) %>%
  group_by(station) %>%
  summarise(reduc = mean((value - qm)/qm,na.rm = TRUE),
            .groups = "keep")


quantile(AAA$reduc,c(0.025,0.5,0.975),na.rm = TRUE)

ggplot(data = df_extremes_months_sum_large) +
  geom_line(aes(x = year + (month - 1/2)/12,
                y = frac, color = type)) +
  theme_bw()

large.stations <- df_extremes_months %>%
  ungroup() %>%
  filter(qm > 5000)  %>%
  pull(station) %>% unique()


ggplot(data = df_extremes_months %>%
         filter(station %in% large.stations)) +
  geom_line(aes(x = year + (month-1/2)/12, y = value,
                color = type,
                group = as.factor(station))) +
  theme_bw() +
  scale_color_manual(values = c("red","black","blue")) +
  guides(color = "none")



basins <- read_sf("~/Downloads/Shapefile/Centroids.shp")
Kin <- basins %>%
  mutate(dist = sqrt( (Ycen + 4.32)**2 + (Xcen - 15.31)**2)) %>%
  arrange(dist) %>%
  slice_head(n = 5)


ggplot(data = df_extremes_months %>%
         filter(station %in% Kin$Mini)) +
  geom_line(aes(x = year + (month-1/2)/12, y = value,
                group = as.factor(station),
                color = type)) +
  theme_bw() +
  facet_wrap(~as.factor(station), scales = "free") +
  scale_color_manual(values = c("red","black","blue")) +
  guides(color = "none")

basins_ext <- basins %>%
  rename(station = Mini) %>%
  left_join(df_extremes_months %>%
              filter(year %in% 2023:2024) %>%
              group_by(station) %>%
              summarise(type = case_when(any(drought) & any(flood) ~ "both",
                                         any(drought) ~ "drought",
                                         any(flood) ~ "flood",
                                         TRUE ~ "none")),
            by = "station")


basin <- read_sf("./data/shapefiles/CongoBasin.shp")

ggplot(data = basins_ext) +
  geom_sf(data = basin,fill = NA, color = "black") +
  geom_point(aes(x = Xcen, y = Ycen,
                 color = type),
             size = 0.2) +
  theme_bw()


basins_ext_drought <-
  df_extremes_months %>%
  filter(year %in% 2023:2024,
         drought | flood) %>%
  mutate(time = as.Date(paste0(year,"/",month,"/01"))) %>%
  left_join(basins %>%
      rename(station = Mini),
      by = "station") %>%
  mutate(q.norm = 100*value/qm)

ggplot(data = basins_ext_drought %>%
         filter(year == 2023)) +
  geom_sf(data = basin,fill = NA, color = "black") +
  geom_point(aes(x = Xcen, y = Ycen,
                 color = q.norm),
             size = 0.2) +
  scale_color_gradient2(midpoint = 100, limits = c(0,200),oob = scales::squish,
                        low = "darkred",high = "darkblue") +
  facet_wrap(year ~ month, ncol = 4) +
  theme_bw()

ggplot(data = basins_ext_drought %>%
         filter(year == 2024)) +
  geom_sf(data = basin,fill = NA, color = "black") +
  geom_point(aes(x = Xcen, y = Ycen,
                 color = q.norm),
             size = 0.2) +
  scale_color_gradient2(midpoint = 100, limits = c(0,200),oob = scales::squish,
                        low = "darkred",high = "darkblue") +
  facet_wrap(year ~ month, ncol = 4) +
  theme_bw()

basins_ext_drought %>%
  filter(q.norm > 100)


hist(basins$Area__km2_)

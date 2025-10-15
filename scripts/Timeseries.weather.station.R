rm(list = ls())

library(dplyr)
library(ggplot2)
library(tidyr)

raw.metadata <- readRDS("./outputs/all.metadata.RDS")
all.metadata <- raw.metadata %>%
  rename(station_id = unified_id) %>%
  dplyr::select(station_id,lon,lat,elevation)

hist(all.metadata %>%
       pull(elevation))

all.metadata %>%
  summarise(frac = sum(is.na(elevation))/n())


world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  guides(fill = "none") +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  geom_point(data = all.metadata %>%
               ungroup() %>%
               dplyr::select(station_id,lon,lat) %>%
               distinct(),
             aes(x = lon, y = lat), size = 0.6, alpha = 0.5) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  theme_bw()

ggplot() +
  stat_bin2d(data = all.metadata %>%
               ungroup() %>%
               dplyr::select(station_id,lon,lat),
             alpha = 0.8,
             aes(lon, lat),
             binwidth = c(2, 2)) +
  labs(fill = "Obs count", x = "Longitude", y = "Latitude") +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  scale_fill_gradient2(mid = "grey") +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  theme_bw()

######################################################################

raw.data <- readRDS("./outputs/all.data.RDS")
all.data <- raw.data %>%
  rename(station_id = unified_id) %>%
  dplyr::select(-ends_with("_source")) %>%
  filter(year <= 2024)


all.data.long <- all.data %>%
  pivot_longer(cols = any_of(c("pre","tas","tasmin","tasmax","dewpoint")),
               names_to = "variable",
               values_to = "value") %>%
  na.omit()




ggplot(data = all.data.long) +
  geom_density(aes(x = value),
               alpha = 0.5, color = NA, fill = "grey") +
  facet_wrap(~ variable, scales = "free") +
  theme_bw()

ggplot(data = all.data.long %>%
         filter(variable == "pre",
                value > 0)) +
  geom_density(aes(x = value),
               alpha = 0.5, color = NA, fill = "grey") +
  scale_x_continuous(limits = c(0,200)) +
  theme_bw() +
  guides(fill = "none")


all.data.long %>%
  filter(variable == "pre",
         value != 0) %>%
  summarise(min = min(value),
            m = mean(value),
            max = max(value),
            qhigh = quantile(value, 0.95),
            med = median(value))

BIN <- 2

all.data.long.wbin <- all.data.long %>%
  left_join(all.metadata,
            by = c("station_id")) %>%
  mutate(lon_bin = round(lon/BIN)*BIN,
         lat_bin = round(lat/BIN)*BIN)

all.data.binned <- all.data.long.wbin  %>%
  group_by(variable,lon_bin,lat_bin,year,month) %>%
  summarise(value.m = mean(value,na.rm = TRUE),
            .groups = "keep")

all.data.long.coord.sum <-
  all.data.binned %>%
  group_by(variable,lon_bin,lat_bin) %>%
  summarise(Nobs = length(value.m),
            .groups = "keep")

ggplot() +
  labs(fill = "Obs count", x = "", y = "") +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  geom_raster(data = all.data.long.coord.sum,
              alpha = 0.8,
              aes(lon_bin, lat_bin,
                  fill = Nobs)) +
  scale_fill_gradient2(mid = "grey") +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  theme_bw() +
  facet_wrap(~ variable)

all.data.sum <-
  all.data.long %>%
  group_by(variable,year) %>%
  summarise(Nstation = length(unique(station_id)),
            Nobs = length((value)),
            .groups = "keep")

ggplot(data = all.data.sum) +
  geom_line(aes(x = year,
                y = Nstation)) +
  theme_bw() +
  facet_wrap(~variable)

ggplot(data = all.data.sum) +
  geom_line(aes(x = year,
                y = Nobs)) +
  theme_bw() +
  facet_wrap(~variable)


all.data.sum <- all.data.binned %>%
  group_by(variable,lon_bin,lat_bin,month) %>%
  summarise(value.m = mean(value.m,na.rm = TRUE),
            .groups = "keep") %>%
  group_by(variable,lon_bin,lat_bin) %>%
  summarise(value = case_when(variable[1] == "pre" ~ sum(value.m),
                              TRUE ~ mean(value.m)),
            .groups = "keep") %>%
  pivot_wider(names_from = variable,
              values_from = value)

hist(all.data.sum$pre)

ggplot() +
  geom_raster(data = all.data.sum %>%
                dplyr::select(lon_bin,lat_bin,pre) %>%
                na.omit(),
              alpha = 0.8,
              aes(lon_bin, lat_bin,
                  fill = pre)) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  scale_fill_gradient2(mid = "grey",midpoint = 1500,
                       limits = c(0,3000),oob = scales::squish) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  theme_bw()

ggplot() +
  geom_raster(data = all.data.sum %>%
                dplyr::select(lon_bin,lat_bin,tas) %>%
                na.omit() %>%
                filter(tas > 0),
              alpha = 0.8,
              aes(lon_bin, lat_bin,
                  fill = tas)) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  scale_fill_gradient2(mid = "grey",midpoint = 23) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  theme_bw()




ggplot() +
  geom_raster(data = all.data.sum %>%
                dplyr::select(lon_bin,lat_bin,tasmin) %>%
                na.omit,
              alpha = 0.8,
              aes(lon_bin, lat_bin,
                  fill = tasmin)) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  scale_fill_gradient2(mid = "grey",midpoint = 17) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  theme_bw()


ggplot() +
  geom_raster(data = all.data.sum %>%
                dplyr::select(lon_bin,lat_bin,tasmax) %>%
                na.omit,
              alpha = 0.8,
              aes(lon_bin, lat_bin,
                  fill = tasmax)) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  scale_fill_gradient2(mid = "grey",midpoint = 30) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  theme_bw()

ggplot() +
  geom_raster(data = all.data.sum %>%
                dplyr::select(lon_bin,lat_bin,dewpoint) %>%
                na.omit,
              alpha = 0.8,
              aes(lon_bin, lat_bin,
                  fill = dewpoint)) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  scale_fill_gradient2(mid = "grey",limits = c(0,25),
                       midpoint = 12.5,
                       oob = scales::squish) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  theme_bw()

MAP <- all.data.binned %>%
  group_by(variable,year,lon_bin,lat_bin,month) %>%
  summarise(value.m = mean(value.m,na.rm = TRUE),
            .groups = "keep") %>%
  group_by(variable,year,lon_bin,lat_bin) %>%
  mutate(N = n()) %>%
  filter(N==12) %>%
  group_by(variable,year,lon_bin,lat_bin) %>%
  summarise(value = case_when(variable[1] == "pre" ~ sum(value.m),
                              TRUE ~ mean(value.m)),
            .groups = "keep")

MAP.lm <- MAP %>%
  filter(year >= 1960) %>%
  group_by(variable,lon_bin,lat_bin) %>%
  mutate(N = n()) %>%
  filter(N > 10) %>%
  summarise(slope = coef(lm(value ~ year))[2],
            p.val = summary(lm(value ~ year))[["coefficients"]][2,4],
            .groups = "keep")

ggplot(data = MAP %>%
         filter(lon_bin == -8,
                lat_bin == 4,
                variable == "pre"),
       aes(x = year, y = value)) +
  geom_point() +
  stat_smooth(method = "lm") +
  theme_bw()

df2plot2test <- all.data.long.wbin %>%
  filter(lon_bin == 0, lat_bin == 4) %>%
  filter(variable == "pre") %>%
  group_by(station_id,year) %>%
  mutate(N = n()) %>%
  ungroup() %>%
  filter(N == 12) %>%
  group_by(station_id,year) %>%
  summarise(MAP = sum(value),
            .groups = "keep")

ggplot(data = df2plot2test %>%
         filter(year >= 1960),
       aes(x = year, y = MAP, color = interaction(station_id))) +
  geom_point() +
  stat_smooth(method = "lm") +
  stat_smooth(color = "black",method = "lm") +
  theme_bw() +
  guides(color = "none", fill = "none")


MAP.lm %>%
  filter(variable == "pre") %>%
  filter(p.val < 0.05,
         slope < 0)

ggplot() +

  geom_raster(data = MAP.lm %>%
                filter(variable == "pre"),
              alpha = 0.8,
              aes(lon_bin, lat_bin,
                  fill = slope*10)) +
  geom_point(data = MAP.lm %>%
               filter(variable == "pre",
                      p.val <= 0.05),
             aes(lon_bin, lat_bin)) +
  scale_fill_gradient2(mid = "grey",
                       limits = c(-1,1)*100,
                       oob = scales::squish) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  theme_bw() +
  facet_wrap(~ variable)

slopes.MAP <-
  as.vector(MAP.lm %>%
              filter(variable == "pre") %>%
              pull(slope)*10)
hist(slopes.MAP)

t.test(slopes.MAP, mu = 0)


ggplot() +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  geom_raster(data = MAP.lm %>%
                filter(variable == "tas"),
              alpha = 0.8,
              aes(lon_bin, lat_bin,
                  fill = slope)) +
  geom_point(data = MAP.lm %>%
               filter(variable == "tas",
                      p.val <= 0.05),
             aes(lon_bin, lat_bin)) +
  scale_fill_gradient2(mid = "grey",limits = c(-1,1)/20,
                       oob = scales::squish) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  theme_bw() +
  facet_wrap(~ variable)



slopes.MAT <-
  as.vector(MAP.lm %>%
              filter(variable == "tas") %>%
              pull(slope)*10)

hist(slopes.MAT)

t.test(slopes.MAT, mu = 0)


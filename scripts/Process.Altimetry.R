rm(list = ls())

library(Drying.CB)
library(dplyr)
library(ggplot2)
library(lubridate)
library(sf)

Dir <- "/Users/felicien/Documents/projects/Drying.CB/data/Altimetry/Alti_River_Water_Level_opperational_20_09_2025/"
files <- list.files(Dir,full.names = TRUE, pattern = "*.txt")

#cfile <- files[1]

all.md <- all.data <-
  data.frame()
for (cfile in files){

  print(basename(cfile))

  L <- read_hydroweb(cfile)
  csplit <- strsplit(basename(cfile),"\\_")[[1]]

  all.md <- bind_rows(all.md,
                      data.frame(file = basename(cfile),
                                 basin = csplit[4],
                                 ID = L$metadata$ID,
                                 lat = L$metadata$REFERENCE.LATITUDE,
                                 lon = L$metadata$REFERENCE.LONGITUDE))
  all.data <- bind_rows(all.data,
                        L$data %>%
                          mutate(ID = L$metadata$ID))

}

kin <- all.md %>%
  filter(basin == "CONGO") %>%
  mutate(dist = sqrt( (lat + 4.32)**2 + (lon - 15.31)**2)) %>%
  arrange(dist) %>%
  slice_head(n = 5)

all.data.month <- all.data %>%
  mutate(year = year(datetime),
         month = month(datetime)) %>%
  group_by(ID,year,month) %>%
  summarise(orthometric_height_m = mean(orthometric_height_m),
            .groups = "keep")

ggplot(data = all.data %>%
         filter(ID %in% kin$ID) %>%
         ungroup()) +
  geom_line(aes(x = as.Date(date), y = orthometric_height_m,
                color = as.factor(ID))) +
  theme_bw()

df2plot <- all.data.month %>%
  filter(ID == 104620) %>%
  ungroup()

df2plot.m <- df2plot %>%
  group_by(month) %>%
  summarise(m = mean(orthometric_height_m),
            S = sd(orthometric_height_m),
            .groups = "keep")

ggplot() +
  geom_ribbon(data = df2plot.m,
            aes(x = month,ymin = m - S, ymax = m + S,
                y = m),color = NA, alpha = 0.4,
            fill = "darkgrey") +
  geom_line(data = df2plot.m,
            aes(x = month,
                y = m),
            color = "darkgrey") +
  geom_line(data = df2plot %>%
              filter(year == 2024),
            aes(x = month,
                y = orthometric_height_m),
            color = "red") +
  labs(x = "", y = "Water level") +
  scale_x_continuous(breaks = 1:12) +
  theme_bw() +
  theme(text = element_text(size = 22))


# Example of one curve

cdf.norm <- all.data.month %>%
  group_by(ID) %>%
  mutate(length(unique(year)) > 5) %>%
  ungroup() %>%
  filter(ID == sample(ID,1)) %>%
  group_by(month) %>%
  mutate(M = mean(orthometric_height_m,na.rm = TRUE),
         SD = sd(orthometric_height_m,na.rm = TRUE)) %>%
  mutate(z = (orthometric_height_m - M)/SD,
         time = year + (month - 1/2)/12)


ggplot(data = cdf.norm) +
  geom_rect(
    data = data.frame(x = 2024, xend = 2025, y = -Inf, yend = Inf),
    aes(xmin = x, xmax = xend, ymin = y, ymax = yend),
    fill = "grey", alpha = 0.4) +
  geom_line(aes(x = time, y = z)) +
  geom_hline(yintercept = 0,linetype = 2) +
  theme_bw()





all.data.anomaly <- all.data.month %>%
  group_by(ID,month) %>%
  mutate(Mean = mean(orthometric_height_m,na.rm = TRUE),
         Q95 = quantile(orthometric_height_m,0.95,na.rm = TRUE),
         Q05 = quantile(orthometric_height_m,0.05,na.rm = TRUE),
         Sd = sd(orthometric_height_m,na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(anomaly = orthometric_height_m - Mean,
         anomaly.norm = (orthometric_height_m - Mean)/Sd,
         type = case_when(orthometric_height_m > Q95 ~ "flood",
                          orthometric_height_m < Q05 ~ "drought",
                          TRUE ~ "else"))

data.years <- all.data.anomaly %>%
  group_by(ID) %>%
  mutate(N = length(which(!is.na(orthometric_height_m)))) %>%
  filter(N >= (10*12)) %>%
  group_by(ID,year) %>%
  summarise(min.anomaly.norm = max(anomaly.norm) - min(anomaly.norm),
            .groups = "keep")

ggplot() +
  geom_density(data = data.years,
               aes(group = as.factor(year), x = min.anomaly.norm), alpha = 0.5,
               size = 0.25) +
  geom_density(data = data.years %>%
                 filter(year %in% c(2015,2024:2025)),
               aes(x = min.anomaly.norm, color = as.factor(year)), size = 1) +
  geom_vline(xintercept = 0,linetype = 2) +
  theme_bw()

all.data.2024 <- all.data.anomaly %>%
  filter(year == 2024)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

df2plot <- all.data.2024 %>%
  left_join(all.md,
            by = "ID")

ggplot(data = df2plot) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  geom_point(aes(x = lon, y = lat,
                 color = anomaly),
             size = 1) +
  coord_sf(xlim = c(5, 35), ylim = c(-15, 10), expand = FALSE) +
  theme_bw() +
  facet_wrap(~ (month)) +
  scale_color_gradient2(limits = c(-1,1)*2,
                       oob = scales::squish) +
  labs(x = "", y = "", color = "") +
  theme(legend.position = "right")


ggplot(data = df2plot) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  geom_point(aes(x = lon, y = lat,
                 color = type)) +
  coord_sf(xlim = c(5, 35), ylim = c(-15, 10), expand = FALSE) +
  theme_bw() +
  facet_wrap(~ (month)) +
  labs(x = "", y = "", color = "") +
  theme(legend.position = "right")


df.ts.type <- all.data.anomaly %>%
  group_by(year,month) %>%
  summarise(drought = 100*sum(type == "drought")/n(),
            flood = 100*sum(type == "flood")/n(),
            .groups = "keep") %>%
  pivot_longer(cols = c(drought,flood),
               values_to = "frac",
               names_to = "type")

ggplot(data = df.ts.type) +
  geom_line(aes(x = year + (month - 1/2)/12,
                y = frac, color = type)) +
  theme_bw()


all.data.2024.ext <- all.data.2024 %>%
  group_by(ID) %>%
  filter(anomaly %in% c(min(anomaly),max(anomaly))) %>%
  left_join(all.md,
            by = "ID") %>%
  mutate(trim = case_when(month %in% c(12,1:2) ~ "DFJ",
                          month %in% c(3:5) ~ "MAM",
                          month %in% c(6:8) ~ "JJA",
                          month %in% c(9:11) ~ "SON",
                          TRUE ~ NA_character_)) %>%
  mutate(trim2 = factor(case_when(month %in% c(1:3) ~ "JFM",
                          month %in% c(4:6) ~ "AMJ",
                          month %in% c(7:9) ~ "JAS",
                          month %in% c(10:12) ~ "OND",
                          TRUE ~ NA_character_),
                        levels = c("JFM","AMJ","JAS","OND")))

ggplot(data = all.data.2024.ext %>%
         group_by(ID) %>%
         filter(anomaly == min(anomaly))) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  geom_point(aes(x = lon, y = lat,
                 color = anomaly)) +
  coord_sf(xlim = c(5, 35), ylim = c(-15, 10), expand = FALSE) +
  theme_bw() +
  facet_wrap(~(trim2)) +
  scale_color_gradient2(limits = c(-1,1)*2,
                        oob = scales::squish) +
  labs(x = "", y = "", color = "") +
  theme(legend.position = "right")


ggplot() +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  geom_point(data = all.data.2024 %>%
               group_by(ID) %>%
               filter(anomaly == min(anomaly)) %>%
               left_join(all.md,
                         by = "ID"),
             aes(x = lon, y = lat,
                 color = anomaly)) +
  coord_sf(xlim = c(5, 35), ylim = c(-15, 10), expand = FALSE) +
  theme_bw() +
  scale_color_gradient2(limits = c(-1,1)*2,
                        oob = scales::squish) +
  labs(x = "", y = "", color = "") +
  theme(legend.position = "right")

station_season_2024 <- all.data.anomaly %>%
  filter(year == 2024) %>%
  mutate(season = factor(case_when(
    month %in% 1:3 ~ "JFM",
    month %in% 4:6 ~ "AMJ",
    month %in% 7:9 ~ "JAS",
    month %in% 10:12 ~ "OND"
  ), levels = c("JFM", "AMJ", "JAS", "OND"))) %>%
  group_by(ID,season) %>%
  summarise(anomaly_season = mean(anomaly, na.rm = TRUE),
            anomaly.norm_season = mean(anomaly.norm, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(all.md,
            by = "ID")

ggplot(data = station_season_2024) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  geom_point(aes(x = lon, y = lat,
                 color = anomaly.norm_season)) +
  coord_sf(xlim = c(5, 35), ylim = c(-15, 10), expand = FALSE) +
  theme_bw() +
  facet_wrap(~(season)) +
  scale_color_gradient2(limits = c(-1,1)*2,
                        oob = scales::squish) +
  labs(x = "", y = "", color = "") +
  theme(legend.position = "right")


################################################################################
sf_use_s2(FALSE)

basins_raw <- read_sf("~/Downloads/hybas_af_lev01-12_v1c/hybas_af_lev05_v1c.shp")
e <- st_bbox(c(xmin = 5, xmax = 35, ymin = -15, ymax = 10), crs = st_crs(basins_raw))
basins_crop <- st_crop(basins_raw, e)

basins_valid <- st_make_valid(basins_crop)
basins <- st_simplify(basins_valid, dTolerance = 0.05)
basins$HYBAS_ID <- as.factor(basins$HYBAS_ID)
plot((basins["HYBAS_ID"]))

stations_sf <- all.md %>%
  dplyr::select(ID,lon,lat) %>%
  distinct() %>%
  st_as_sf(coords = c("lon", "lat"),
           crs = 4326, remove = FALSE)

# Make sure CRS matches
basins <- st_transform(basins, st_crs(stations_sf))

stations_basin <- st_join(stations_sf, basins, left = TRUE)

# Back to plain data.frame for joins with time series
all.md.basin <- stations_basin %>%
  st_drop_geometry() %>%
  distinct(ID, .keep_all = TRUE)

idx <- st_nearest_feature(stations_sf, basins)
all.md.basin$HYBAS_ID <- basins$HYBAS_ID[idx]

group_var <- "HYBAS_ID"

# Keep only stations that were matched to a basin/subbasin
selected <- all.md.basin %>%
  filter(!is.na(.data[[group_var]])) %>%
  pull(ID)

Sum <- all.data.anomaly %>%
  left_join(all.md.basin, by = "ID") %>%
  filter(ID %in% selected) %>%
  filter(!is.na(.data[[group_var]])) %>%
  group_by(.data[[group_var]], year, month) %>%
  summarise(
    anomaly.m   = mean(anomaly, na.rm = TRUE),
    anomaly.med = median(anomaly, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(group = all_of(group_var))

ggplot(data = all.md.basin %>%
         filter(HYBAS_ID %in% c("1051163840",
                                ""))) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  geom_sf(data = basins %>%
            filter(HYBAS_ID %in% unique(all.md.basin$HYBAS_ID)), fill = NA, color = "black", linewidth = 1) +
  geom_point(aes(x = lon,y = lat,
                 color = as.factor(HYBAS_ID))) +
  coord_sf(xlim = c(5, 35), ylim = c(-15, 10), expand = FALSE) +
  theme_bw() +
  guides(color = "none")

n_stations <- all.md.basin %>%
  filter(!is.na(.data[[group_var]])) %>%
  distinct(ID, .data[[group_var]]) %>%
  count(.data[[group_var]], name = "n") %>%
  rename(group = all_of(group_var))

Sum <- all.data.anomaly %>%
  left_join(all.md.basin, by = "ID") %>%
  filter(ID %in% selected) %>%
  filter(!is.na(.data[[group_var]])) %>%
  group_by(.data[[group_var]], year, month) %>%
  summarise(
    anomaly.m   = mean(anomaly, na.rm = TRUE),
    anomaly.med = median(anomaly, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(group = all_of(group_var)) %>%
  left_join(n_stations, by = "group") %>%
  mutate(group_lab = paste0(group, " (n=", n, ")"))

Lines <- all.data.anomaly %>%
  left_join(all.md.basin, by = "ID") %>%
  filter(ID %in% selected) %>%
  filter(!is.na(.data[[group_var]])) %>%
  mutate(group = .data[[group_var]]) %>%
  left_join(n_stations, by = "group") %>%
  mutate(group_lab = paste0(group, " (n=", n, ")"))


# ------------------------------------------------------------------------------
# 5. Plot: one panel per basin/subbasin
# ------------------------------------------------------------------------------

ggplot() +
  geom_rect(
    data = data.frame(x = 2024, xend = 2025, y = -Inf, yend = Inf),
    aes(xmin = x, xmax = xend, ymin = y, ymax = yend),
    fill = "grey", alpha = 0.4) +
  geom_line(
    data = Lines,
    aes(x = year + (month - 0.5) / 12,
        y = anomaly,
        group = ID,
        color = as.factor(ID)),
    alpha = 0.35) +
  geom_line(
    data = Sum,
    aes(x = year + (month - 0.5) / 12,
        y = anomaly.m),
    color = "black",
    linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~group_lab, scales = "free_y") +
  scale_y_continuous(limits = c(-1,1)*3) +
  theme_bw() +
  guides(color = "none") +
  labs(x = "", y = "Water level anomaly (m)")


select.station <- all.md %>%
  mutate(dist = (lat - (-4.322447)**2+
                   (lon - 15.307045)**2)) %>%
  arrange(dist) %>%
  slice_head(n = 1)

all.data %>%
  filter(ID == select.station[["ID"]])

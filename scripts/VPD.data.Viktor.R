rm(list = ls())

library(dplyr)
library(ggplot2)
library(Drying.CB)
library(scales)

raw.data <- readRDS("./outputs/all.data.RDS") %>%
  filter(!is.na(dewpoint) & !is.na(tas)) %>%
  filter(year >= 2006 & year <= 2025)

unique(raw.data$dewpoint_source)             # Sources
length(unique(raw.data$unified_id))          # Stations
nrow(raw.data)                               # Monthly obs.
print(paste0(min(raw.data$year),"-",max(raw.data$year)))

data.years <- raw.data %>%
  group_by(year) %>%
  summarise(N = n(),.groups = "keep")

ggplot(data = data.years) +
  geom_line(aes(x = year, y = N)) +
  theme_bw()

raw.data.vpd <- raw.data %>%
  mutate(vpd = vpd_kPa(tas,dewpoint))

data.vpd.sum <- raw.data.vpd %>%
  group_by(unified_id) %>%
  summarise(vpd.m = mean(vpd,na.rm = TRUE),
            .groups = "keep")

# group_by(station)
hist(data.vpd.sum$vpd.m)


raw.metadata <- readRDS("./outputs/all.metadata.RDS")
all.metadata <- raw.metadata %>%
  filter(unified_id %in% unique(raw.data$unified_id)) %>%
  rename(station_id = unified_id) %>%
  dplyr::select(station_id,lon,lat,elevation)

data.vpd.sum.pos <- data.vpd.sum %>%
  rename(station_id = unified_id) %>%
  left_join(all.metadata,
            by = "station_id")

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  guides(fill = "none") +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  geom_point(data = data.vpd.sum.pos %>%
               ungroup() %>%
               dplyr::select(station_id,lon,lat,vpd.m) %>%
               distinct(),
             aes(x = lon, y = lat,
                 color = vpd.m), size = 0.6, alpha = 0.5) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_color_gradient2(midpoint = 2,
                        low = muted("blue"),
                        high = muted('red')) +
  theme_bw()





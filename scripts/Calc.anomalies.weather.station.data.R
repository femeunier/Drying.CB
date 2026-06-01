rm(list = ls())

library(dplyr)
library(ggplot2)
library(tidyr)
library(sf)
library(slider)

raw.metadata <- readRDS("./outputs/all.metadata.RDS")
all.metadata <- raw.metadata %>%
  rename(station_id = unified_id) %>%
  dplyr::select(station_id,lon,lat,elevation)

year.min <- 2000 ; year.max <- 2024

######################################################################

raw.data <- readRDS("./outputs/all.data.RDS")
all.data <- raw.data %>%
  rename(station_id = unified_id) %>%
  dplyr::select(-ends_with("_source")) %>%
  filter(year>= year.min & year <= year.max)


Mask <- read_sf("./data/shapefiles/Rainforests.shp")

pts <- st_as_sf(all.metadata,
                coords = c("lon", "lat"), crs = 4326, remove = FALSE)

mask <- st_make_valid(Mask)
pts  <- st_transform(pts, st_crs(mask))
pts_in <- st_filter(pts, mask, .predicate = st_within)
stations_in_mask <- st_drop_geometry(pts_in)

all.data.long <- all.data %>%
  filter(station_id %in% c(stations_in_mask$station_id)) %>%
  pivot_longer(cols = any_of(c("pre","tas","tasmin","tasmax","dewpoint")),
               names_to = "variable",
               values_to = "value") %>%
  na.omit()


world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  guides(fill = "none") +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  geom_sf(data = Mask, fill = NA, color = "black", linewidth = 1) +
  geom_point(data = stations_in_mask %>%
               ungroup() %>%
               dplyr::select(station_id,lon,lat) %>%
               distinct(),
             aes(x = lon, y = lat), size = 0.6, alpha = 0.5) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  theme_bw()

all.data.long.compl <- all.data.long %>%
  complete(station_id = unique(all.data.long$station_id),
           year = year.min:year.max,
           month = 1:12,
           variable = unique(all.data.long$variable),
           fill = list(value = NA)) %>%
  arrange(station_id,variable,year,month) %>%
  group_by(station_id) %>%
  mutate(Ntot = length(which(!is.na(value)))) %>%
  group_by(station_id,variable) %>%
  mutate(roll12_mean =
           slide_dbl(value, mean, .before = 11, .complete = TRUE, na.rm = TRUE),
         N = length(which(!is.na(value))),
         .groups = "keep") %>%
  ungroup()


ggplot(data = all.data.long.compl %>%
         filter(Ntot == max(Ntot)),
       aes(x = year + (month - 1/2)/12,
           y = roll12_mean)) +
  geom_line() +
  facet_wrap(~ variable, scales = "free") +
  stat_smooth(method = "lm",
              se = FALSE) +
  theme_bw()

slopes <- all.data.long.compl %>%
  filter(N > 120) %>%
  mutate(time = year + (month - 1/2)/12) %>%
  group_by(station_id,variable) %>%
  summarise(slope = coef(lm(roll12_mean ~ time))[2],
            p.val = summary(lm(roll12_mean ~ time))[["coeficients"]][2,4],
            .groups = "keep") %>%
  left_join(all.metadata %>%
              dplyr::select(lon,lat,station_id),
            by = c("station_id"))

ggplot(data = slopes) +
  geom_density(aes(x = slope*10)) +
  geom_vline(xintercept = 0, linetype = 2) +
  facet_wrap(~variable, scales = "free") +
  theme_bw()


ggplot() +
  guides(fill = "none") +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  # geom_sf(data = Mask, fill = NA, color = "black", linewidth = 1) +
  geom_point(data = slopes %>%
               filter(variable == "pre"),
             aes(x = lon, y = lat,
                 fill = slope),
             color = "black", shape = 21,
             size = 1, alpha = 1) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2() +
  facet_wrap(~ variable) +
  theme_bw()


ggplot() +
  guides(fill = "none") +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  # geom_sf(data = Mask, fill = NA, color = "black", linewidth = 1) +
  geom_point(data = slopes %>%
               filter(variable == "tas"),
             aes(x = lon, y = lat,
                 fill = slope),
             color = "black", shape = 21,
             size = 1, alpha = 1) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2() +
  facet_wrap(~ variable) +
  theme_bw()

data2save <- all.data.long.compl %>%
  filter(!is.na(value)) %>%
  filter(variable %in% c("tas","tasmin","tasmax","pre")) %>%
  group_by(year,month,variable) %>%
  summarise(N = length(value),
            .groups = "keep")

plot(data2save$N)

data2save %>%
  filter(year >= 2020)

saveRDS(data2save,
        "./outputs/Weather.overview.RDS")


# # Anomalies
#
# selected.station <- all.data.long %>%
#   group_by(station_id,variable) %>%
#   summarise(Nhistorical = length(which(year %in% 1961:1990)),
#             Nrecent = length(which(year %in% 2024)),
#             .groups = "keep") %>%
#   filter(Nrecent > 0,
#          Nhistorical >= 180) %>%
#   mutate(station_var = paste0(station_id,"_",variable))
#
# all.data.long.selected <- all.data.long %>%
#   mutate(station_var = paste0(station_id,"_",variable)) %>%
#   filter(station_var %in% c(selected.station[["station_var"]]))
#
#
# all.data.long.selected.anomalies <-
#   all.data.long.selected %>%
#   filter(year %in% c(1961:1990,2024)) %>%
#   mutate(period = case_when(year == 2024 ~ "recent",
#                             TRUE ~ "historical")) %>%
#   group_by(station_id,period,variable,month) %>%
#   summarise(value.m = mean(value,na.rm = TRUE),
#             .groups = "keep") %>%
#   pivot_wider(names_from = "period",
#               values_from = "value.m")

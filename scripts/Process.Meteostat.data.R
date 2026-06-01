rm(list = ls())

library(jsonlite)
library(dplyr)
library(terra)
library(ggplot2)
library(tidyr)
library(raster)

# Read JSON file
stations <- fromJSON("./data/station_data/Meteostat/lite.json", flatten = TRUE)

# Extract lat/lon
coords <- stations %>%
  transmute(id,
            latitude = location.latitude,
            longitude = location.longitude,
            elevation = location.elevation) %>%
  mutate(elevation = case_when(elevation == -9999 ~ NA,
                               TRUE ~ elevation))

Central.African.stations <- coords %>%
  filter(latitude >= -25, latitude <= 25,
         longitude >= -30, longitude <= 65)
r <- rast("/home/femeunier/Downloads/treecover2000_005.tif")  # % tree cover (0-100)
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
r_plot <- r %>% terra::aggregate(4, fun = "mean")
df_r <- as.data.frame(r_plot, xy = TRUE, na.rm = TRUE) %>%
  mutate(mask = case_when(treecover2000_vrt > 50 ~ 1,
                          TRUE ~ 0))
names(df_r)[3] <- "tree"


all.df <- data.frame()
for (irow in seq(1,nrow(Central.African.stations))){
  cid <- Central.African.stations[["id"]][irow]

  op.file <- paste0("./data/station_data/Meteostat/",cid,".csv.gz")
  if (!file.exists(op.file)){
    system2("curl",
            c(paste0("https://data.meteostat.net/monthly/",
                     cid,".csv.gz"),
              "--output",
              op.file))
  }


  cdf <- tryCatch(read.csv(gzfile(paste0("./data/station_data/Meteostat/",cid,".csv.gz"))),
                  error = function(e) NULL)

  if (all(c("year","prcp") %in% colnames(cdf))){
    all.df <- bind_rows(all.df,
                        cdf %>%
                          mutate(id = cid))
  }

}

all.df.selected <- all.df %>%
  dplyr::select(id,year,month,temp,tmin,tmax,prcp)

all.df.selected.long <- all.df.selected %>%
  pivot_longer(cols = -c(id,year,month),
               names_to = "variable",
               values_to = "value") %>%
  na.omit()

stations <- all.df.selected.long %>%
  group_by(id,variable) %>%
  summarise(N = n(),
            Nmonth = length(unique(month)),
            .groups = "keep")

ggplot(data = stations,
       aes(x = N, fill = variable)) +
  geom_density(alpha = 0.5, color = NA) +
  theme_bw()

stations.months <- all.df.selected.long %>%
  group_by(id,variable,month) %>%
  summarise(N = n(),
            .groups = "keep") %>%
  ungroup() %>%
  complete(id = unique(id),
           variable = unique(variable),
           month = 1:12,
           fill = list(N = 0)) %>%
  group_by(id,variable) %>%
  summarise(Nmin = min(N),
            .groups = "keep")


selected.ids <- sort(unique(stations.months$id))
selected.stations <- Central.African.stations %>%
  filter(id %in% selected.ids)

ggplot() +
  geom_raster(data = df_r %>%
                filter(mask == 1),
              aes(x = x, y = y, fill = as.factor(mask)), alpha = 0.5) +
  guides(fill = "none") +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  geom_point(data = selected.stations,
             aes(x = longitude, y = latitude),
             color = "red",
             size = 0.6, alpha = 0.5) +
  coord_sf(xlim = c(-30, 65), ylim = c(-25, 25), expand = FALSE) +
  theme_bw()


saveRDS(selected.stations %>%
          rename(station_id = id,
                 lon = longitude,
                 lat = latitude),
        "./data/station_data/Meteostat/metadata.Meteostat.RDS")

data2save <- all.df.selected %>%
  rename(station_id = id,
         tas = temp,
         tasmin = tmin,
         tasmax = tmax,
         pre = prcp)


saveRDS(data2save,
        "./data/station_data/Meteostat/data.Meteostat.RDS")


SC <- all.df.selected.long %>%
  filter(id %in% selected.ids) %>%
  left_join(selected.stations,
            by = "id") %>%
  mutate(hemisphere = case_when(latitude>=0 ~ "N",
                                TRUE ~ "S")) %>%
  group_by(hemisphere,id,month,variable) %>%
  summarise(value.m = mean(value,na.rm = TRUE),
            .groups = "keep")


SC.m <- SC %>%
  group_by(hemisphere,month,variable) %>%
  summarise(value.m = mean(value.m,na.rm = TRUE),
            .groups = "keep")

ggplot() +
  geom_line(data = SC.m,
            aes(x = month, y = value.m)) +
  theme_bw() +
  facet_grid(variable ~ hemisphere, scales = "free")

SC %>%
  group_by(variable,hemisphere) %>%
  summarise(N = length(unique(id)))



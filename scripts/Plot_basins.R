rm(list = ls())

library(sf)
library(ggplot2)

Basins.main <- st_read("./data/shapefiles/hybas_lake_af_lev01-12_v1c/hybas_lake_af_lev02_v1c.shp")
# "1020018110

ggplot(Basins.main) +
  geom_sf(aes(fill = as.factor(HYBAS_ID)), color = "grey20", linewidth = 0.2) +
  coord_sf(expand = FALSE) +
  scale_fill_viridis_d(name = "HYBAS_ID") +
  theme_minimal()

Basins <- st_read("./data/shapefiles/hybas_lake_af_lev01-12_v1c/hybas_lake_af_lev04_v1c.shp")

selected <- as.character(Basins.main$PFAF_ID)[c(3)]
selected <- c("132","133")

subbasins_big_pfaf <- Basins %>%
  dplyr::filter(
    Reduce(`|`, lapply(selected, function(p)
      startsWith(as.character(PFAF_ID), p)))
  )

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

Mask <- read_sf("./data/shapefiles/Rainforests.dbf")

ggplot(subbasins_big_pfaf) +
  geom_sf(aes(fill = as.factor(HYBAS_ID)), color = "grey20", linewidth = 0.2) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  geom_sf(data = Mask, fill = NA, color = "red", linewidth = 0.2) +
  coord_sf(expand = FALSE) +
  scale_fill_viridis_d(name = "HYBAS_ID") +
  coord_sf(xlim = c(-15, 55), ylim = c(-25, 15), expand = FALSE) +
  theme_bw()

ggplot(Basins) +
  geom_sf(aes(fill = as.factor(HYBAS_ID)), color = "grey20", linewidth = 0.2) +
  coord_sf(expand = FALSE) +
  scale_fill_viridis_d(name = "HYBAS_ID") +
  theme_minimal() +
  guides(fill = "none")


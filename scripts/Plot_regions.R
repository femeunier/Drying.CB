library(terra)
library(sf)
library(ggplot2)
library(rnaturalearth)

CATL <- ext(-23, 3, -2.4, 5.5)
CC   <- ext(22, 29, -10, 4)

# function to convert ext to sf polygon
ext_to_sf <- function(e, name){
  as.polygons(e) |>
    st_as_sf() |>
    st_set_crs(4326) |>
    dplyr::mutate(region = name)
}

catl_sf <- ext_to_sf(CATL, "CATL")
cc_sf   <- ext_to_sf(CC, "CC")

rects <- rbind(catl_sf, cc_sf)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
zoom <- st_bbox(c(xmin = -25, xmax = 55, ymin = -15, ymax = 10), crs = 4326)

ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey40") +
  geom_sf(data = rects, color = "darkgrey", fill = NA, linewidth = 1.2) +
  coord_sf(xlim = c(zoom["xmin"], zoom["xmax"]),
           ylim = c(zoom["ymin"], zoom["ymax"]),
           expand = FALSE) +
  theme_minimal() +
  labs(color = "Region")

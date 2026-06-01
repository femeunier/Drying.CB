rm(list = ls())

library(dplyr)
library(slider)
library(tidyr)
library(scales)
library(ggplot2)

system2("scp",
        c("hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/All.climatevars.Amazon.RDS",
          "./outputs/"))

All.anomalies.Amazon <- readRDS("./outputs/All.anomalies.Amazon.RDS")

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")


ggplot() +
  geom_raster(data = All.anomalies.Amazon %>%
                dplyr::filter(var == "pre",
                              type == "Observational") %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = mean)) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-85, -35), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(limits = c(-1,1)*100,
                       oob = scales::squish) +
  facet_wrap(~ model) +
  theme_bw()


ggplot() +
  geom_raster(data = All.anomalies.Amazon %>%
                dplyr::filter(var == "tas",
                              type == "Observational") %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = mean)) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-85, -35), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(limits = c(-1,1)*2,high = "darkred",low = "darkblue",
                       oob = scales::squish) +
  facet_wrap(~ model) +
  theme_bw()

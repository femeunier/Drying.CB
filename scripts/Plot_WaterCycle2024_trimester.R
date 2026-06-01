rm(list = ls())

library(dplyr)
library(slider)
library(tidyr)
library(scales)
library(raster)
library(ggplot2)
library(corrplot)

system2("scp",
        c("hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/All.anomalies.CA.RDS",
          "./outputs/"))


All.anomalies.CA <- readRDS("./outputs/All.anomalies.CA.RDS")

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  geom_raster(data = All.anomalies.CA %>%
                dplyr::filter(var == "pre",
                              year == 2024,
                              trimester == 3,
                              type == "Observational") %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = value)) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(limits = c(-1,1)*100,
                       oob = scales::squish) +
  facet_wrap(~ model) +
  theme_bw() +
  labs(x = "", y = "", fill = "Anomaly (mm/month)") +
  theme(legend.position = c(0.9,0.1),
        panel.spacing = unit(2,"lines"),
        text = element_text(size = 20))

ggplot(data = All.anomalies.CA %>%
         dplyr::filter(var == "pre",
                       type == "Observational") %>%
         ungroup()) +
  geom_density(aes(x = value,
                   fill = model),
               alpha = 0.5, color = NA) +
  theme_bw() +
  facet_grid(year ~ trimester) +
  # scale_x_continuous(limits = c(-1,1)*100) +
  geom_vline(xintercept = 0,
             linetype = 2) +
  labs(x = "Anomaly (mm/month)", fill = "",
       y = "") +
  theme(legend.position = "right",
        text = element_text(size = 18))



ggplot() +
  geom_raster(data = All.anomalies.CA %>%
                dplyr::filter(var == "tas",
                              year == 2024, trimester == 3,
                              type == "Observational") %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = value)) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(limits = c(-1,2),
                       low = muted("blue"),
                       high = muted("red"),
                       oob = scales::squish) +
  facet_wrap(~ model) +
  theme_bw() +
  labs(x = "", y = "", fill = "Anomaly (°C)") +
  theme(legend.position = c(0.9,0.1),
        panel.spacing = unit(2,"lines"),
        text = element_text(size = 20))

ggplot(data = All.anomalies.CA %>%
         dplyr::filter(var == "tas",
                       type == "Observational") %>%
         ungroup()) +
  geom_density(aes(x = value,
                   fill = model),
               alpha = 0.5, color = NA) +
  theme_bw() +
  scale_x_continuous(limits = c(-1,2.5)) +
  geom_vline(xintercept = 0,
             linetype = 2) +
  facet_grid(year ~ trimester) +
  labs(x = "Anomaly (°C)", fill = "",
       y = "") +
  theme(legend.position = "right",
        text = element_text(size = 18))


#################################################################################

ggplot() +
  geom_raster(data = All.anomalies.CA %>%
                dplyr::filter(var == "SMrz") %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = value)) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(limits = c(-1,1)*0.1,
                       oob = scales::squish) +
  facet_grid(year ~ trimester) +
  theme_bw() +
  labs(x = "", y = "", fill = "Anomaly (mm/month)") +
  theme(legend.position = "right",
        panel.spacing = unit(2,"lines"),
        text = element_text(size = 20))




A <- readRDS("./outputs/All.timeseries.CA.RDS") %>%
  filter(var == "Cbackscatter") %>%
  arrange(model, var, type, scenario, time) %>%
  group_by(model, var, type, scenario) %>%
  mutate(
    roll3_mean = slide_dbl(value, mean, .before = 3, .complete = TRUE, na.rm = TRUE)
  )

plot(A$roll3_mean,type = "l")

ggplot() +
  geom_raster(data = All.anomalies.CA %>%
                dplyr::filter(var == "Cbackscatter") %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = value)) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(limits = c(-1,1)*0.1,
                       oob = scales::squish) +
  facet_grid(year ~ trimester) +
  theme_bw() +
  labs(x = "", y = "", fill = "Anomaly (mm/month)") +
  theme(legend.position = "right",
        panel.spacing = unit(2,"lines"),
        text = element_text(size = 20))


A <- readRDS("./outputs/All.timeseries.CA.RDS") %>%
  filter(var == "tws") %>%
  arrange(model, var, type, scenario, time) %>%
  group_by(model, var, type, scenario) %>%
  mutate(
    roll3_mean = slide_dbl(value, mean, .before = 3, .complete = TRUE, na.rm = TRUE)
  )

plot(A$roll3_mean,type = "l")

ggplot() +
  geom_raster(data = All.anomalies.CA %>%
                dplyr::filter(var == "tws") %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = value)) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(limits = c(-1,1)*40,
                       oob = scales::squish) +
  facet_grid(year ~ trimester) +
  theme_bw() +
  labs(x = "", y = "", fill = "Anomaly (mm/month)") +
  theme(legend.position = "right",
        panel.spacing = unit(2,"lines"),
        text = element_text(size = 20))


A <- readRDS("./outputs/All.timeseries.CA.RDS") %>%
  filter(var == "EVI") %>%
  arrange(model, var, type, scenario, time) %>%
  group_by(model, var, type, scenario) %>%
  mutate(
    roll3_mean = slide_dbl(value, mean, .before = 3, .complete = TRUE, na.rm = TRUE)
  )

plot(A$roll3_mean/10000/10000,type = "l")

ggplot() +
  geom_raster(data = All.anomalies.CA %>%
                dplyr::filter(var == "EVI") %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = value)) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(limits = c(-1,1)*0.1,
                       oob = scales::squish) +
  facet_grid(year ~ trimester) +
  theme_bw() +
  labs(x = "", y = "", fill = "Anomaly") +
  theme(legend.position = "right",
        panel.spacing = unit(2,"lines"),
        text = element_text(size = 20))

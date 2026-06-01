rm(list = ls())

library(dplyr)
library(slider)
library(tidyr)
library(scales)
library(raster)
library(ggplot2)
library(corrplot)

system2("scp",
        c("hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/All.slopes.CA.RDS",
          "./outputs/"))

All.slopes.CA <- readRDS("./outputs/All.slopes.CA.RDS") %>%
  filter(!(model %in% c("NCEP","RFE2","CIESM")))

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")


selected <- All.slopes.CA %>%
  filter(var == "pre",
         type == "CMIP6",
         period == "historical") %>%
  group_by(model) %>%
  summarise(N = length(which(slope_per_year < 0))) %>%
  filter(N > 200)



ggplot() +
  geom_raster(data = All.slopes.CA %>%
                filter(var == "pre",
                       type == "Observational",
                       period == "historical") %>%
                rename(origin = type) %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = slope_per_year*10)) +
  geom_sf(data = world, fill = NA,
          color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(limits = c(-1,1)*25,
                       oob = scales::squish) +
  facet_wrap(~ model) +
  theme_bw() +
  labs(x = "", y = "", fill = "MAP (mm/decade)") +
  theme(legend.position = "right",
        panel.spacing = unit(2,"lines"),
        text = element_text(size = 20))


ggplot() +
  geom_raster(data = All.slopes.CA %>%
                filter(var == "pre",
                       type == "CMIP6",
                       model %in% c(selected$model),
                       period == "historical") %>%
                rename(origin = type) %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = slope_per_year*10)) +
  geom_sf(data = world, fill = NA,
          color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(limits = c(-1,1)*25,
                       oob = scales::squish) +
  facet_wrap(~ model) +
  theme_bw() +
  labs(x = "", y = "", fill = "MAP (mm/decade)") +
  theme(legend.position = "right",
        panel.spacing = unit(2,"lines"),
        text = element_text(size = 20))

ggplot() +
  geom_density(data = All.slopes.CA %>%
                 filter(var == "pre") %>%
                 filter((type == "CMIP6" & period == "historical") |
                          (type == "Observational")) %>%
                 rename(origin = type) %>%
                 ungroup(),
               aes(x = slope_per_year*10,
                   fill = model),
               alpha = 0.5, color = NA) +
  facet_wrap(~ origin) +
  theme_bw()+
  guides(fill = "none")


ggplot() +
  geom_raster(data = All.slopes.CA %>%
                filter(var == "tas",
                       type == "Observational",
                       period == "historical") %>%
                rename(origin = type) %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = intercept_t0)) +
  geom_sf(data = world, fill = NA,
          color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(oob = scales::squish) +
  facet_wrap(~ model) +
  theme_bw() +
  labs(x = "", y = "", fill = "MAT (°C)") +
  theme(legend.position = c(0.9,0.1),
        panel.spacing = unit(2,"lines"),
        text = element_text(size = 20))


ggplot() +
  geom_density(data = All.slopes.CA %>%
                 filter(var == "tas") %>%
                 filter((type == "CMIP6" & period == "historical") |
                          (type == "Observational")) %>%
                 rename(origin = type) %>%
                 ungroup(),
               aes(x = intercept_t0,
                   fill = model),
               alpha = 0.5, color = NA) +
  facet_wrap(~ origin) +
  theme_bw() +
  guides(fill = "none")



################################################################################
# MEM

All.slopes.CA.select <- All.slopes.CA %>%
  filter((period == "historical" & type == "Observational") |
           (period == "historical"),
         var %in% c("tas","pre")) %>%
  group_by(var,type,lon,lat) %>%
  summarise(slope.m = mean(slope_per_year,na.rm = TRUE),
            intercept.m = mean(intercept_t0,na.rm = TRUE),
            .groups = "keep")

################################################################################
# MAP

ggplot() +
  geom_raster(data = All.slopes.CA.select %>%
                filter(var == "pre") %>%
                rename(origin = type) %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = slope.m*10)) +
  geom_sf(data = world, fill = NA,
          color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(limits = c(-1,1)*25,
                       oob = scales::squish) +
  facet_grid( ~ origin) +
  theme_bw() +
  labs(x = "", y = "", fill = "slope (mm/decade)") +
  theme(legend.position = "right",
        panel.spacing = unit(2,"lines"),
        text = element_text(size = 20))


ggplot() +
  geom_density(data = All.slopes.CA.select %>%
                 filter(var == "pre") %>%
                 rename(origin = type) %>%
                 ungroup(),
               aes(x = slope.m*10,
                   fill = origin),
               alpha = 0.5, color = NA) +
  theme_bw()


ggplot() +
  geom_raster(data = All.slopes.CA.select %>%
                filter(var == "pre") %>%
                rename(origin = type) %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = intercept.m*12)) +
  geom_sf(data = world, fill = NA,
          color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(oob = scales::squish) +
  facet_grid( ~ origin) +
  theme_bw() +
  labs(x = "", y = "", fill = "MAP (mm/year)") +
  theme(legend.position = c(0.9,0.1),
        panel.spacing = unit(2,"lines"),
        text = element_text(size = 20))


ggplot() +
  geom_density(data = All.slopes.CA.select %>%
                 filter(var == "pre") %>%
                 rename(origin = type) %>%
                 ungroup(),
               aes(x = intercept.m*12,
                   fill = origin),
               alpha = 0.5, color = NA) +
  theme_bw()

################################################################################
# tas

ggplot() +
  geom_raster(data = All.slopes.CA.select %>%
                filter(var == "tas") %>%
                rename(origin = type) %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = slope.m*10)) +
  geom_sf(data = world, fill = NA,
          color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(limits = c(-0.5,1)*0.5,
                       oob = scales::squish) +
  facet_grid( ~ origin) +
  theme_bw() +
  labs(x = "", y = "", fill = "slope (mm/decade)") +
  theme(legend.position = "right",
        panel.spacing = unit(2,"lines"),
        text = element_text(size = 20))


ggplot() +
  geom_density(data = All.slopes.CA.select %>%
                 filter(var == "tas") %>%
                 rename(origin = type) %>%
                 ungroup(),
               aes(x = slope.m*10,
                   fill = origin),
               alpha = 0.5, color = NA) +
  theme_bw()



ggplot() +
  geom_raster(data = All.slopes.CA.select %>%
                filter(var == "tas") %>%
                rename(origin = type) %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = intercept.m)) +
  geom_sf(data = world, fill = NA,
          color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(oob = scales::squish) +
  facet_grid( ~ origin) +
  theme_bw() +
  labs(x = "", y = "", fill = "MAT (°C)") +
  theme(legend.position = c(0.9,0.1),
        panel.spacing = unit(2,"lines"),
        text = element_text(size = 20))


ggplot() +
  geom_density(data = All.slopes.CA.select %>%
                 filter(var == "tas") %>%
                 rename(origin = type) %>%
                 ungroup(),
               aes(x = intercept.m,
                   fill = origin),
               alpha = 0.5, color = NA) +
  theme_bw()


################################################################################
# RMSE of CMIP6 models

CMIP6.vs.Obs <- All.slopes.CA %>%
  filter(type == "CMIP6",
         period == "historical") %>%
  left_join(All.slopes.CA.select %>%
              ungroup() %>%
              filter(type == "Observational") %>%
              dplyr::select(-type) %>%
              rename(slope_ref = slope.m,
                     intercept_ref = intercept.m),
            by = c("lon","lat","var"))

ggplot(data = CMIP6.vs.Obs %>%
         filter(var == "pre"),
       aes(x = slope_ref*10, y = slope_per_year*10,
           color = model)) +
  geom_point() +
  stat_smooth(method = "lm",
              se = FALSE) +
  geom_abline(slope = 1, intercept = 0,
              linetype = 2) +
  facet_wrap(~ var) +
  theme_bw() +
  guides(color = "none")

RMSE <- CMIP6.vs.Obs %>%
  filter(var == "pre") %>%
  group_by(var,model) %>%
  summarise(rmse = 1/n()*sqrt(sum((slope_ref - slope_per_year)**2,na.rm = TRUE)),
            .groups = "keep") %>%
  arrange(rmse) %>%
  group_by(var) %>%
  mutate(w = (1/rmse)**10) %>%
  mutate(w = w/sum(w))

hist(RMSE$w)


ggplot() +
  geom_raster(data = All.slopes.CA %>%
                filter(var == "pre",
                       type == "CMIP6",
                       model %in% (RMSE %>%
                                     ungroup() %>%
                                     slice_head(n = 5) %>%
                                     pull(model)),
                       period == "historical") %>%
                rename(origin = type) %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = slope_per_year*10)) +
  geom_sf(data = world, fill = NA,
          color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(limits = c(-1,1)*25,
                       oob = scales::squish) +
  facet_wrap(~ model) +
  theme_bw() +
  labs(x = "", y = "", fill = "MAP (mm/decade)") +
  theme(legend.position = c(0.9,0.1),
        panel.spacing = unit(2,"lines"),
        text = element_text(size = 20))


All.slopes.CA.select.comp <- bind_rows(All.slopes.CA.select %>%
                                         filter(type == "Observational",
                                                var == "pre") %>%
                                         rename(origin = type),
                                       All.slopes.CA %>%
  filter((period == "historical"),
         var %in% c("pre")) %>%
  left_join(RMSE,
            by = c("model","var")) %>%
  group_by(var,type,lon,lat) %>%
  summarise(slope.m = weighted.mean(slope_per_year,w = w,
                                    na.rm = TRUE),
            origin = "CMIP6",
            intercept.m = mean(intercept_t0,na.rm = TRUE),
            .groups = "keep"))


ggplot() +
  geom_raster(data = All.slopes.CA.select.comp %>%
                filter(var == "pre") %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = slope.m*10)) +
  geom_sf(data = world, fill = NA,
          color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(limits = c(-1,1)*5,
                       oob = scales::squish) +
  theme_bw() +
  facet_wrap(~origin) +
  labs(x = "", y = "", fill = "MAP (mm/decade)") +
  theme(legend.position = c(0.9,0.1),
        panel.spacing = unit(2,"lines"),
        text = element_text(size = 20))


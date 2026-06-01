rm(list = ls())

library(dplyr)
library(ggplot2)
library(Drying.CB)
library(slider)
library(terra)
library(sf)


system2("scp",
        c("hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/All.timeseries.CA.RDS",
          "./outputs/"))

system2("scp",
        c("hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/All.slopes.CA.RDS",
          "./outputs/"))

CATL <- ext(-23,3,-2.4,5.5)
CC <- ext(22,29,-10,4)

raw <- readRDS("./outputs/All.timeseries.CA.RDS") %>%
  filter(!(model %in% c("RFE2","NCEP")))
A <- bind_rows(raw,
               raw %>%
                 filter(type == "Observational",
                        var == "pre",
                        scenario == "historical") %>%
                 group_by(time,year,month,var,type,scenario) %>%
                 summarise(value = mean(value,na.rm = TRUE),
                           model = "MEM",
                           .groups = "keep")) %>%
  arrange(model, var, type, scenario, time) %>%
  group_by(model, var, type, scenario) %>%
  mutate(
    # 12-month window ending at the current month
    roll12_sum  = slide_dbl(value, sum,  .before = 11, .complete = TRUE, na.rm = TRUE),
    roll12_mean = slide_dbl(value, mean, .before = 11, .complete = TRUE, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(roll12 = case_when(var == "pre" ~ roll12_sum,
                            TRUE ~ roll12_mean)) %>%
  dplyr::select(-c(roll12_sum,roll12_mean)) %>%
  filter(model != "CIESM")

A.sum <- A %>%
  filter(var == "pre",
         model == "MEM") %>%
  filter((type == "CMIP6" & scenario == "ssp245") |
           (type == "Observational")) %>%
  group_by(year,month,type) %>%
  summarise(roll12 = mean(roll12),
            .groups = "keep")

A.ind <-  A %>%
  filter(var == "pre",
         model != "MEM") %>%
  filter((type == "CMIP6" & scenario == "ssp245")) %>%
  group_by(year,month,model,type) %>%
  summarise(roll12 = mean(roll12),
            .groups = "keep")

slopes <- A.ind %>%
  mutate(time = year + (month - 1/2)/12) %>%
  filter(year >= 1980, year <= 2014) %>%
  group_by(model,type) %>%
  summarise(slope = coef(lm(roll12 ~ time))[2]*10,
            signif = summary(lm(roll12 ~ time))[["coefficients"]][2,4],
            .groups = "keep")

selected <- slopes %>%
  filter(slope < 0) %>%
  pull(model)
# selected <- "KACE-1-0-G"
hist(slopes$slope)

ggplot(data = A.sum,
       aes(x = year + (month - 1/2)/12,
           y = roll12,
           color = type)) +
  geom_line(data = A.ind %>%
              filter(model %in% selected),
            aes(group = model),linewidth = 0.25) +
  stat_smooth(data = A.ind %>%
              filter(model %in% selected),
              method = "lm",se = FALSE,
            aes(group = model),linewidth = 0.25) +
  geom_line(linewidth = 1) +
  stat_smooth(method = "lm", linewidth = 1,
              se = FALSE) +
  scale_x_continuous(limits = c(1980,2025)) +
  theme_bw()


ggplot(data = A.sum,
       aes(x = year + (month - 1/2)/12,
           y = roll12,
           color = type)) +
  geom_line(data = A.ind %>%
              filter(model %in% selected),
            aes(group = model),linewidth = 0.25) +
  stat_smooth(data = A.ind %>%
                filter(!(model %in% selected)),
              method = "lm",se = FALSE,
              aes(group = model),linewidth = 0.25) +
  geom_line(linewidth = 1) +
  stat_smooth(method = "lm", linewidth = 1,
              se = FALSE) +
  scale_x_continuous(limits = c(1980,2025)) +
  theme_bw()

################################################################################

raw.slopes <- readRDS("./outputs/All.slopes.CA.RDS") %>%
  filter(!(model %in% c("RFE2","NCEP")))

Slopes.Obs <- bind_rows(raw.slopes  %>%
                          filter(var == "pre",
                                 type == "Observational"),
                        raw.slopes %>%
                          filter(var == "pre",
                                 type == "Observational") %>%
                          group_by(lon,lat,var,type,period) %>%
                          summarise(slope_per_year = mean(slope_per_year,na.rm = TRUE),
                                    model = "MEM",
                                    .groups = "keep")) %>%
  group_by(model) %>%
  mutate(frac.neg = length(which(slope_per_year < 0))/n()*100) %>%
  mutate(frac.neg.cat = case_when(frac.neg > 50 ~ "50+%",
                                  frac.neg > 25 ~ "25-50%",
                                  frac.neg > 10  ~ "10-25%",
                                  frac.neg > 0 ~ "0-10%",
                                  TRUE ~ "None")) %>%
  group_by(frac.neg.cat) %>%
  mutate(frac.neg.cat = paste0(frac.neg.cat,", N = ",
                               length(unique(model)),
                               " products"))


ggplot(data = Slopes.Obs) +
  geom_density(aes(x = slope_per_year,
                   fill = model),
               color = NA,
               alpha = 0.5) +
  facet_wrap(~ frac.neg.cat, nrow = 1) +
  theme_bw() +
  guides(fill = "none")

All.slopes.CA <- raw.slopes %>%
  ungroup() %>%
  dplyr::filter(type == "CMIP6",
                period == "ssp245",
                var == "pre") %>%
  filter(model != "CIESM") %>%
  group_by(model) %>%
  mutate(frac.neg = length(which(slope_per_year < 0))/n()*100) %>%
  mutate(frac.neg.cat = case_when(frac.neg > 50 ~ "50+%",
                                  frac.neg > 25 ~ "25-50%",
                                  frac.neg > 10  ~ "10-25%",
                                  frac.neg > 0 ~ "0-10%",
                                  TRUE ~ "None")) %>%
  group_by(frac.neg.cat) %>%
  mutate(frac.neg.cat = paste0(frac.neg.cat,", N = ",
                               length(unique(model)),
                               " models"))

ggplot(data = All.slopes.CA) +
  geom_density(aes(x = slope_per_year,
                   fill = model),
               color = NA,
               alpha = 0.5) +
  facet_wrap(~ frac.neg.cat, nrow = 1) +
  theme_bw() +
  guides(fill = "none")

selected.models <- All.slopes.CA %>%
  group_by(model) %>%
  summarise(N = length(which(slope_per_year < 0)),
            .groups = "keep") %>%
  dplyr::filter(N > 200) %>%
  pull(model)

selected.models <- All.slopes.CA %>%
  filter(frac.neg > 40) %>%
  pull(model) %>%
  unique()

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  geom_raster(data = Slopes.Obs %>%
                filter(model %in% c("ERA5")),
              aes(x = lon, y = lat,
                  fill = slope_per_year*10)) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2() +
  facet_wrap(~ model) +
  theme_bw() +
  labs(x = "", y = "", fill = "slope (mm/decade)") +
  theme(legend.position = "right",
        panel.spacing = unit(2,"lines"),
        text = element_text(size = 20))

ggplot() +
  geom_raster(data = All.slopes.CA %>%
                filter(model %in% selected.models),
              aes(x = lon, y = lat,
                  fill = slope_per_year*10)) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  geom_sf(data =  vect(CATL, crs = "EPSG:4326") |> st_as_sf(),
          fill = NA, colour = "red", linewidth = 0.7) +
  geom_sf(data =  vect(CC, crs = "EPSG:4326") |> st_as_sf(),
          fill = NA, colour = "red", linewidth = 0.7) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2() +
  facet_wrap(~ model) +
  theme_bw() +
  labs(x = "", y = "", fill = "slope (mm/decade)") +
  theme(legend.position = "right",
        panel.spacing = unit(2,"lines"),
        text = element_text(size = 20))

ggplot() +
  geom_raster(data = Slopes.Obs ,
              aes(x = lon, y = lat,
                  fill = slope_per_year*10)) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2() +
  facet_wrap(~ model) +
  theme_bw() +
  labs(x = "", y = "", fill = "slope (mm/decade)") +
  theme(legend.position = "right",
        panel.spacing = unit(2,"lines"),
        text = element_text(size = 20))



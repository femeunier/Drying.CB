rm(list = ls())

library(dplyr)
library(slider)
library(tidyr)
library(scales)
library(ggplot2)

system2("scp",
        c("hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/All.timeseries.CA.RDS",
        "./outputs/"))

All.timeseries.CA <- readRDS("./outputs/All.timeseries.CA.RDS") %>%
  filter(model != "CIESM")

A <- All.timeseries.CA %>%
  filter(var == "pre",
         year >= 1961,
         scenario  == "historical",
         model == "AWI-ESM-1-REcoM")

plot(A$value,type = "l")

coef(lm(data = A,
           value ~ time))


a <- All.timeseries.CA %>%
  filter(var == "pre",
         type == "Observational") %>%
  group_by(model) %>%
  filter(row_number() <= 1 | row_number() > n() - 1) %>%
  summarise(t.extent =
              paste0(year[year == min(year)],"/",month[year == min(year)],
                     "-",
                     year[year == max(year)],"/",month[year == max(year)]))

A <- bind_rows(All.timeseries.CA,
               All.timeseries.CA %>%
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
  dplyr::select(-c(roll12_sum,roll12_mean))


ggplot(A %>%
         filter(year %in% 1981:2010)) +
  geom_density(aes(x = value, fill = model),
               alpha = 0.5, color = NA) +
  facet_wrap(~ interaction(type,var), scales = "free") +
  theme_bw() +
  guides(fill = "none")

ggplot(A %>%
         filter(year %in% 1981:2010,
                var == "tasmax",
                type == "Observational")) +
  geom_density(aes(x = value, fill = model),
               alpha = 0.5, color = NA) +
  facet_wrap(~ interaction(type,var), scales = "free") +
  theme_bw()

A %>%
  filter(year %in% 1981:2010,
         var == "pre",
         type == "Observational") %>%
  group_by(model) %>%
  summarise(m = mean(value))

ggplot(A %>%
         filter(year %in% 1981:2010)) +
  geom_density(aes(x = value, fill = model),
               alpha = 0.5, color = NA) +
  facet_wrap(~ interaction(type,var), scales = "free") +
  theme_bw() +
  guides(fill = "none")

A.anomaly <- A %>%
  group_by(var,model,type,scenario,month) %>%
  mutate(m.month = mean(value[year %in% c(1981:2010)],na.rm = TRUE),
         sd.month = sd(value[year %in% c(1981:2010)],na.rm = TRUE)) %>%
  group_by(var,model) %>%
  mutate(anomaly = value - m.month,
         anomaly.norm = (value - m.month)/(sd.month)) %>%
  mutate(anomaly.rm =
           slide_dbl(anomaly, mean, .before = 11, .complete = TRUE, na.rm = TRUE),
         anomaly.norm.rm =
           slide_dbl(anomaly.norm, mean, .before = 11, .complete = TRUE, na.rm = TRUE))

cvar = "tas"
ggplot(data = A.anomaly %>%
         filter(year >= 1961,
                type == "Observational",var == cvar),
       aes(x = year + (month - 1/2)/12,
           y = roll12)) +
  geom_line(aes(color = model)) +
  geom_line(data = A.anomaly %>%
              filter(year >= 1961,
                     type == "Observational",var == cvar,
                     model == "MEM") ,
            aes(x = year + (month - 1/2)/12,
                y = roll12),
            color = "black") +
  stat_smooth(data = A.anomaly %>%
                filter(year >= 1961,var == cvar,
                       type == "Observational") ,
              color = "black",
              # method = "loess",
              se = FALSE) +
  facet_wrap(~ interaction(type,
                           var),
             scales = "free") +
  scale_x_continuous(limits = c(1980,2027)) +
  theme_bw()

cvar = "pre"
ggplot(data = A.anomaly %>%
         filter(year >= 1961,
                type == "Observational",var == cvar),
       aes(x = year + (month - 1/2)/12,
           y = roll12)) +
  geom_line(aes(color = model)) +
  geom_line(data = A.anomaly %>%
            filter(year >= 1961,
                   type == "Observational",var == cvar,
                   model == "MEM") ,
            aes(x = year + (month - 1/2)/12,
                y = roll12),
            color = "black") +
  stat_smooth(data = A.anomaly %>%
                filter(year >= 1961,var == cvar,
                       type == "Observational") ,
              color = "black",
              # method = "loess",
              se = FALSE) +
  facet_wrap(~ interaction(type,
                           var),
             scales = "free") +
  scale_x_continuous(limits = c(1980,2027)) +
  theme_bw()

stop()

A.anomaly %>%
  filter(year >= 1961,
         type == "Observational",var == "tasmax", model == "Berk") %>%
  View()


ggplot(data = A.anomaly %>%
         filter(year >= 1961,
                scenario == "historical",
                type == "CMIP6"),
       aes(x = year + (month - 1/2)/12,
           y = roll12)) +
  geom_line(aes(color = model)) +
  geom_line(data = A.anomaly %>%
              filter(year >= 1961,
                     type == "CMIP6",
                     scenario == "historical",
                     model == "MEM") ,
            aes(x = year + (month - 1/2)/12,
                y = roll12),
            color = "black") +
  stat_smooth(data = A.anomaly %>%
                filter(year >= 1961,
                       model == "MEM",
                       scenario == "historical",
                       type == "CMIP6") ,
              color = "black", method = "lm", se = FALSE) +
  facet_wrap(~ interaction(type,
                           var),
             scales = "free") +
  scale_x_continuous(limits = c(1961,2015)) +
  theme_bw() +
  guides(color = "none")


A.anomaly %>%
  filter(year >= 1961,
         type == "CMIP6",
         var == "pre",
         scenario == "historical") %>%
  group_by(model) %>%
  summarise(slope = coef(lm(roll12 ~ time))[2])

ggplot(data = A.anomaly %>%
         filter(year >= 1961,
                scenario == "historical",
                model == "E3SM-2-0",
                type == "CMIP6"),
       aes(x = year + (month - 1/2)/12,
           y = roll12)) +
  geom_line(aes(color = model)) +
  geom_line(data = A.anomaly %>%
              filter(year >= 1961,
                     type == "CMIP6",
                     model == "E3SM-2-0",
                     scenario == "historical",
                     model == "MEM") ,
            aes(x = year + (month - 1/2)/12,
                y = value),
            color = "black") +
  stat_smooth(data = A.anomaly %>%
                filter(year >= 1961,
                       model == "E3SM-2-0",
                       scenario == "historical",
                       type == "CMIP6") ,
              color = "black", method = "lm", se = FALSE) +
  facet_wrap(~ interaction(type,
                           var),
             scales = "free") +
  scale_x_continuous(limits = c(1961,2014)) +
  theme_bw() +
  guides(color = "none")

ggplot(data = A.anomaly %>%
         filter(year >= 1961,
                type == "Observational"),
       aes(x = year + (month - 1/2)/12,
           y = roll12)) +
  geom_line(aes(color = model)) +
  geom_line(data = A.anomaly %>%
              filter(year >= 1961,
                     type == "Observational",
                     model == "MEM") ,
            aes(x = year + (month - 1/2)/12,
                y = roll12),
            color = "black") +
  facet_wrap(~ interaction(type,
                           var),
             scales = "free") +
  scale_x_continuous(limits = c(2020,2025)) +
  theme_bw() +
  guides(color = "none")




ggplot(data = A.anomaly %>%
         filter(year >= 1961,
                var == "pre") %>%
         filter((type == "Observational" & scenario == "historical") |
                  (scenario == "ssp245")),
       aes(x = year + (month - 1/2)/12,
           y = roll12,
           color = model)) +
  geom_line() +
  stat_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ type, scales = "free") +
  scale_x_continuous(limits = c(2000,2025)) +
  theme_bw() +
  guides(color = "none")


stop()

ggplot(data = A.anomaly %>%
         filter(year >= 1961,
                var == "pre") %>%
         filter((type == "Observational" & scenario == "historical") |
                  (scenario == "ssp245")),
       aes(x = year + (month - 1/2)/12,
           y = roll12,
           color = model)) +
  geom_line() +
  stat_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ type, scales = "free") +
  scale_x_continuous(limits = c(2000,2025)) +
  theme_bw()




ggplot(data = A.anomaly %>%
         filter(year >= 1961,
                var == "pre",
                type == "CMIP6"),
       aes(x = year + (month - 1/2)/12,
           y = roll12,
           color = model)) +
  geom_line() +
  stat_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ scenario,
             scales = "free") +
  scale_x_continuous(limits = c(1960,2100)) +
  theme_bw() +
  guides(color = "none")


A.anomaly %>%
  filter(model == "MEM",
         year >= 1961) %>%
  mutate(date = year + (month - 1/2)/12) %>%
  group_by(type,var) %>%
  summarise(slope = coef(lm(roll12 ~ date))[2]*10, # per decade
            p.val = summary(lm(roll12 ~ date))[["coefficients"]][2,4],
            r.sq = summary(lm(roll12 ~ date))[["r.squared"]],
            .groups = "keep")

ggplot(data = A.anomaly) +
  geom_rect(data = data.frame(x = 2024,xend = 2025,
                              y = -Inf,yend = Inf),
            aes(xmin = x, xmax = xend,
                ymin = y, ymax = yend),
            fill = "grey",
            alpha = 0.5) +
  geom_line(aes(x = year + (month - 1/2)/12,
                y = roll12,
                color = interaction(scenario,model))) +
  geom_line(data = A.anomaly %>%
              filter(model == "MEM"),
            aes(x = year + (month - 1/2)/12,
                y = roll12,
                group = scenario),
            color = "black") +
  facet_wrap(~ interaction(type,var),
             scales = "free") +
  scale_x_continuous(limits = c(2000,2025)) +
  theme_bw() +
  guides(color = "none")

ggplot(data = A.anomaly) +
  geom_rect(data = data.frame(x = 2024,xend = 2025,
                              y = -Inf,yend = Inf),
            aes(xmin = x, xmax = xend,
                ymin = y, ymax = yend),
            fill = "grey",
            alpha = 0.5) +
  geom_line(aes(x = year + (month - 1/2)/12,
                y = anomaly,
                color = interaction(scenario,model))) +
  geom_line(data = A.anomaly %>%
              filter(model == "MEM"),
            aes(x = year + (month - 1/2)/12,
                y = anomaly,
                group = scenario),
            color = "black") +
  facet_wrap(~ interaction(type,var),
             scales = "free") +
  scale_x_continuous(limits = c(2020,2025)) +
  theme_bw() +
  geom_hline(yintercept = 0,
             linetype = 2) +
  guides(color = "none")


ggplot(data = A.anomaly) +
  geom_rect(data = data.frame(x = 2024,xend = 2025,
                              y = -Inf,yend = Inf),
            aes(xmin = x, xmax = xend,
                ymin = y, ymax = yend),
            fill = "grey",
            alpha = 0.5) +
  geom_line(aes(x = year + (month - 1/2)/12,
                y = anomaly.norm,
                color = interaction(scenario,model))) +
  geom_line(data = A.anomaly %>%
              filter(model == "MEM"),
            aes(x = year + (month - 1/2)/12,
                y = anomaly.norm,
                group = scenario),
            color = "black") +
  facet_wrap(~ interaction(type,var),
             scales = "free") +
  scale_x_continuous(limits = c(2020,2025)) +
  theme_bw() +
  geom_hline(yintercept = 0,
             linetype = 2) +
  guides(color = "none")

################################################################

system2("scp",
        c("hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/All.slopes.CA.RDS",
          "./outputs/"))

All.slopes.CA <- readRDS("./outputs/All.slopes.CA.RDS")

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  geom_raster(data = All.slopes.CA %>%
               dplyr::filter(var == "pre",
                             type == "Observational") %>%
               ungroup(),
             aes(x = lon, y = lat,
                 fill = slope_per_year*10)) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(limits = c(-1,1)*25,
                       oob = scales::squish) +
  facet_wrap(~ model) +
  theme_bw()

ggplot(data = All.slopes.CA %>%
         dplyr::filter(var == "pre",
                       type == "Observational")) +
  geom_density(aes(x = slope_per_year*10,
                   fill = model),
               alpha = 0.5, color = NA) +
  theme_bw() +
  scale_x_continuous(limits = c(-1,1)*50) +
  geom_vline(xintercept = 0,
             linetype = 2)

ggplot() +
  geom_raster(data = All.slopes.CA %>%
                dplyr::filter(var == "tas",
                              type == "Observational") %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = slope_per_year*10)) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(limits = c(-1,1)*0.5,
                       oob = scales::squish) +
  facet_wrap(~ model) +
  theme_bw()


ggplot(data = All.slopes.CA %>%
         dplyr::filter(var == "tas",
                       type == "Observational")) +
  geom_density(aes(x = slope_per_year*10,
                   fill = model),
               alpha = 0.5, color = NA) +
  theme_bw() +
  geom_vline(xintercept = 0,
             linetype = 2)

ggplot() +
  geom_raster(data = All.slopes.CA %>%
                dplyr::filter(var == "pre",
                              model != "CIESM",
                              period == "ssp585",
                              type == "CMIP6") %>%
                dplyr::filter(grepl("EC-Earth",model)) %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = slope_per_year*10)) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(limits = c(-1,1)*10,
                       oob = scales::squish) +
  facet_wrap(~ model) +
  theme_bw()

ggplot(data = All.slopes.CA %>%
         dplyr::filter(var == "pre",
                       model != "CIESM",
                       period == "historical",
                       type == "CMIP6")) +
  geom_density(aes(x = slope_per_year*10,
                   fill = model),
               alpha = 0.5, color = NA) +
  theme_bw() +
  geom_vline(xintercept = 0,
             linetype = 2) +
  guides(fill = "none")

df.slope.comp <- All.slopes.CA %>%
  dplyr::filter(var == "pre",
                model != "CIESM",
                type == "CMIP6") %>%
  pivot_wider(names_from = period,
              values_from = slope_per_year) %>%
  pivot_longer(cols = starts_with("ssp"),
               names_to = "scenario",
               values_to = "slope")

ggplot(data = df.slope.comp,
       aes(x = historical*10,
           y = slope*10,
           color = model,
           shape = scenario)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  stat_smooth(method = "lm", se = FALSE) +
  geom_abline(slope = 1, intercept = 0,
              color = "black",linetype = 2) +
  theme_bw()

################################################################

system2("scp",
        c("hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/All.anomalies.CA.RDS",
          "./outputs/"))

All.anomalies.CA <- readRDS("./outputs/All.anomalies.CA.RDS")

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

ggplot(data = A.anomaly %>%
         filter(var == "pre",
                type == "Observational"),
       aes(x = year + (month - 1/2)/12,
           y = roll12,
           group = model)) +
  geom_line() +
  stat_smooth(method = "lm", se = FALSE, color = "darkgrey") +
  scale_x_continuous(limits = c(1961,2025)) +
  facet_wrap(~ model) +
  theme_bw()


ggplot(data = A.anomaly %>%
         filter(var == "pre",
                type == "Observational"),
       aes(x = year + (month - 1/2)/12,
           y = anomaly,
           group = model)) +
  geom_hline(yintercept = 0,linetype = 2) +
  geom_line() +
  stat_smooth(method = "lm", se = FALSE, color = "darkgrey") +
  scale_x_continuous(limits = c(2020,2025)) +
  facet_wrap(~ model) +
  theme_bw()


ggplot() +
  geom_raster(data = All.anomalies.CA %>%
                dplyr::filter(var == "pre",
                              type == "Observational") %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = mean)) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(limits = c(-1,1)*100,
                       oob = scales::squish) +
  facet_wrap(~ model) +
  theme_bw()

All.anomalies.CA %>%
  dplyr::filter(var == "pre",
                type == "Observational",
                model == "NCEP") %>%
  pull(mean) %>% hist()

ggplot(data = All.anomalies.CA %>%
         dplyr::filter(var == "pre",
                       type == "Observational")) +
  geom_density(aes(x = mean,
                   fill = model),
               alpha = 0.5, color = NA) +
  theme_bw() +
  geom_vline(xintercept = 0,
             linetype = 2)


ggplot() +
  geom_raster(data = All.anomalies.CA %>%
                dplyr::filter(var == "tas",
                              type == "Observational") %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = mean)) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(
    limits = c(-1/3,1)*3,
                       low = muted("blue"),
                       high = muted("red"),
                       oob = scales::squish) +
  facet_wrap(~ model) +
  theme_bw()

ggplot(data = All.anomalies.CA %>%
         dplyr::filter(var == "tas",
                       type == "Observational")) +
  geom_density(aes(x = mean,
                   fill = model),
               alpha = 0.5, color = NA) +
  theme_bw() +
  geom_vline(xintercept = 0,
             linetype = 2)


ggplot() +
  geom_raster(data = All.anomalies.CA %>%
                dplyr::filter(var == "GPP") %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = mean)) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(oob = scales::squish) +
  facet_wrap(~ model) +
  theme_bw()

ggplot() +
  geom_raster(data = All.anomalies.CA %>%
                dplyr::filter(var == "Cbackscatter") %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = mean)) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(oob = scales::squish) +
  facet_wrap(~ model) +
  theme_bw()

ggplot() +
  geom_raster(data = All.anomalies.CA %>%
                dplyr::filter(var %in%
                                c("NDVI","EVI")) %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = mean)) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(oob = scales::squish,
                       limits = c(-1,1)*0.5e7) +
  facet_wrap(~ var) +
  theme_bw()


ggplot() +
  geom_raster(data = All.anomalies.CA %>%
                dplyr::filter(var %in%
                                c("tws")) %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = mean)) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(oob = scales::squish) +
  facet_wrap(~ var) +
  theme_bw()


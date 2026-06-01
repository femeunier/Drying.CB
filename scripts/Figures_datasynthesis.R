rm(list = ls())

library(dplyr)
library(slider)
library(tidyr)
library(scales)
library(raster)
library(ggplot2)
library(corrplot)

system2("scp",
        c("hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/All.timeseries.CA.RDS",
          "./outputs/"))

system2("scp",
        c("hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/All.slopes.CA.RDS",
          "./outputs/"))

system2("scp",
        c("hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/All.anomalies.CA.RDS",
          "./outputs/"))

system2("scp",
        c("hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/All.Zanomalies.CA.RDS",
          "./outputs/"))

system2("scp",
        c("hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/All.slopesZanomalies.CA.RDS",
          "./outputs/"))

system2("scp",
        c("hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/All.slopesanomalies.CA.RDS",
          "./outputs/"))


All.anomalies.CA <- readRDS("./outputs/All.anomalies.CA.RDS")
All.Zanomalies.CA <- readRDS("./outputs/All.Zanomalies.CA.RDS")
All.slopes.CA <- readRDS("./outputs/All.slopes.CA.RDS")
All.slopesZanomalies.CA <- readRDS("./outputs/All.slopesZanomalies.CA.RDS")
All.slopesanomalies.CA <- readRDS("./outputs/All.slopesanomalies.CA.RDS")




A <- bind_rows(readRDS("./outputs/All.timeseries.CA.RDS"),
               readRDS("./outputs/All.timeseries.CA.RDS") %>%
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


A.anomaly <- A %>%
  group_by(var,model,type,scenario,month) %>%
  mutate(m.month = mean(value[year %in% c(2000:2024)],na.rm = TRUE),
         sd.month = sd(value[year %in% c(2000:2024)],na.rm = TRUE)) %>%
  group_by(var,model) %>%
  mutate(anomaly = value - m.month,
         anomaly.norm = (value - m.month)/(sd.month)) %>%
  mutate(anomaly.rm =
           slide_dbl(anomaly, mean, .before = 11, .complete = TRUE, na.rm = TRUE),
         anomaly.norm.rm =
           slide_dbl(anomaly.norm, mean, .before = 11, .complete = TRUE, na.rm = TRUE))

##############################################################################################"


ggplot(data = A.anomaly %>%
         filter(year >= 2000,
                var %in% c("pre","tas"),
                type == "Observational"),
       aes(x = year + (month - 1/2)/12,
           y = roll12)) +
  geom_line(aes(color = model)) +
  geom_line(data = A.anomaly %>%
              filter(year >= 2000,
                     var %in% c("pre","tas"),
                     type == "Observational",
                     model == "MEM") ,
            aes(x = year + (month - 1/2)/12,
                y = roll12),
            color = "black") +
  stat_smooth(data = A.anomaly %>%
                filter(year >= 2000,
                       var %in% c("pre","tas"),
                       type == "Observational"),
              aes(color = model),
              method = "lm", se = FALSE) +
  stat_smooth(data = A.anomaly %>%
                filter(year >= 1961,
                       var %in% c("pre","tas"),
                       model == "MEM",
                       type == "Observational"),
              color = "black",
              method = "lm", se = FALSE) +
  facet_wrap(~ interaction(var),
             scales = "free") +
  scale_x_continuous(limits = c(2000,2025)) +
  theme_bw() +
  guides(color = "none") +
  theme(text = element_text(size = 18),
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(5,"lines")) +
  labs(x = "",y = "")

A.anomaly %>%
  filter(year >= 2000,
         var %in% c("pre","tas"),
         type == "Observational") %>%
  mutate(time = year + (month - 1/2)/12) %>%
  group_by(model,var) %>%
  summarise(slope = coef(lm(roll12 ~ time))[2]*10,
            .groups = "keep") %>%
  group_by(var) %>%
  summarise(Min = min(slope),
            Max = max(slope),
            MEM = slope[model == "MEM"])


world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

#############################################################################################################
# Pre

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
  theme_bw() +
  labs(x = "", y = "", fill = "slope (mm/decade)") +
  theme(legend.position = c(0.9,0.1),
        panel.spacing = unit(2,"lines"),
        text = element_text(size = 20))

ggplot(data = All.slopes.CA %>%
         dplyr::filter(var == "pre",
                       type == "Observational")) +
  geom_density(aes(x = slope_per_year*10,
                   fill = model),
               alpha = 0.5, color = NA) +
  theme_bw() +
  scale_x_continuous(limits = c(-1,1)*50) +
  geom_vline(xintercept = 0,
             linetype = 2) +
  labs(x = "slope (mm/decade)", fill = "",
       y = "") +
  theme(legend.position = c(0.1,0.75),
        text = element_text(size = 18))

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
  theme_bw() +
  labs(x = "", y = "", fill = "Anomaly (mm/month)") +
  theme(legend.position = c(0.9,0.1),
        panel.spacing = unit(2,"lines"),
        text = element_text(size = 20))

ggplot(data = All.anomalies.CA %>%
         dplyr::filter(var == "pre",
                       type == "Observational") %>%
         ungroup()) +
  geom_density(aes(x = mean,
                   fill = model),
               alpha = 0.5, color = NA) +
  theme_bw() +
  # scale_x_continuous(limits = c(-1,1)*100) +
  geom_vline(xintercept = 0,
             linetype = 2) +
  labs(x = "Anomaly (mm/month)", fill = "",
       y = "") +
  theme(legend.position = c(0.7,0.65),
        text = element_text(size = 18))



#############################################################################################################
# tas

ggplot() +
  geom_raster(data = All.slopes.CA %>%
                dplyr::filter(var == "tas",
                              type == "Observational") %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = slope_per_year*10)) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(limits =  c(-1,1),
                       low = muted("blue"),
                       high = muted("red"),
                       oob = scales::squish) +
  facet_wrap(~ model) +
  theme_bw() +
  labs(x = "", y = "", fill = "slope (°C/decade)") +
  theme(legend.position = c(0.9,0.1),
        panel.spacing = unit(2,"lines"),
        text = element_text(size = 20))

ggplot(data = All.slopes.CA %>%
         dplyr::filter(var == "tas",
                       type == "Observational")) +
  geom_density(aes(x = slope_per_year*10,
                   fill = model),
               alpha = 0.5, color = NA) +
  theme_bw() +
  scale_x_continuous(limits = c(-0.5,1)*2) +
  geom_vline(xintercept = 0,
             linetype = 2) +
  labs(x = "slope (°C/decade)", fill = "",
       y = "") +
  theme(legend.position = c(0.1,0.75),
        text = element_text(size = 18))


ggplot() +
  geom_raster(data = All.anomalies.CA %>%
                dplyr::filter(var == "tas",
                              type == "Observational") %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = mean)) +
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
  geom_density(aes(x = mean,
                   fill = model),
               alpha = 0.5, color = NA) +
  theme_bw() +
  scale_x_continuous(limits = c(-1,2.5)) +
  geom_vline(xintercept = 0,
             linetype = 2) +
  labs(x = "Anomaly (°C)", fill = "",
       y = "") +
  theme(legend.position = c(0.1,0.65),
        text = element_text(size = 18))

################################################################
# Vegetation metrics variables

ggplot(data = A.anomaly %>%
         filter(year >= 2000,
                var %in% c("GPP","EVI","NDVI","Cbackscatter"),
                type == "Observational"),
       aes(x = year + (month - 1/2)/12,
           y = roll12)) +
  geom_line(aes(color = model)) +
  geom_line(data = A.anomaly %>%
              filter(year >= 2000,
                     var %in% c("GPP","EVI","NDVI",'Cbackscatter'),
                     type == "Observational",
                     model == "MEM") ,
            aes(x = year + (month - 1/2)/12,
                y = roll12),
            color = "black") +
  stat_smooth(data = A.anomaly %>%
                filter(year >= 2000,
                       var %in% c("GPP","EVI","NDVI","Cbackscatter"),
                       type == "Observational") ,
              color = "black", method = "lm", se = FALSE) +
  facet_wrap(~ interaction(var),
             scales = "free") +
  scale_x_continuous(limits = c(2000,2025)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        strip.text = element_blank(),
        legend.position = c(0.1,0.85),
        strip.background = element_blank(),
        panel.spacing = unit(3.5,"lines")) +
  labs(x = "",y = "", color = "") +
  guides(color = "none")


ggplot() +
  geom_raster(data = All.slopes.CA %>%
                dplyr::filter(var %in% c("GPP","EVI","NDVI","Cbackscatter"),
                              type == "Observational") %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = rel_slope*100*10)) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(limits = c(-1,1)*5,
                       oob = scales::squish) +
  facet_wrap(~ interaction(model,var)) +
  theme_bw() +
  labs(x = "", y = "", fill = "slope (%/decade)") +
  theme(legend.position = "right",
        panel.spacing = unit(2,"lines"),
        text = element_text(size = 20))

ggplot(data = All.slopes.CA %>%
         dplyr::filter(var %in% c("GPP","EVI","NDVI","Cbackscatter"),
                       type == "Observational")) +
  geom_density(aes(x = rel_slope*100*10,
                   fill = interaction(model,var)),
               alpha = 0.5, color = NA) +
  theme_bw() +
  scale_x_continuous(limits = c(-1,1)*10) +
  geom_vline(xintercept = 0,
             linetype = 2) +
  labs(x = "slope (%/decade)", fill = "",
       y = "") +
  theme(legend.position = c(0.15,0.75),
        text = element_text(size = 18))

ggplot() +
  geom_raster(data = All.Zanomalies.CA %>%
                dplyr::filter(var %in% c("GPP","EVI","NDVI","Cbackscatter"),
                              type == "Observational") %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = mean)) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(oob = scales::squish,
                       limits = c(-2,2)) +
  facet_wrap(~ interaction(model,var)) +
  theme_bw() +
  labs(x = "", y = "", fill = "Z-score (-)") +
  theme(legend.position = "right",
        panel.spacing = unit(2,"lines"),
        text = element_text(size = 20))

ggplot(data = All.Zanomalies.CA %>%
         dplyr::filter(var %in% c("GPP","EVI","NDVI","Cbackscatter"),
                       type == "Observational") %>%
         ungroup()) +
  geom_density(aes(x = mean,
                   fill = interaction(model,var)),
               alpha = 0.5, color = NA) +
  theme_bw() +
  scale_x_continuous(limits = c(-2,3)) +
  geom_vline(xintercept = 0,
             linetype = 2) +
  labs(x = "Z-score (-)", fill = "",
       y = "") +
  theme(legend.position = c(0.7,0.8),
        text = element_text(size = 18))

################################################################
# GLEAM variables


ggplot(data = A.anomaly %>%
         filter(year >= 2000,
                model %in% c("GLEAM","GRACE"),
                !(var %in% c("leakage","std","model")),
                type == "Observational"),
       aes(x = year + (month - 1/2)/12,
           y = roll12)) +
  geom_line(data = A.anomaly %>%
              filter(year >= 2000,
                     model %in% c("GLEAM","GRACE"),
                     !(var %in% c("leakage","std","model")),
                     type == "Observational"),
            aes(x = year + (month - 1/2)/12,
                y = value),
            size = 0.1) +
  geom_line() +
  geom_line(data = A.anomaly %>%
              filter(year >= 2000,
                     model %in% c("GLEAM","GRACE"),
                     !(var %in% c("leakage","std","model")),
                     type == "Observational",
                     model == "MEM") ,
            aes(x = year + (month - 1/2)/12,
                y = roll12),
            color = "black") +
  stat_smooth(data = A.anomaly %>%
                filter(year >= 2000,
                       model %in% c("GLEAM","GRACE"),
                       !(var %in% c("leakage","std","model")),
                       type == "Observational") ,
              color = "black", method = "lm", se = FALSE) +
  facet_wrap(~ interaction(var),
             scales = "free") +
  scale_x_continuous(limits = c(2022,2025)) +
  theme_bw() +
  # theme(text = element_text(size = 18),
  #       strip.text = element_blank(),
  #       legend.position = c(0.1,0.85),
  #       strip.background = element_blank(),
  #       panel.spacing = unit(5,"lines")) +
  labs(x = "",y = "", color = "") +
  guides(color = "none")


ggplot(data = A.anomaly %>%
         filter(year >= 2000,
                model %in% c("GLEAM","GRACE"),
                !(var %in% c("leakage","std","model")),
                type == "Observational"),
       aes(x = year + (month - 1/2)/12,
           y = roll12)) +
  geom_line() +
  geom_line(data = A.anomaly %>%
              filter(year >= 2000,
                     model %in% c("GLEAM","GRACE"),
                     !(var %in% c("leakage","std","model")),
                     type == "Observational",
                     model == "MEM") ,
            aes(x = year + (month - 1/2)/12,
                y = roll12),
            color = "black") +
  stat_smooth(data = A.anomaly %>%
                filter(year >= 2000,
                       model %in% c("GLEAM","GRACE"),
                       !(var %in% c("leakage","std","model")),
                       type == "Observational") ,
              color = "black", method = "lm", se = FALSE) +
  facet_wrap(~ interaction(var),
             scales = "free") +
  scale_x_continuous(limits = c(2000,2025)) +
  theme_bw() +
  # theme(text = element_text(size = 18),
  #       strip.text = element_blank(),
  #       legend.position = c(0.1,0.85),
  #       strip.background = element_blank(),
  #       panel.spacing = unit(5,"lines")) +
  labs(x = "",y = "", color = "") +
  guides(color = "none")


ggplot() +
  geom_raster(data = All.slopes.CA %>%
                dplyr::filter(model %in% c("GLEAM","GRACE"),
                              !(var %in% c("leakage","std","model")),
                              type == "Observational") %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = rel_slope*100*10)) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(limits = c(-1,1)*10,
                       oob = scales::squish) +
  facet_wrap(~ interaction(model,var)) +
  theme_bw() +
  labs(x = "", y = "", fill = "slope (%/decade)") +
  theme(legend.position = c(0.9,0.1),
        panel.spacing = unit(2,"lines"),
        text = element_text(size = 20))

ggplot(data = All.slopes.CA %>%
         dplyr::filter(model %in% c("GLEAM","GRACE"),
                       !(var %in% c("leakage","std","model")),
                       type == "Observational")) +
  geom_density(aes(x = rel_slope*100*10,
                   fill = interaction(model,var)),
               alpha = 0.5, color = NA) +
  theme_bw() +
  scale_x_continuous(limits = c(-1,1)*10) +
  geom_vline(xintercept = 0,
             linetype = 2) +
  labs(x = "slope (%/decade)", fill = "",
       y = "") +
  theme(legend.position = c(0.1,0.75),
        text = element_text(size = 18))

ggplot() +
  geom_raster(data = All.Zanomalies.CA %>%
                dplyr::filter(model %in% c("GLEAM","GRACE"),
                              !(var %in% c("leakage","std","model")),
                              type == "Observational") %>%
                ungroup(),
              aes(x = lon, y = lat,
                  fill = mean)) +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  scale_fill_gradient2(oob = scales::squish,
                       limits = c(-1,1)*5) +
  facet_wrap(~ interaction(model,var)) +
  theme_bw() +
  labs(x = "", y = "", fill = "Z-score (-)") +
  theme(legend.position = c(0.9,0.1),
        panel.spacing = unit(2,"lines"),
        text = element_text(size = 20))

ggplot(data = All.Zanomalies.CA %>%
         dplyr::filter(model == "GLEAM",
                       type == "Observational") %>%
         ungroup()) +
  geom_density(aes(x = mean,
                   fill = interaction(model,var)),
               alpha = 0.5, color = NA) +
  theme_bw() +
  scale_x_continuous(limits = c(-1,1)*20) +
  geom_vline(xintercept = 0,
             linetype = 2) +
  labs(x = "Z-score (-)", fill = "",
       y = "") +
  theme(legend.position = c(0.8,0.85),
        text = element_text(size = 18))


# We resample to get a correlation plot

all.products <- All.Zanomalies.CA %>%
  filter(!(var %in% c("std","model","leakage"))) %>%
  ungroup() %>%
  dplyr::select(model,period,var,type) %>%
  distinct()

coord <- expand.grid(lon = seq(-179.75,179.75,0.5),
                     lat = seq(-30.25,30.25,0.5)) %>%
  mutate(value = 1)
craster <- raster::rasterFromXYZ(coord)

all.df <- data.frame()
for (iproduct in seq(1,nrow(all.products))){

  print(iproduct/nrow(all.products))

  cmodel <- all.products$model[iproduct]
  cvar <- all.products$var[iproduct]
  cperiod <- all.products$period[iproduct]
  ctype <- all.products$type[iproduct]


  cdf <- All.Zanomalies.CA %>%
    filter(model == cmodel,
           var == cvar,
           period == cperiod,
           type == ctype)

  cr <- rasterFromXYZ(cdf %>%
                        dplyr::select(lon,lat,mean))
  cr.rspld <- resample(cr,craster)

  all.df <- bind_rows(all.df,
                      as.data.frame(cr.rspld,
                                    xy = TRUE) %>%
                        rename(lon = x,
                               lat = y) %>%
                        mutate(model = cmodel,
                               var = cvar,
                               period = cperiod,
                               type = ctype))
}

head(all.df)


dat <- bind_rows(all.df %>%
                   filter(!(var %in% c("tas","pre","tasmin","tasmax"))),
                 all.df %>%
                   filter(var %in% c("tas","pre","tasmin","tasmax")) %>%
                   group_by(lon,lat,var,period,type) %>%
                   summarise(mean = mean(mean),
                             model = "MEM",
                             .groups = "keep")) %>%
  dplyr::filter(var %in%
                  c("tas",'pre','GPP',"Cbackscatter",
                    "Ep","S","SMrz","tws")) %>%
  mutate(series = paste(var, model, sep = " | ")) %>%
  dplyr::select(lon, lat, series, mean) %>%
  tidyr::pivot_wider(id_cols = c(lon, lat), names_from = series, values_from = mean) %>%
  dplyr::select(-lon, -lat)


rc <- Hmisc::rcorr(as.matrix(dat), type = "pearson")  # or "spearman"
corrplot(rc$r,
         p.mat = rc$P,
         sig.level = 0.05,
         insig = "blank",
         method = "color",
         type   = "upper",
         order  = "hclust",
         addCoef.col   = "black",  # <- add Pearson r as text
         number.digits = 2,        # digits to show
         number.cex    = 0.7,      # text size
         tl.col = "black",
         tl.cex = 0.6,
         diag   = FALSE)


ggplot(data = dat,
       aes(y = `GPP | FLUXSAT`,
           x = `tas | MEM`)) +
  geom_point(size = 0.1) +
  stat_smooth(method = "lm") +
  theme_bw()

ggplot(data = dat,
       aes(x = `tws | GRACE`,
           y = `tas | MEM`)) +
  geom_point(size = 0.1) +
  stat_smooth(method = "lm") +
  theme_bw()

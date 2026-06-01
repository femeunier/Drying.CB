# rm(list = ls())
#
# library(terra)
# library(lubridate)
# library(dplyr)
#
# files <- list.files("/data/gent/vo/000/gvo00074/ED_common_data/met/Tropics/",
#                     pattern = "ERA5_SST_monthly*",full.names = TRUE)
# sst <- rast(files)
# time(sst) <- as.Date(depth(sst)/86400,origin = "1970-01-01")
#
# CATL <- ext(-23,3,
#             -2.4,5.5)
# sst.sel <- crop(sst,CATL)
#
# raw.data <- data.frame(time = time(sst),
#                                  var = "sst",
#                                  value = as.numeric(
#                                    global(sst.sel,mean,na.rm = TRUE)[["mean"]]))
# saveRDS(raw.data,
#         "./outputs/SST_ERA5.RDS")
#
# library(terra)
# library(lubridate)
# library(dplyr)
# library(tidyr)
#
# files <- list.files("/data/gent/vo/000/gvo00074/ED_common_data/met/Tropics/",
#                     pattern = "ERA5_omega_monthly*",full.names = TRUE)
# Omega <- rast(files)
# time(Omega) <- as.Date(depth(Omega)/86400,origin = "1970-01-01")
#
# CATL <- ext(-23,3,-2.4,5.5)
# CC <- ext(22,29,-10,4)
#
# Omega.CATL <- crop(Omega,CATL)
# Omega.CC <- crop(Omega,CC)
#
# raw.data <- bind_rows(data.frame(time = time(Omega),
#                        var = "ws",
#                        region = "CATL",
#                        value = as.numeric(
#                                    global(Omega.CATL,
#                                           mean,na.rm = TRUE)[["mean"]])),
#                       data.frame(time = time(Omega),
#                                  var = "ws",
#                                  region = "CC",
#                                  value = as.numeric(
#                                    global(Omega.CC,
#                                           mean,na.rm = TRUE)[["mean"]]))) %>%
#   pivot_wider(names_from = "region",
#               values_from = "value")
#
# saveRDS(raw.data,
#         "./outputs/Omega_ERA5.RDS")


#
# library(terra)
# library(lubridate)
# library(dplyr)
# library(tidyr)
#
# files <- "/data/gent/vo/000/gvo00074/felicien/R/outputs/all.climate/ERA5_pre_all.years.tif"
# precip <- rast(files)
#
# CATL <- ext(-23,3,-2.4,5.5)
# CC <- ext(22,29,-10,4)
#
# precip.CATL <- crop(precip,CATL)
# precip.CC <- crop(precip,CC)
#
# raw.data <- bind_rows(data.frame(time = time(precip),
#                        var = "precip",
#                        region = "CATL",
#                        value = as.numeric(
#                                    global(precip.CATL,
#                                           mean,na.rm = TRUE)[["mean"]])),
#                       data.frame(time = time(precip),
#                                  var = "precip",
#                                  region = "CC",
#                                  value = as.numeric(
#                                    global(precip.CC,
#                                           mean,na.rm = TRUE)[["mean"]]))) %>%
#   pivot_wider(names_from = "region",
#               values_from = "value")
#
# saveRDS(raw.data,
#         "./outputs/precip_ERA5.RDS")
#

###############################

rm(list = ls())

library(tidyr)
library(dplyr)
library(cowplot)
library(ggplot2)
library(slider)

system2("scp",
        c("hpc:/kyukon/data/gent/vo/000/gvo00074/felicien/R/outputs/SST_ERA5.RDS",
        "./outputs/SST_ERA5.RDS"))
system2("scp",
        c("hpc:/kyukon/data/gent/vo/000/gvo00074/felicien/R/outputs/Omega_ERA5.RDS",
          "./outputs/Omega_ERA5.RDS"))
system2("scp",
        c("hpc:/kyukon/data/gent/vo/000/gvo00074/felicien/R/outputs/precip_ERA5.RDS",
          "./outputs/precip_ERA5.RDS"))

base_period <- c(1961, 1990)

SST <- readRDS("./outputs/SST_ERA5.RDS") %>%
  mutate(year = year(time),
         month = month(time)) %>%
  pivot_wider(values_from = "value",
              names_from = "var") %>%
  mutate(sst_av = slide_dbl(sst,mean, .before = 2, complete = TRUE, na.rm = TRUE)) %>%
  group_by(month) %>%
  mutate(sst.m = mean(sst[year %in% c(1961:1990)])) %>%
  mutate(anomaly.sst = sst - sst.m) %>%
  ungroup() %>%
  mutate(anomaly.sst_av = slide_dbl(anomaly.sst, mean, .before = 2, complete = TRUE, na.rm = TRUE))

ggplot(data = SST %>%
         filter(month == 11) %>%
         filter(year %in% base_period[1]:base_period[2]),
       aes(x = year + (month - 1/2)/12,
           y = anomaly.sst_av)) +
  geom_line() +
  geom_point() +
  stat_smooth(method = "lm",
              se = FALSE) +
  theme_bw()


B <- readRDS("./outputs/SST_ERA5.RDS") %>%
  mutate(year = year(time),
         month = month(time)) %>%
  pivot_wider(values_from = "value",
              names_from = "var") %>%
  mutate(sst_av = slide_dbl(sst,mean, .before = 2, complete = TRUE, na.rm = TRUE)) %>%
  group_by(month) %>%
  mutate(sst.m = mean(sst[year %in% c(1960:1990)])) %>%
  mutate(anomaly.sst = sst - sst.m) %>%
  ungroup() %>%
  mutate(anomaly.sst_av = slide_dbl(anomaly.sst, mean, .before = 2, complete = TRUE, na.rm = TRUE))

plotA <- ggplot(data = B %>%
         filter(month == 11) %>%
         filter(year <= 2024)) +
  geom_rect(data = data.frame(x = 1960, xend = 1990,
                       y = -Inf, yend = Inf),
            aes(xmin = x, xmax = xend, ymin = y, ymax = yend),
            fill = "grey", alpha = 0.5) +
  geom_line(aes(x = year + (month - 1/2)/12,
                y = anomaly.sst_av)) +
  geom_point(aes(x = year + (month - 1/2)/12,
                 y = anomaly.sst_av)) +
  stat_smooth(aes(x = year + (month - 1/2)/12,
                  y = anomaly.sst_av),
              se = FALSE) +
  labs(x = "") +
  theme_bw()



Omega <- readRDS("./outputs/Omega_ERA5.RDS") %>%
  mutate(year = year(time),
         month = month(time))


clim <- Omega %>%
    filter(year >= base_period[1], year <= base_period[2]) %>%
    group_by(month) %>%
    summarise(
      CATL_clim = mean(CATL, na.rm = TRUE),
      CC_clim   = mean(CC,   na.rm = TRUE),
      .groups   = "drop")

plot(clim$CATL_clim,type = 'l',ylim = c(-1,1)*0.1)
lines(clim$CC_clim, lty = 2)

Omega_anom <- Omega %>%
    left_join(clim, by = "month") %>%
    mutate(
      CATL_anom = CATL - CATL_clim,
      CC_anom   = CC   - CC_clim,

      omega_index = (CATL_anom - CC_anom)
    )  %>%
  mutate(omega_index_av = slide_dbl(omega_index, mean, .before = 2, complete = TRUE, na.rm = TRUE))

plotB <- ggplot() +
  geom_rect(data = data.frame(x = 1960, xend = 1990,
                              y = -Inf, yend = Inf),
            aes(xmin = x, xmax = xend, ymin = y, ymax = yend),
            fill = "grey", alpha = 0.5) +
  geom_line(data = Omega_anom %>%
              filter(month == 11),
            aes(x = year + (month - 1/2)/12,
                y = omega_index_av)) +
  geom_point(data = Omega_anom %>%
               filter(month == 11),
             aes(x = year + (month - 1/2)/12,
                 y = omega_index_av)) +

  stat_smooth(data = Omega_anom %>%
                filter(month == 11),
              aes(x = year + (month - 1/2)/12,
                  y = omega_index_av),
              se = FALSE) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = "") +
  stat_smooth(method = "lm", se = FALSE) +
  theme_bw()

df.all <- SST %>%
  left_join(Omega_anom %>%
              dplyr::select(-var),
            by = c("time","year","month"))

plotC <- ggplot(data = df.all %>%
         filter(month == 11),
       aes(y = anomaly.sst_av, x = omega_index_av)) +
  geom_point() +
  stat_smooth(method = "lm", se = FALSE) +
  theme_bw()

sqrt(summary(lm(data = df.all %>%
                  filter(month == 11),
   formula = anomaly.sst_av ~ omega_index))[["r.squared"]])

plot_grid(plotA,plotB,plotC,align = "hv",ncol = 1)


precip <- readRDS("./outputs/precip_ERA5.RDS") %>%
  mutate(year = year(time),
         month = month(time)) %>%
  mutate(precip_av = slide_dbl(CC, mean, .before = 2, complete = TRUE, na.rm = TRUE)) %>%
  group_by(month) %>%
  mutate(precip_av.m = mean(precip_av[year %in% c(1961:1990)])) %>%
  mutate(anomaly.precip = precip_av - precip_av.m) %>%
  ungroup() %>%
  mutate(anomaly.precip_av = slide_dbl(anomaly.precip, mean, .before = 2, complete = TRUE, na.rm = TRUE))



plotA2 <- ggplot(data = precip %>%
                  filter(month == 11) %>%
                  filter(year <= 2024)) +
  geom_rect(data = data.frame(x = 1960, xend = 1990,
                              y = -Inf, yend = Inf),
            aes(xmin = x, xmax = xend, ymin = y, ymax = yend),
            fill = "grey", alpha = 0.5) +
  geom_line(aes(x = year + (month - 1/2)/12,
                y = anomaly.precip_av)) +
  geom_point(aes(x = year + (month - 1/2)/12,
                 y = anomaly.precip_av)) +
  stat_smooth(aes(x = year + (month - 1/2)/12,
                  y = anomaly.precip_av),
              se = FALSE) +
  labs(x = "") +
  theme_bw()


all.combined <- df.all %>%
  left_join(precip %>%
              dplyr::select(-c(CC,CATL)),
            by = c("year","month","time"))


plotB2 <- ggplot(data = all.combined %>%
                   filter(month == 11),
                 aes(y = anomaly.precip_av, x = anomaly.sst_av)) +
  geom_point() +
  stat_smooth(method = "lm", se = FALSE) +
  theme_bw()


sqrt(summary(lm(data = all.combined %>%
                  filter(month == 11),
                formula = anomaly.precip_av ~ anomaly.sst_av))[["r.squared"]])

plotC2 <- ggplot(data = all.combined %>%
                  filter(month == 11),
                aes(y = anomaly.precip_av, x = omega_index_av)) +
  geom_point() +
  stat_smooth(method = "lm", se = FALSE) +
  theme_bw()


sqrt(summary(lm(data = all.combined %>%
                  filter(month == 11),
                formula = anomaly.precip_av ~ omega_index_av))[["r.squared"]])

plot_grid(plotA2,plotB2,plotC2,align = "hv",ncol = 1)


saveRDS(all.combined,
        "./outputs/ERA5.teleconnection.RDS")

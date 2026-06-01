rm(list = ls())


library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(slider)
library(cowplot)

# system2("scp",
#         c("hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/CMIP6.summary.RDS",
#           "./outputs/"))

CMIP6.summ <- readRDS("./outputs/CMIP6.summary.RDS") %>%
  ungroup() %>%
  filter(!(model %in% c("CIESM","FGOALS-g3","MCM-UA-1-0","SAM0-UNICON"))) %>%
  mutate(year = year(time),
         month = month(time)) %>%
  dplyr::select(-file) %>%
  arrange(time) %>%
  mutate(value = case_when(var == "pr" ~ value*86400*365/12,
                           TRUE ~value)) %>%
  filter((time < as.Date("2014/01/01") & scenario == "historical") |
           (time >= as.Date("2014/01/01") & scenario != "historical"))

CMIP6.long <- CMIP6.summ %>%
  distinct() %>%
  group_by(month,model,scenario,variant,var, region) %>%
  mutate(value.m = mean(value[year %in% 1981:2020],na.rm = TRUE)) %>%
  mutate(anomaly = value - value.m)

SST <- CMIP6.long %>%
  mutate(year = year(time),
         month = month(time)) %>%
  mutate(time = as.Date(paste0(year,"/",month,"/01"))) %>%
  filter(var == "tos") %>%
  ungroup() %>%
  dplyr::select(time, year, month, scenario, model, variant, region, value.m, anomaly) %>%
  group_by(model, variant, time,year, month,region, scenario) %>%
  summarise(value.m = mean(value.m,na.rm = TRUE),
            anomaly = mean(anomaly,na.rm = TRUE),
            .groups = "keep") %>%
  group_by(model, variant, region, scenario) %>%
  mutate(
    value.m.ma = slide_dbl(value.m, mean, .before = 2, complete = FALSE, na.rm = TRUE),
    anomaly.ma = slide_dbl(anomaly, mean, .before = 2, complete = FALSE, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(year >= 1960)

SST.m <- SST %>%
  group_by(time,variant, region, scenario) %>%
  summarise(
    value.m.ma = mean(value.m.ma,na.rm = TRUE),
    anomaly.ma = mean(anomaly.ma,na.rm = TRUE),
    .groups = "keep") %>%
  group_by(variant, region, scenario) %>%
  mutate(anomaly.ma.long = slide_dbl(anomaly.ma, mean, .before = 4, complete = FALSE, na.rm = TRUE)) %>%
  mutate(anomaly.ma.long = case_when(scenario == "historical" ~ anomaly.ma.long,
                                     TRUE ~ anomaly.ma.long + 0.4)) %>%
  mutate(anomaly.ma = case_when(scenario == "historical" ~ anomaly.ma,
                                     TRUE ~ anomaly.ma + 0.4))


ggplot(data = SST) +
  geom_hline(yintercept = 0,linetype = 2) +
  # geom_line(aes(x = time, y = anomaly.ma,
  #               color = scenario,
  #               group = interaction(model,scenario)),
  #           size = 0.1) +
  geom_line(data = SST.m,
            aes(x = time, y = anomaly.ma.long,
                color = scenario)) +
  scale_color_manual(values = c("black","#263b5d","#8b9bac","#b48a40","#6a2d31")) +
  theme_bw() +
  labs(x = "", y = "SST anomaly (°C)", color = "") +
  theme(text = element_text(size = 20),
        legend.position = c(0.12,0.75))

Omega <- CMIP6.long %>%
  filter(var == "wap") %>%
  ungroup() %>%
  dplyr::select(time, year, month, scenario, model, variant, region, anomaly) %>%
  pivot_wider(names_from = region, values_from = anomaly) %>%
  arrange(model, variant, time) %>%
  group_by(model, variant, scenario) %>%
  mutate(
    omega = CATL - CC,
    omega_av = slide_dbl(omega, mean, .before = 2, complete = FALSE, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(year >= 1960)

Omega.MEM <- Omega %>%
  group_by(year,month,scenario) %>%
  summarise(omega_av_MEM = mean(omega_av,na.rm = TRUE),
            .groups = "keep") %>%
  filter((scenario == "historical" & year <= 2014) | (scenario != "historical"))

ggplot(data = Omega.MEM %>%
         filter(month == 11),
       aes(x = year, y = omega_av_MEM,color = scenario)) +
  geom_line() +
  stat_smooth(se = FALSE) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_color_manual(values = c("black","#263b5d","#8b9bac","#b48a40","#6a2d31")) +
  theme_bw() +
  labs(x = "", y = "") +
  theme(legend.position = c(0.15,0.15))


precip <- CMIP6.long %>%
  filter(year >= 1960) %>%
  group_by(model,scenario,variant) %>%
  arrange(time) %>%
  filter(var == "pr",
         region == "CC") %>%
  mutate(pr_av = slide_dbl(value, mean, .before = 2, complete = TRUE, na.rm = TRUE),
         pr_anomaly_av = slide_dbl(anomaly, mean, .before = 2, complete = TRUE, na.rm = TRUE))

precip.MEM <- precip %>%
  group_by(year,month,scenario) %>%
  summarise(pr_av_MEM = mean(pr_av,na.rm = TRUE),
            pr_anomaly_av_MEM = mean(pr_anomaly_av,na.rm = TRUE),
            .groups = "keep") %>%
  filter((scenario == "historical" & year <= 2014) | (scenario != "historical"))

ggplot(data = precip.MEM %>%
         filter(month == 11),
       aes(x = year, y = pr_anomaly_av_MEM,color = scenario)) +
  geom_line() +
  stat_smooth(se = FALSE) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_color_manual(values = c("black","#263b5d","#8b9bac","#b48a40","#6a2d31")) +
  theme_bw() +
  labs(x = "", y = "") +
  theme(legend.position = c(0.12,0.9))


combined <- precip.MEM %>%
  left_join(Omega.MEM,
            by = c("year","month","scenario")) %>%
  left_join(SST.m %>%
             ungroup() %>%
             mutate(year = year(time),
                    month = month(time)) %>%
             rename(SST.anomaly = anomaly.ma) %>%
             dplyr::select(year,month,scenario,SST.anomaly) %>%
              distinct()) %>%
  group_by(scenario) %>%
  mutate(omega_dt = resid(lm(omega_av_MEM ~ year)),
         pr_dt    = resid(lm(pr_anomaly_av_MEM ~ year)),
         SST_dt    = resid(lm(SST.anomaly ~ year)))

ggplot(data = combined %>%
         filter(month == 11),
       aes(x = omega_av_MEM,
           y = pr_anomaly_av_MEM,
           color = scenario)) +
  geom_point() +
  stat_smooth(se = FALSE, method = "lm") +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_color_manual(values = c("black","#263b5d","#8b9bac","#b48a40","#6a2d31")) +
  theme_bw() +
  labs(x = "", y = "") +
  theme(legend.position = c(0.12,0.9))


combined %>%
  filter(month == 11, year >= 1960) %>%
  group_by(scenario) %>%
  summarise(s = coef(lm(pr_anomaly_av_MEM ~ omega_av_MEM))[2],
            r2 = summary(lm(pr_anomaly_av_MEM ~ omega_av_MEM))[["r.squared"]],
            r = sqrt(summary(lm(pr_anomaly_av_MEM ~ omega_av_MEM))[["r.squared"]])*sign(s),
            .groups = "keep")


ggplot(data = combined %>%
         filter(month == 11),
       aes(x = SST.anomaly,
           y = pr_anomaly_av_MEM,
           color = scenario)) +
  geom_point() +
  stat_smooth(se = FALSE, method = "lm") +

  stat_smooth(data = combined %>%
                filter(month == 11),
              inherit.aes = FALSE, aes(x = SST.anomaly,
                                       y = pr_anomaly_av_MEM),
              color = "grey", linetype = 2,
              se = FALSE, method = "lm") +

  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_color_manual(values = c("black","#263b5d","#8b9bac","#b48a40","#6a2d31")) +
  theme_bw() +
  labs(x = "SST anomaly (°C)",
       y = "CC Precip anomaly (mm)", color = "") +
  theme(legend.position = c(0.1,0.75),
        text = element_text(size = 20))

combined %>%
  filter(month == 11, year >= 1960) %>%
  ungroup() %>%
  # group_by(scenario) %>%
  summarise(s = coef(lm(pr_anomaly_av_MEM ~ SST.anomaly))[2],
            r2 = summary(lm(pr_anomaly_av_MEM ~ SST.anomaly))[["r.squared"]],
            r = sqrt(summary(lm(pr_anomaly_av_MEM ~ SST.anomaly))[["r.squared"]])*sign(s),
            .groups = "keep")

ggplot(data = combined %>%
         filter(month == 11),
       aes(x = SST_dt,
           y = pr_dt,
           color = scenario)) +
  geom_point() +
  stat_smooth(se = FALSE, method = "lm") +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_color_manual(values = c("black","#263b5d","#8b9bac","#b48a40","#6a2d31")) +
  theme_bw() +
  labs(x = "", y = "") +
  theme(legend.position = c(0.12,0.9))

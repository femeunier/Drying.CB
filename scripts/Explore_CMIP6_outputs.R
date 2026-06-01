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
#
# system2("scp",
#         c("hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/CMIP6.checks.RDS",
#           "./outputs/"))


CMIP6.checks <- terra::unwrap(readRDS("./outputs/CMIP6.checks.RDS"))
selected <- 150
names(CMIP6.checks)[selected]

CMIP6.summ <- readRDS("./outputs/CMIP6.summary.RDS") %>%
  ungroup() %>%
  filter(scenario == "historical") %>%
  filter(!(model %in% c("CIESM","FGOALS-g3","MCM-UA-1-0","SAM0-UNICON"))) %>%
  mutate(year = year(time),
         month = month(time)) %>%
  dplyr::select(-file) %>%
  arrange(time) %>%
  mutate(value = case_when(var == "pr" ~ value*86400*365/12,
                           TRUE ~value))

CMIP6.long <- CMIP6.summ %>%
  group_by(month,model,scenario,variant,var, region) %>%
  mutate(value.m = mean(value[year %in% 1981:2010],na.rm = TRUE)) %>%
  mutate(anomaly = value - value.m)

# CMIP6.long %>%
#   ungroup() %>%
#   filter(var == "pr",
#          scenario == "historical") %>%
#   filter(year == year[1]) %>%
#   group_by(region,model) %>%
#   summarise(Max.Pr = max(value.m,na.rm = TRUE)) %>%
#   arrange((Max.Pr))

CMIP6.long.m <- CMIP6.long %>%
  filter(year == year[1]) %>%
  group_by(var,region,scenario,month) %>%
  summarise(value.m.MEM = mean(value.m,na.rm = TRUE),
            value.m.sd = sd(value.m,na.rm = TRUE),
            .groups = "keep")


ggplot(data = CMIP6.long.m) +
  geom_ribbon(aes(x = month,
                  ymin = value.m.MEM - value.m.sd,
                  ymax = value.m.MEM + value.m.sd,
                  group = region),
              alpha = 0.5, fill = "grey") +
  geom_line(aes(x = month, y = value.m.MEM,
                linetype = region)) +
  facet_wrap(~var,ncol = 1,
             scales = "free_y") +
  scale_x_continuous(breaks = 1:12, labels= c("J","F","M","A","M","J",
                                              "J","A","S","O","N","D")) +
  theme_bw() +
  labs(linetype = "",x = "",y = "") +
  guides(color = "none") +
  theme(legend.position = c(0.1,0.94))


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
  ungroup()

Omega.MEM <- Omega %>%
  group_by(year,month,scenario) %>%
  summarise(omega_av_MEM = mean(omega_av,na.rm = TRUE),
            omega_av_sd = sd(omega_av,na.rm = TRUE),
            .groups = "keep")

ggplot(data = Omega.MEM %>%
         filter(month == 11,
                scenario == "historical",
                year >= 1980,
                year <= 2014),
       aes(x = year + (month - 1/2)/12,
           y = omega_av_MEM)) +
  # geom_line(aes(color = model),size = 0.1) +
  # geom_point(aes(color = model),size = 0.2) +

  # facet_wrap(~ scenario) +
  # stat_smooth(aes(color = model),
  #             se = FALSE,method = "lm",
  #             size = 0.5) +

  geom_line(color = "black") +

  geom_point(color = "black") +

  stat_smooth(data = Omega.MEM %>%
                filter(month == 11,
                       scenario == "historical",
                       year >= 1980,
                       year <= 2014),
              aes(x = year + (month - 1/2)/12,
                  y = omega_av_MEM), color = "black",
              se = FALSE,method = "lm",
              size = 0.5) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = 2) +
  guides(color = "none")


plotA <- ggplot() +
  geom_rect(data = data.frame(x = 1981, xend = 2010,
                              y = -Inf, yend = Inf),
            aes(xmin = x, xmax = xend, ymin = y, ymax = yend),
            fill = "grey", alpha = 0.5) +
  geom_line(data = Omega.MEM %>%
              filter(year <= 2014) %>%
              filter(month == 11,scenario == "historical") ,
            aes(x = year + (month - 1/2)/12,
                y = omega_av_MEM)) +
  geom_point(data = Omega.MEM %>%
               filter(year <= 2014) %>%
               filter(month == 11,scenario == "historical") ,
             aes(x = year + (month - 1/2)/12,
                 y = omega_av_MEM)) +
  geom_hline(yintercept = 0, linetype = 2) +
  stat_smooth(data = Omega.MEM %>%
                filter(year >= 1960,
                       year <= 2014) %>%
                filter(month == 11,scenario == "historical") ,
              aes(x = year + (month - 1/2)/12,
                  y = omega_av_MEM),
              method = "lm", color = "black", fill = "darkgrey",
              se = TRUE) +
  scale_x_continuous(limits = c(1960,2014)) +
  theme_bw() +
  labs(x = "", y = "") +
  guides(color = "none")

plotA

summary(lm(data = Omega.MEM %>%
          filter(year >= 1960,
                 year <= 2014) %>%
          filter(month == 11,scenario == "historical"),
          formula = omega_av_MEM~year))[["r.squared"]]


ggplot() +
  geom_rect(data = data.frame(x = 1991, xend = 2010,
                              y = -Inf, yend = Inf),
            aes(xmin = x, xmax = xend, ymin = y, ymax = yend),
            fill = "grey", alpha = 0.5) +
  geom_line(data = Omega %>%
              filter(year <= 2014) %>%
              filter(month == 11,scenario == "historical") ,
            aes(x = year + (month - 1/2)/12,
                y = omega_av,
                color = model)) +
  # geom_point(data = Omega %>%
  #              filter(year <= 2014) %>%
  #              filter(month == 11,scenario == "historical") ,
  #            aes(x = year + (month - 1/2)/12,
  #                y = omega, color = model)) +
  geom_hline(yintercept = 0, linetype = 2) +
  stat_smooth(data = Omega %>%
                filter(year <= 2014) %>%
                filter(month == 11,scenario == "historical") ,
              aes(x = year + (month - 1/2)/12,
                  y = omega_av,
                  color = model),
              method = "lm",
              se = FALSE) +
  scale_x_continuous(limits = c(1980,2014)) +
  theme_bw() +
  labs(x = "") +
  guides(color = "none")


Omega %>%
  filter(model == "NorCPM1",
         month == 11) %>%
  pull(omega_av) %>%
  plot()


SST <- SST <- CMIP6.long %>%
  filter(var == "tos", region == "CATL") %>%
  ungroup() %>%
  arrange(model, variant, scenario, time) %>%
  group_by(model, variant, scenario) %>%
  mutate(
    sst_av = slide_dbl(value, mean, .before = 2, complete = TRUE, na.rm = TRUE),
    sst_anomaly_av = slide_dbl(anomaly, mean, .before = 2, complete = TRUE, na.rm = TRUE)
  ) %>%
  ungroup()


ggplot() +

  geom_rect(data = data.frame(x = 1981, xend = 2010,
                              y = -Inf, yend = Inf),
            aes(xmin = x, xmax = xend, ymin = y, ymax = yend),
            fill = "grey", alpha = 0.5) +
  geom_line(data = SST %>%
              filter(month == 11,
                     scenario == "historical",
                     year <= 2014) ,
            aes(x = year + (month - 1/2)/12,
                y = sst_anomaly_av,
                color = model), size = 0.1) +
  geom_point(data = SST %>%
               filter(month == 11,
                      scenario == "historical",
                      year <= 2014) ,
             aes(x = year + (month - 1/2)/12,
                 y = sst_anomaly_av,
                 color = model), size = 0.2) +
  facet_wrap(~scenario) +
  geom_hline(yintercept = 0, linetype = 2) +
  stat_smooth(data = SST %>%
                filter(month == 11,
                       scenario == "historical",
                       year <= 2014) ,
              aes(x = year + (month - 1/2)/12,
                  y = sst_anomaly_av,
                  color = model),
              method = "loess", se = FALSE) +
  theme_bw() +
  scale_x_continuous(limits = c(1960,2014)) +
  guides(color = "none")


SST.MEM <- SST %>%
  group_by(year,month,scenario) %>%
  summarise(sst_anomaly_av_MEM = mean(sst_anomaly_av,na.rm = TRUE),
            sst_anomaly_av_sd = sd(sst_anomaly_av,na.rm = TRUE),
            .groups = "keep")


plotB <- ggplot() +
  geom_rect(data = data.frame(x = 1981, xend = 2010,
                              y = -Inf, yend = Inf),
            aes(xmin = x, xmax = xend, ymin = y, ymax = yend),
            fill = "grey", alpha = 0.5) +
  geom_line(data = SST.MEM %>%
              filter(year <= 2014) %>%
              filter(month == 11,scenario == "historical") ,
            aes(x = year + (month - 1/2)/12,
                y = sst_anomaly_av_MEM)) +
  geom_point(data = SST.MEM %>%
               filter(year <= 2014) %>%
               filter(month == 11,scenario == "historical") ,
             aes(x = year + (month - 1/2)/12,
                 y = sst_anomaly_av_MEM)) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = "", y = "") +
  stat_smooth(data = SST.MEM %>%
                filter(year <= 2014) %>%
                filter(month == 11,scenario == "historical") ,
              aes(x = year + (month - 1/2)/12,
                  y = sst_anomaly_av_MEM), color = "black", fill = "darkgrey",
              se = TRUE) +
  scale_x_continuous(limits = c(1960,2014)) +
  theme_bw() +
  guides(color = "none")

plotB

combined <- SST.MEM %>%
  filter(year <= 2014) %>%
  filter(month == 11,scenario == "historical")  %>%
  left_join(Omega.MEM %>%
              filter(year <= 2014) %>%
              filter(month == 11,scenario == "historical"),
            by = c("year","month","scenario"))

combined.models <- SST %>%
  ungroup() %>%
  rename(tos = value,
         tos.anomaly = anomaly) %>%
  left_join(Omega %>%
              ungroup() %>%
              dplyr::select(-c(CATL,CC)),
            by = c("time","year","month","scenario",
                   "model","variant")) %>%
  filter(!is.na(omega_av) & !is.na(sst_anomaly_av))

ggplot(data = combined.models %>%
         filter(month == 11),
       aes(x = sst_anomaly_av,
           y = omega_av,
           color = model)) +
  geom_point(size = 0.1) +
  stat_smooth(method = "lm",se = FALSE) +
  stat_smooth(data = combined.models %>%
                filter(month == 11),
              aes(x = sst_anomaly_av,
                  y = omega_av), color = "black",
              method = "lm",se = FALSE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_bw()

df.slopes <- combined.models %>%
  filter(month == 11) %>%
  group_by(model) %>%
  summarise(s = coef(lm(omega_av ~sst_anomaly_av))[2],
            r2 = summary(lm(omega_av ~sst_anomaly_av))[["r.squared"]],
            r = as.numeric(sign(coef(lm(omega_av ~sst_anomaly_av))[2]) * sqrt(r2)))

hist(df.slopes$r)

plotC <- ggplot(data = combined %>%
                  filter(year >= 1960),
                aes(x = sst_anomaly_av_MEM, y = omega_av_MEM)) +
  geom_point() +
  stat_smooth(method = "lm", se = TRUE, color = "black", fill = "darkgrey") +
  theme_bw()

plotC

sqrt(summary(lm(data = combined %>%
                  filter(year >= 1960),
                formula = sst_anomaly_av_MEM ~ omega_av_MEM))[["r.squared"]])


plot_grid(plotA,plotB,plotC,align = "hv",ncol = 1)

precip <- CMIP6.long %>%
  group_by(model,scenario,variant) %>%
  arrange(time) %>%
  filter(var == "pr",
         region == "CC") %>%
  mutate(pr_av = slide_dbl(value, mean, .before = 2, complete = TRUE, na.rm = TRUE),
         pr_anomaly_av = slide_dbl(anomaly, mean, .before = 2, complete = TRUE, na.rm = TRUE))

precip.MEM <- precip %>%
  group_by(year,month,scenario) %>%
  summarise(pr_anomaly_av_MEM = mean(pr_anomaly_av,na.rm = TRUE),
            pr_anomaly_av_sd = sd(pr_anomaly_av,na.rm = TRUE),
            .groups = "keep")

df2plot <- precip %>%
  filter(scenario == "historical",
         month == 11,
         year >= 1960,
         year <= 2014) %>%
  ungroup() %>%
  mutate(time = year + (month - 1/2)/12)

plotCbis <- ggplot() +
  geom_rect(data = data.frame(x = 1981, xend = 2010,
                              y = -Inf, yend = Inf),
            aes(xmin = x, xmax = xend, ymin = y, ymax = yend),
            fill = "grey", alpha = 0.5) +
  # geom_line() +
  # geom_point() +
  # stat_smooth(data = df2plot,
  #             aes(x = time,
  #                 y = pr_anomaly_av,
  #                 color = model),size = 0.5,
  #             method = "lm", se = FALSE) +
  geom_line(data = precip.MEM %>%
               filter(scenario == "historical",
                      month == 11,
                      year >= 1960,
                      year <= 2014),
             aes(x = year + (month - 1/2)/12, y = pr_anomaly_av_MEM),
             color = "black") +
  geom_point(data = precip.MEM %>%
              filter(scenario == "historical",
                     month == 11,
                     year >= 1960,
                     year <= 2014),
            aes(x = year + (month - 1/2)/12, y = pr_anomaly_av_MEM),
            color = "black") +
  stat_smooth(data = precip.MEM %>%
                filter(scenario == "historical",
                       month == 11,
                       year >= 1960,
                       year <= 2014),
              aes(x = year + (month - 1/2)/12,
                  y = pr_anomaly_av_MEM),color = "black",fill = "darkgrey",
              method = "lm", se = TRUE) +
  geom_hline(yintercept = 0,linetype = 2) +
  theme_bw() +
  labs(x = "", y = "") +
  guides(color = "none")

plotCbis


plot_grid(plotA,plotB,plotCbis,
          align = "hv",ncol = 1)

df2plot %>%
  group_by(model) %>%
  summarise(s = coef(lm(pr_anomaly_av ~ time))[2]) %>%
  pull(s) %>%
  hist()


all.combined <- combined %>%
  left_join(precip.MEM,
            by = c("year","month","scenario")) %>%
  filter(year >= 1960, month == 11)

plot(all.combined$omega_av_MEM)
plot(all.combined$pr_anomaly_av_MEM)

all.combined.models <- combined.models %>%
  left_join(precip %>%
              ungroup() %>%
              dplyr::select(-c(var,region)) %>%
              rename(pr = value,
                     anomaly_pr = anomaly),
            by = c("time","year","month","scenario",
                   "model","variant")) %>%
  filter(!is.na(pr))


ggplot(data = all.combined.models %>%
         filter(month == 11,
                year >= 1960),
       aes(y = pr_anomaly_av,
           x = omega_av,
           color = model)) +
  geom_point(size = 0.1) +
  stat_smooth(method = "lm",se = FALSE) +
  stat_smooth(data = all.combined.models %>%
                filter(month == 11),
              aes(y = pr_anomaly_av,
                  x = omega_av), color = "black",
              method = "lm",se = FALSE) +
  theme_bw()

df.slopes2 <- all.combined.models %>%
  filter(month == 11,
         year >= 1960) %>%
  group_by(model) %>%
  summarise(s = coef(lm(sst_anomaly_av ~ omega_av))[2],
            r2 = summary(lm(sst_anomaly_av ~ omega_av))[["r.squared"]],
            r = as.numeric(sign(coef(lm(sst_anomaly_av ~ omega_av))[2]) * sqrt(r2)))

hist(df.slopes2$r)


ggplot(data = all.combined %>%
         filter(month == 11,
                year >= 1960),
       aes(x = omega_av_MEM, y = pr_anomaly_av_MEM)) +
  geom_point() +
  geom_hline(yintercept = 0,linetype = 2) +
  geom_vline(xintercept = 0,linetype = 2) +
  stat_smooth(method = "lm", se = TRUE,
              color = "black", fill = "grey") +
  labs(x="",y = "") +
  theme_bw()

(summary(lm(data = all.combined %>%
                  filter(month == 11,
                         year >= 1960),
                formula = pr_anomaly_av_MEM ~ omega_av_MEM))[["r.squared"]])

all.combined <- all.combined %>%
  ungroup() %>%
  mutate(
    sst_dt = resid(lm(sst_anomaly_av_MEM ~ year)),
    omega_dt = resid(lm(omega_av_MEM ~ year)),
    pr_dt    = resid(lm(pr_anomaly_av_MEM ~ year))
  )

ggplot(all.combined %>%
         filter(year >= 1960,
                month == 11),
       aes(x = omega_dt , y = pr_dt)) +
  geom_point() +
  stat_smooth(method = "lm",se = TRUE,
              color = "black", fill = "grey") +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x="",y = "") +
  theme_bw()


(summary(lm(data = all.combined %>%
                  filter(month == 11,
                         year >= 1960),
                formula = pr_dt ~ omega_dt))[["r.squared"]])

slopes <- precip %>%
  ungroup() %>%
  filter(year >= 1980,year <= 2014) %>%
  mutate(year = year + (month - 1/2)/12) %>%
  group_by(model) %>%
  summarise(s = coef(lm(pr_av ~time))[2])

hist(slopes$s*10*12)

system2("scp",
        c("hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/CMIP6.precip.RDS",
          "./outputs/"))


A <- readRDS("./outputs/CMIP6.precip.RDS") %>%
  filter(!(model %in% c("CIESM","FGOALS-g3","MCM-UA-1-0"))) %>%
  mutate(value = value*86400*365/12,
         month = month(time),
         year = year(time)) %>%
  group_by(month,model,scenario,variant,var, region) %>%
  mutate(value.m = mean(value[year %in% 1981:2010],na.rm = TRUE)) %>%
  mutate(anomaly = value - value.m)

ggplot(A %>%
         group_by(model,variant) %>%
         filter(year == year[1],
                scenario == "historical")) +
  geom_line(aes(x = month, y = value.m,
                color = model)) +
  facet_wrap(~region, nrow = 1) +
  geom_hline(yintercept = 100, linetype = 2) +
  theme_bw() +
  guides(color = "none")

precip <- A %>%
  group_by(model,scenario,variant,region) %>%
  arrange(time) %>%
  mutate(pr_av = slide_dbl(value, mean, .before = 11, complete = TRUE, na.rm = TRUE),
         pr_anomaly_av = slide_dbl(anomaly, mean, .before = 2, complete = TRUE, na.rm = TRUE))

precip.MEM <- A %>%
  group_by(year,month,scenario,variant,region) %>%
  summarise(pr.m = mean(value,na.rm = TRUE),
            .groups = "keep") %>%
  group_by(scenario,variant,region) %>%
  arrange(year,month) %>%
  mutate(pr.m_av = slide_dbl(pr.m, mean, .before = 11, complete = TRUE, na.rm = TRUE))


ggplot() +
  geom_line(data = precip %>%
              filter(scenario == "historical",
                     year <= 2014),
            aes(x = year + (month - 1/2)/12,
                y = pr_av,
                color = model)) +
  stat_smooth(data = precip %>%
                filter(scenario == "historical",
                       year <= 2014),
              aes(x = year + (month - 1/2)/12,
                  y = pr_av,
                  color = model),
              method = "lm", se = FALSE) +
  geom_line(data = precip.MEM %>%
              filter(scenario == "historical",
                     month == 11,
                     year <= 2014),
            aes(x = year + (month - 1/2)/12,
                y = pr.m_av),
            color = "black") +
  stat_smooth(data = precip.MEM %>%
                filter(scenario == "historical",
                       month == 11,
                       year <= 2014),
              aes(x = year + (month - 1/2)/12,
                  y = pr.m_av), color = "black",
              method = "loess", se = FALSE) +
  scale_x_continuous(limits = c(1980,2014)) +
  facet_wrap(~region, nrow = 1,
             scales = "free_y") +
  theme_bw() +
  guides(color = "none")



ggplot() +
  # geom_line(data = precip %>%
  #             filter(scenario != "historical",
  #                    year >= 2014),
  #           aes(x = year + (month - 1/2)/12,
  #               y = pr_av,
  #               color = model)) +
  # stat_smooth(data = precip %>%
  #               filter(scenario != "historical",
  #                      year >= 2014),
  #             aes(x = year + (month - 1/2)/12,
  #                 y = pr_av,
  #                 color = model),
  #             method = "lm", se = FALSE) +
  geom_line(data = precip.MEM %>%
              filter(month == 11,
                     scenario != "historical",
                     year >= 2014),
            aes(x = year + (month - 1/2)/12,
                y = pr.m_av*12),
            color = "black") +
  stat_smooth(data = precip.MEM %>%
                filter(scenario != "historical",
                       month == 11,
                       year >= 2014),
              aes(x = year + (month - 1/2)/12,
                  y = pr.m_av*12), color = "black",
              method = "loess", se = FALSE) +
  scale_x_continuous(limits = c(2014,2100)) +
  facet_grid(region ~ scenario,
             scales = "free_y") +
  theme_bw() +
  guides(color = "none")



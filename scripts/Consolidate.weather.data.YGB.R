rm(list = ls())

library(ggplot2)
library(tidyr)
library(dplyr)

A <- readRDS("/home/femeunier/Documents/projects/Precip.Africa/data/YGB/climate.YGB.RDS")
unique(A$year)

A.long <- A %>%
  mutate(tmp = (tmin + tmax)/2) %>%
  dplyr::select(-c(tmin,tmax)) %>%
  pivot_longer(cols = c(tmp,Pmm),
               names_to = "variable",
               values_to = "value")
A.sum <- A.long %>%
  group_by(month,variable) %>%
  summarise(value.m = mean(value,na.rm = TRUE),
            .groups = "keep")

ggplot(data = A.sum) +
  geom_line(aes(x = month, y = value.m)) +
  facet_wrap(~ variable, scales = "free") +
  theme_bw()

raw <- read.csv("/home/femeunier/Downloads/ICOSETC_CD-Ygb_METEOSENS_L2/ICOSETC_CD-Ygb_METEOSENS_L2.csv")
B <- raw %>%
  mutate(year = as.numeric(substr(TIMESTAMP_START,1,4)),
         month = as.numeric(substr(TIMESTAMP_START,5,6)),
         day = as.numeric(substr(TIMESTAMP_START,7,8)),
         hour = as.numeric(substr(TIMESTAMP_START,9,10)),
         min = as.numeric(substr(TIMESTAMP_START,11,12))) %>%
  dplyr::select("year","month","day","hour","min",
                c(starts_with("P_1_1_1"),
                  starts_with("TA_1_2_1")))
B.long <- B %>%
  pivot_longer(cols = -c("year","month","day",
                         "hour","min"),
               names_to = "variable",
               values_to = "value") %>%
  filter(value != -9999)

B.long.pp <- B.long %>%
  group_by(year,month,day,variable) %>%
  summarise(value.m = case_when(variable[1] == "P_1_1_1" ~ sum(value,na.rm = TRUE),
                                TRUE ~ mean(value,na.rm = TRUE)),
            N = length(which(!is.na(value))),
            .groups = "keep") %>%
  group_by(year,month,variable) %>%
  summarise(value.m = case_when(variable[1] == "P_1_1_1" ~ sum(value.m,
                                                               na.rm = TRUE),
                                TRUE ~ mean(value.m,
                                            na.rm = TRUE)),
            N = sum(N,na.rm = TRUE),
            .groups = "keep") %>%
  mutate(Nth = lubridate::days_in_month(as.Date(
    paste0(year,"/",month,"/01")))*48) %>%
  filter((variable == "TA_1_2_1" & N > 0.7*Nth) |
           (variable == "P_1_1_1" & N > 0.95*Nth))


B.sum <- B.long.pp %>%
  group_by(month,variable) %>%
  summarise(value.m = mean(value.m,na.rm = TRUE),
            .groups = "keep")

all <- bind_rows(A.long %>%
                   mutate(origin = "Kasongo") %>%
                   rename(value.m = value),
                 B.long.pp %>%
                   mutate(variable = case_when(variable == "P_1_1_1" ~ "Pmm",
                                               variable == "TA_1_2_1" ~ "tmp")) %>%
                   mutate(origin = "FT"))

ggplot(data = all) +
  geom_line(aes(x = year + (month -1/2)/12, y = value.m,
                color = origin)) +
  facet_wrap(~ variable, scales = "free") +
  theme_bw()


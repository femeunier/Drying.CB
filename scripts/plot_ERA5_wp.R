rm(list = ls())

library(lubridate)
library(dplyr)
library(tidyr)


system2("scp",
        c("hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/Monthly.maps.Congo.wp.RDS",
          "./outputs/"))

system2("scp",
        c("hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/ts.Congo.wp.RDS",
          "./outputs/"))

CATL <- ext(-23,3,-2.4,5.5)
CC <- ext(22,29,-10,4)

Monthly.maps.Congo.wp <- readRDS("./outputs/Monthly.maps.Congo.wp.RDS")
times <- as.Date(time(Monthly.maps.Congo.wp))


df.all <- data.frame()
for(itime in seq(1,length(times))){

  cR <- Monthly.maps.Congo.wp[[itime]]

  cR.CATL <- crop(cR,CATL)
  cR.CC <- crop(cR,CC)



  raw <- bind_rows(
    data.frame(time = time(cR),
               region = "CATL",
               value = as.numeric(global(cR.CATL, mean, na.rm = TRUE)[["mean"]])),
    data.frame(time = time(cR),
               region = "CC",
               value = as.numeric(global(cR.CC, mean, na.rm = TRUE)[["mean"]])))

  df.all <- bind_rows(df.all,
                      raw)


}

df.seasonal <- df.all %>%
  mutate(month = month(time)) %>%
  group_by(region, month) %>%
  summarise(value.m = mean(value,na.rm = TRUE),
            .groups = "keep")

df.seasonal.diff <- df.seasonal %>%
  pivot_wider(names_from = region,
              values_from = value.m)

ggplot(data = df.seasonal) +
  geom_line(aes(x = month, y = value.m, color = region)) +
  geom_hline(linetype = 2, yintercept = 0) +
  theme_bw()



ts.Congo.wp <- readRDS("./outputs/ts.Congo.wp.RDS")
ts.Congo.wp.m <- ts.Congo.wp %>%
  group_by(month) %>%
  summarise(value.m = mean(value))


A <- readRDS("./outputs/CMIP6.summary.RDS")

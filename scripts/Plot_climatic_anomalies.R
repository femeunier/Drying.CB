rm(list = ls())

library(terra)
library(lubridate)
library(slider)

system2("rsync",
        c("-avz","hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/Drying.CB/*_pre_Zanomalies_Afrorainforests.tif",
          "./outputs/"))
#
# system2("rsync",
#         c("-avz","hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/Drying.CB/MEM_pre_anomalies_Afrorainforests.tif",
#           "./outputs/"))
#
# system2("rsync",
#         c("-avz","hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/all.climate/MEM_pre_all.years.tif",
#           "./outputs/"))


Mask <- read_sf("data/shapefiles/Rainforests.shp")

E <- ext(-180,180,-25,25)

Zscores <- crop(mask(rast("./outputs/MEM_pre_Zanomalies_Afrorainforests.tif"),Mask),
                E)
Anomalies <- crop(mask(rast("./outputs/MEM_pre_Zanomalies_Afrorainforests.tif"),Mask),
                  E)
Abs <- crop(mask(rast("./outputs/MEM_pre_all.years.tif"),Mask),
            E)

pos <- which(year(time(Zscores)) == 2024 &
               month(time(Zscores)) == 1)

df <- data.frame(time = time(Zscores),
                 pre = global(Abs,mean,na.rm = TRUE)[["mean"]],
                 anomaly = global(Anomalies,mean,na.rm = TRUE)[["mean"]],
                 zscore = global(Zscores,mean,na.rm = TRUE)[["mean"]]) %>%
  mutate(year = year(time),
         month = month(time)) %>%
  arrange(time) %>%
  mutate(
    pre_cumsum_12 = slide_dbl(pre, ~sum(.x, na.rm = TRUE), .before = 11, .complete = TRUE),
    anomaly_cumsum_12 = slide_dbl(anomaly, ~sum(.x, na.rm = TRUE), .before = 11, .complete = TRUE),
    zscore_mean_12 = slide_dbl(zscore, ~mean(.x, na.rm = TRUE), .before = 11, .complete = TRUE))


ggplot(data = df %>%
         filter(year >= 1960)) +
  geom_rect(inherit.aes = FALSE,
            aes(xmin = as.Date("1985/01/01"),xmax = as.Date("2015/01/01"),
                ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.4) +
  geom_line(aes(x = time, y = pre),
            size = 0.1) +
  geom_line(aes(x = time, y = pre_cumsum_12/12)) +
  theme_bw()

ggplot(data = df %>%
         filter(year >= 1960 & year <= 2024)) +
  geom_rect(inherit.aes = FALSE,
            aes(xmin = as.Date("1985/01/01"),xmax = as.Date("2015/01/01"),
                ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.4) +
  geom_line(aes(x = time, y = anomaly),
            size = 0.1) +
  geom_line(aes(x = time, y = anomaly_cumsum_12/12)) +
  geom_point(data = df %>%
               filter(year == 2023, month %in% c(11,12)),
             aes(x = time, y = anomaly), color = "darkblue") +
  geom_line(data = df %>%
               filter(year >= 2024 & year <= 2024),
             aes(x = time, y = anomaly), color = "darkred") +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw()

ggplot(data = df %>%
         filter(year >= 1980 & year <= 2024)) +
  geom_rect(inherit.aes = FALSE,
            aes(xmin = as.Date("1985/01/01"),xmax = as.Date("2015/01/01"),
                ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.4) +
  geom_line(aes(x = time, y = zscore),
            size = 0.1) +
  geom_line(aes(x = time, y = zscore_mean_12)) +
  geom_line(data = df %>%
              filter(year >= 2024 & year <= 2024),
            aes(x = time, y = zscore_mean_12),
            color = "darkred") +
  geom_point(data = df %>%
               filter(year == 2023, month %in% c(11,12)),
             aes(x = time, y = zscore), color = "darkblue") +
  geom_line(data = df %>%
              filter(year >= 2024 & year <= 2024),
            aes(x = time, y = zscore), color = "darkred", size = 0.1) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw()


ggplot(data = df %>%
         filter(year >= 2023 & year <= 2024)) +
  geom_line(aes(x = time, y = zscore),
            size = 0.1) +
  geom_line(aes(x = time, y = zscore_mean_12)) +
  geom_point(data = df %>%
               filter(year == 2023, month %in% c(11,12)),
             aes(x = time, y = zscore), color = "darkblue") +
  geom_line(data = df %>%
              filter(year >= 2024 & year <= 2024),
            aes(x = time, y = zscore), color = "darkred") +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw()


ggplot(data = df %>%
         filter(year >= 1960)) +
  geom_line(aes(x = month, y = pre,
                group = as.factor(year)),
            size = 0.1) +
  geom_line(data = df %>%
              filter(year %in% 1985:2014) %>%
              group_by(month) %>%
              summarise(pre = mean(pre,na.rm = TRUE),
                        .groups = "keep"),
            aes(x = month, y = pre), size = 1.4) +
  geom_line(data = df %>%
              filter(year %in% c(2020,2023:2024)),
            aes(x = month, y = pre,
                color = as.factor(year)),size = 1) +
  theme_bw()

ref <- df %>%
  filter(year %in% 1985:2014) %>%
  group_by(month) %>%
  summarise(pre = mean(pre,na.rm = TRUE),
            .groups = "keep") %>%
  pull(pre)

df.2024 <- df %>%
  filter(year %in% 2024) %>%
  pull(pre)

mean((df.2024 - ref)/ref)

ggplot(data = df %>%
         filter(year >= 1960)) +
  geom_line(aes(x = month, y = anomaly,
                group = as.factor(year)),
            size = 0.1) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(data = df %>%
              filter(year %in% 2023:2024),
            aes(x = month, y = (anomaly),
                color = as.factor(year)),size = 1) +
  theme_bw()

df.cumsum <- df %>%
  group_by(year) %>%
  arrange(month) %>%
  mutate(anomaly.cumsum = cumsum(anomaly))

ggplot(data = df.cumsum %>%
         filter(year >= 1960)) +
  geom_line(aes(x = month, y = anomaly.cumsum,
                group = as.factor(year)),
            size = 0.1) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(data = df.cumsum %>%
              ungroup() %>%
              filter(year %in% c(df.cumsum %>%
                                   ungroup() %>%
                                   filter(year >= 1960) %>%
                                   filter(month == 12) %>%
                                   arrange(anomaly.cumsum) %>%
                                   slice_head(n = 10) %>%
                                   pull(year))),
            aes(x = month, y = anomaly.cumsum,
                group = year,
                color = as.factor(year)),size = 1) +
  theme_bw()


ggplot(data = df.cumsum %>%
         filter(year >= 1960)) +
  geom_line(aes(x = month, y = anomaly.cumsum,
                group = as.factor(year)),
            size = 0.1) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(data = df.cumsum %>%
              ungroup() %>%
              filter(year %in% c(df.cumsum %>%
                                   ungroup() %>%
                                   filter(year %in% c(2023:2024)) %>%
                                   filter(month == 12) %>%
                                   arrange(anomaly.cumsum) %>%
                                   slice_head(n = 10) %>%
                                   pull(year))),
            aes(x = month, y = anomaly.cumsum,
                group = year,
                color = as.factor(year)),size = 1) +
  theme_bw()




pos <- which(year(time(Zscores)) == 2024)

df.map <- as.data.frame(mean(Zscores[[pos]]),
                        xy = TRUE) %>%
  rename(lon = x,
         lat = y)


world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
CB <- read_sf("./data/shapefiles/CongoBasin.shp")

ggplot() +
  geom_tile(data = df.map,
             aes(x = lon, y = lat,
                 fill = mean)) +
  geom_sf(data = world, fill = NA, color = "black", linewidth = 0.3) +
  theme_bw() +
  coord_sf(xlim = c(-5, 35), ylim = c(-10, 10), expand = FALSE) +
  scale_fill_gradient2(limits = c(-1,1)*1.5,oob = scales::squish,
                       low = "darkred",high = "darkblue") +
  labs(
    fill = "Aveg. Zscore (2024)",
    x = "",
    y = ""
  ) +
  theme_void()



thr_dry <- -1
thr_wet <-  0.5
window_size <- 6   # months

whiplash_in_window <- df %>%
  arrange(time) %>%
  group_modify(~{

    x <- .x
    n <- nrow(x)
    out <- vector("list", n)

    for (k in 1:n) {

      # rolling window centered on k
      i1 <- max(1, k - floor(window_size))
      i2 <- min(n, k)

      w <- x[i1:i2, ]

      # all pairs within the window
      pairs <- expand.grid(i = 1:nrow(w), j = 1:nrow(w)) %>%
        filter(j > i) %>%
        mutate(
          start_time = w$time[i],
          end_time   = w$time[j],
          start_val  = w$zscore[i],
          end_val    = w$zscore[j],
          lag_months = j - i,
          type = case_when(
            start_val <= thr_dry & end_val >= thr_wet ~ "dry_to_wet",
            start_val >= thr_wet & end_val <= thr_dry ~ "wet_to_dry",
            TRUE ~ NA_character_
          ),
          magnitude = case_when(
            !is.na(type) ~ abs(end_val - start_val),
            TRUE ~ NA_real_
          )
        ) %>%
        filter(!is.na(type))

      if (nrow(pairs) == 0) {
        out[[k]] <- data.frame(
          time = x$time[k],
          window_start = x$time[i1],
          window_end = x$time[i2],
          whiplash_max = NA_real_,
          event_start = as.Date(NA),
          event_end = as.Date(NA),
          start_val = NA_real_,
          end_val = NA_real_,
          lag_months = NA_real_,
          type = NA_character_
        )
      } else {
        best <- pairs %>%
          slice_max(magnitude, n = 1, with_ties = FALSE)

        out[[k]] <- data.frame(
          time = x$time[k],
          window_start = x$time[i1],
          window_end = x$time[i2],
          whiplash_max = best$magnitude,
          event_start = best$start_time,
          event_end = best$end_time,
          start_val = best$start_val,
          end_val = best$end_val,
          lag_months = best$lag_months,
          type = best$type
        )
      }
    }

    bind_rows(out)
  }) %>%
  ungroup()

ggplot(data = whiplash_in_window %>%
         filter(year(time) >= 0),
       aes(x = time, y = whiplash_max)) +
  geom_line() +
  # stat_smooth(method = "lm",se = FALSE) +
  theme_bw()



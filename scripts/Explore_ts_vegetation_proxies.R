rm(list = ls())

library(slider)
library(ggplot2)
library(sf)
library(dplyr)
library(tidyr)
library(lubridate)
library(raster)
library(terra)

Ext <- ext(-180,180,-25,25)

IF <- crop(rast("./data/shapefiles/Intact_mask_mod.tif"),
           Ext)

Radar.data <- rast("./outputs/Radar_all.years.tif")
# IF <- rasterize(vect(read_sf("~/Documents/projects/Drying.CB/data/shapefiles/Rainforests.shp")),
#                 Radar.data, field = "LC")
IF <- rasterize(vect(read_sf("~/Documents/projects/Drying.CB/data/shapefiles/Rainforests.shp")), Radar.data, field = "LC")

plot(IF)

files <- list.files("./data/GPP",pattern = "gppnormanomaly*",full.names = TRUE)


df.models <- data.frame()
all.r <- list()
for (ifile in seq(1,length(files))){
  cfile <- files[ifile]

  cproduct <- strsplit(tools::file_path_sans_ext(basename(cfile)),"\\.")[[1]][2]
  R <- rast(cfile) # This is Z-anomaly already
  times <-as.Date(paste0(gsub("_", "-", names(R)), "-01"))

  crast.r <- resample(R, IF, method = "bilinear")
  ccrast.r <- mask(crast.r,IF,maskvalues=c(0,NA))

  terra::time(ccrast.r) <- as.Date(times)
  names(ccrast.r) <- paste0(cproduct,"_",1:nlyr(ccrast.r))

  all.r[[ifile]] <- ccrast.r

  terra::time(ccrast.r) <- as.Date(times)

  df.models <- bind_rows(df.models,
                         data.frame(model = cproduct,
                                    time = times,
                                    var = "gpp",
                                    anom.norm = global(ccrast.r,mean,na.rm = TRUE)[["mean"]]))



}

all.norm.anom <- rast(all.r)

pos <- (year(time(all.norm.anom)) == 2024 &
          month(time(all.norm.anom)) >= 6)

r.sel <- all.norm.anom[[pos]]

model.names <- as.factor(sapply(strsplit(names(r.sel), "_"), "[[", 1))
mdls <- as.numeric(model.names)

anomalies.per.model <- tapp(r.sel, index = mdls, fun = mean, na.rm = TRUE)
names(anomalies.per.model) <- levels(model.names)

df.anomalies.per.model <- as.data.frame(anomalies.per.model, xy = TRUE) %>%
  rename(lon = x, lat = y) %>%
  pivot_longer(cols = -c(lon, lat),
               names_to = "model")

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  geom_tile(data = df.anomalies.per.model,
            aes(x = lon, y = lat, fill = value)) +
  geom_sf(data = world, fill = NA, color = "black", linewidth = 0.3) +
  scale_fill_gradient2(low = "darkred",mid = "white",high = "darkblue",
                       midpoint = 0, oob = scales::squish) +
  theme_bw() +
  facet_wrap(~model) +
  coord_sf(xlim = c(5, 35), ylim = c(-10, 5), expand = FALSE) +
  labs(
    fill = "Z-anomalies"
  )


Nsmooth = 6
df.models <- df.models %>%
  mutate(anom.norm_ma = slide_dbl(anom.norm, mean, .before = (Nsmooth-1),
                               .complete = TRUE, na.rm = TRUE))

ggplot(data = df.models %>%
         filter(year(time) < 2025),
       aes(x = year(time) + (month(time) - 1/2)/12, y = anom.norm,
           color = model)) +
  geom_line(size = 0.5) +
  geom_line(data = df.models %>%
              filter(year(time) < 2025) %>%
              group_by(time) %>%
              summarise(anom.norm = mean(anom.norm,na.rm = TRUE),
                        model = "MEM",
                        .groups = "keep"),
            color = "black") +
  scale_x_continuous(limits = c(2022,2026)) +
  geom_hline(yintercept = 0,linetype = 2) +
  facet_wrap(~var) +
  labs(x = "") +
  theme_bw()

df.models %>%
  group_by(time) %>%
  filter(year(time) < 2025) %>%
  summarise(anom.norm = mean(anom.norm,na.rm = TRUE),
            model = "MEM",
            .groups = "keep") %>%
  ungroup() %>%
  filter(anom.norm == min(anom.norm))

ggplot(data = df.models,
       aes(x = year(time) + (month(time) - 1/2)/12, y = anom.norm_ma,
           color = model)) +
  geom_line(size = 0.5) +
  geom_line(data = df.models %>%
              group_by(time) %>%
              summarise(anom.norm_ma = mean(anom.norm_ma,na.rm = TRUE),
                        model = "MEM",
                        .groups = "keep"),
            color = "black") +
  scale_x_continuous(limits = c(2022,2025)) +
  geom_hline(yintercept = 0,linetype = 2) +
  facet_wrap(~var) +
  labs(x = "") +
  theme_bw()


thr_dry <- 1
thr_wet <-  0.5
window_size <- 12   # months

whiplash_in_window <- df.models %>%
  arrange(model, var, time) %>%
  group_by(model, var) %>%
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
          start_val  = w$anom.norm[i],
          end_val    = w$anom.norm[j],
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

ggplot(data = whiplash_in_window,
       aes(x = time, y = whiplash_max, color = model)) +
  geom_line() +
  # stat_smooth(method = "lm",se = FALSE) +
  theme_bw()

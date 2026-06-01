rm(list = ls())

library(terra)
library(lubridate)
library(tidyr)

E <- ext(-180,180,-25,25)

Zscores <- crop(rast("./outputs/chirpsv3_pre_Zanomalies_CA.tif"),
                E)

ERA5_pre_Zanomalies_CA.tif

pos <- which(year(time(Zscores)) >= 1980)
r <- Zscores[[pos]]
names(r) <- format(time(r), "%Y_%m")

df <- as.data.frame(r, xy = TRUE)

# convert to long format
df_long <- df %>%
  pivot_longer(
    cols = -c(x, y),
    names_to = "date",
    values_to = "value"
  ) %>%
  separate(date, into = c("year", "month"), sep = "_") %>%
  mutate(
    year = as.integer(year),
    month = as.integer(month)
  ) %>%
  rename(
    lon = x,
    lat = y
  )

df2plot <- df_long %>%
  dplyr::filter((month > 6 & year == 2023) | year == 2024) %>%
  mutate(
    time = as.Date(paste0(year, "-", month, "-01"))
  ) %>%
  arrange(time)


world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  geom_tile(data = df2plot,
            aes(x = lon, y = lat, fill = value)) +
  geom_sf(data = world, fill = NA, color = "black", linewidth = 0.3) +
  scale_fill_gradient2(low = "darkred",mid = "white",high = "darkblue",limits = c(-3,3),
                       midpoint = 0, oob = scales::squish) +
  theme_bw() +
  facet_wrap(~time) +
  coord_sf(xlim = c(5, 35), ylim = c(-10, 5), expand = FALSE) +
  labs(
    fill = "Z-anomalies"
  )

plot(time(r),global(r,mean,na.rm = TRUE)[["mean"]],type = "l")

vals <- values(r)               # cells x time
dates <- time(r)
nT <- nlyr(r)

pos.thr <- 2
maxlag  <- 5

best.mag  <- rep(NA, nrow(vals))
best.time <- rep(NA, nrow(vals))

for(i in 1:nrow(vals)) {

  x <- vals[i, ]

  if(all(is.na(x))) next

  best_delta <- Inf
  best_date  <- NA

  for(t in 1:(nT-maxlag)) {

    if(!is.na(x[t]) && x[t] >= pos.thr) {

      fut <- x[(t+1):(t+maxlag)]

      if(any(!is.na(fut))) {

        j <- which.min(fut)   # driest future month
        delta <- fut[j] - x[t]

        if(delta < best_delta) {
          best_delta <- delta
          best_date  <- dates[t + j]
        }
      }
    }
  }

  if(is.finite(best_delta)) {
    best.mag[i]  <- best_delta
    best.time[i] <- year(best_date)
  }
}

# rasters
WhiplashMag  <- setValues(r[[1]], best.mag)
plot(WhiplashMag)

WhiplashDate <- setValues(r[[1]], best.time)
plot(WhiplashDate)

is.2024 <- (WhiplashDate > 2023)

plot(is.2024)

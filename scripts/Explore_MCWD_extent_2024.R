rm(list = ls())

library(terra)
library(lubridate)
library(rnaturalearth)
library(sf)
library(dplyr)
library(ggplot2)
library(tidyr)


system2("rsync",
        c("-avz",
          "hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/all.climate/*MCWD*",
          "./outputs/MCWD/"))

# ---------------------------
# 1. Inputs and masks
# ---------------------------

Extent <- ext(-15, 65, -15, 10)

Radar.data <- rast("./outputs/Radar_all.years.tif")
Radar.data <- crop(Radar.data[[4:nlyr(Radar.data)]], Extent)

IF <- crop(rast("./data/shapefiles/Intact_mask_mod.tif"), Extent)
IF <- rasterize(vect(read_sf("~/Documents/projects/Drying.CB/data/shapefiles/Rainforests.shp")), Radar.data, field = "LC")

IF.crop <- crop(IF, ext(Radar.data))
IF.crop[is.na(IF.crop)] <- 0
IF.r <- resample(IF.crop, Radar.data, method = "bilinear")
IF.r <- round(IF.r)

MCWD.allfiles <- list.files("./outputs/MCWD/", pattern = "\\.tif$", full.names = TRUE)
MCWD.files <- MCWD.allfiles[grepl("100mm", MCWD.allfiles)]
MCWD.files <- MCWD.allfiles[grepl("GLEAM", MCWD.allfiles)]

# ---------------------------
# 2. Read all products
# ---------------------------

df.ts <- data.frame()
all.MCWD <- list()

for (ifile in seq_along(MCWD.files)) {

  cMCWD.file <- MCWD.files[ifile]
  cproduct <- strsplit(basename(cMCWD.file), "_")[[1]][1]

  if (cproduct %in% c("CRUJRA","GLDAS")){
    next()
  }

  MCWD <- rast(cMCWD.file)
  MCWD <- crop(MCWD, ext(Radar.data))
  MCWD <- resample(MCWD, Radar.data, method = "bilinear")
  MCWD.IF <- mask(MCWD, IF.r, maskvalues = c(0, NA))

  all.MCWD[[cproduct]] <- MCWD.IF

  df.ts <- bind_rows(
    df.ts,
    data.frame(
      time = as.Date(time(MCWD.IF)),
      product = cproduct,
      MCWD.m = global(MCWD.IF, "mean", na.rm = TRUE)[, 1]
    )
  )
}

# Optional check
layer_info <- data.frame(
  product = names(all.MCWD),
  nlayers = sapply(all.MCWD, nlyr),
  start = sapply(all.MCWD, function(r) as.character(min(as.Date(time(r))))),
  end   = sapply(all.MCWD, function(r) as.character(max(as.Date(time(r)))))
)

print(layer_info)

ggplot(data = df.ts,
       aes(x = time, y = MCWD.m, color = product)) +
  geom_line() +
  theme_bw()


# ---------------------------
# 3. Restrict to dates common to all products
# ---------------------------

dates_list <- lapply(all.MCWD, function(r) as.Date(time(r)))
common_dates <- Reduce(intersect, dates_list)
common_dates <- sort(common_dates)

all.MCWD.common <- lapply(all.MCWD, function(r) {
  r[[match(common_dates, as.Date(time(r)))]]
})

print(sapply(all.MCWD.common, nlyr))

ntime <- length(common_dates)

# ---------------------------
# 4. Mean and SD maps across products for each common date
# ---------------------------

MCWD_mean <- rast(lapply(seq_len(ntime), function(i) {
  r_i <- rast(lapply(all.MCWD.common, function(r) r[[i]]))
  app(r_i, mean, na.rm = TRUE)
}))

MCWD_sd <- rast(lapply(seq_len(ntime), function(i) {
  r_i <- rast(lapply(all.MCWD.common, function(r) r[[i]]))
  app(r_i, sd, na.rm = TRUE)
}))

terra::time(MCWD_mean) <- common_dates
terra::time(MCWD_sd) <- common_dates

names(MCWD_mean) <- paste0("mean_", seq_len(ntime))
names(MCWD_sd) <- paste0("sd_", seq_len(ntime))

dates <- as.Date(time(MCWD_mean))

# ---------------------------
# 5. Monthly climatology excluding 2024
# ---------------------------

idx_ref <- which(year(dates) < 2024)
dates_ref <- dates[idx_ref]
mon_ref <- month(dates_ref)

MCWD_ref <- MCWD_mean[[idx_ref]]

MCWD_clim_mean <- tapp(MCWD_ref, index = mon_ref, fun = mean, na.rm = TRUE)
MCWD_clim_sd   <- tapp(MCWD_ref, index = mon_ref, fun = sd, na.rm = TRUE)

names(MCWD_clim_mean) <- month.abb
names(MCWD_clim_sd)   <- month.abb

# ---------------------------
# 6. 2024 monthly maps and anomalies
# ---------------------------

idx_2024 <- which(year(dates) >=2022)
dates_2024 <- dates[idx_2024]
mon_2024 <- month(dates_2024)

MCWD_2024 <- MCWD_mean[[idx_2024]]
terra::time(MCWD_2024) <- dates_2024
names(MCWD_2024) <- format(dates_2024, "%Y-%m")

MCWD_anom_2024 <- rast(lapply(seq_along(idx_2024), function(j) {
  MCWD_2024[[j]] - MCWD_clim_mean[[mon_2024[j]]]
}))

terra::time(MCWD_anom_2024) <- dates_2024
names(MCWD_anom_2024) <- format(dates_2024, "%Y-%m")

MCWD_z_2024 <- rast(lapply(seq_along(idx_2024), function(j) {
  (MCWD_2024[[j]] - MCWD_clim_mean[[mon_2024[j]]]) / MCWD_clim_sd[[mon_2024[j]]]
}))

terra::time(MCWD_z_2024) <- dates_2024
names(MCWD_z_2024) <- format(dates_2024, "%Y-%m")

# ---------------------------
# 7. Time series summaries
# ---------------------------

df_clim_ts <- data.frame(
  month = 1:12,
  clim = sapply(1:12, function(m) {
    global(MCWD_clim_mean[[m]], "mean", na.rm = TRUE)[1, 1]
  }),
  clim_sd = sapply(1:12, function(m) {
    global(MCWD_clim_sd[[m]], "mean", na.rm = TRUE)[1, 1]
  })
)

df_2024_ts <- data.frame(
  time = dates_2024,
  month = month(dates_2024),
  value_2024 = sapply(seq_along(idx_2024), function(j) {
    global(MCWD_2024[[j]], "mean", na.rm = TRUE)[1, 1]
  }),
  anomaly_2024 = sapply(seq_along(idx_2024), function(j) {
    global(MCWD_anom_2024[[j]], "mean", na.rm = TRUE)[1, 1]
  }),
  z_2024 = sapply(seq_along(idx_2024), function(j) {
    global(MCWD_z_2024[[j]], "mean", na.rm = TRUE)[1, 1]
  })
) %>%
  left_join(df_clim_ts, by = "month")

print(df_2024_ts)

# ---------------------------
# 8. Maps
# ---------------------------

world <- ne_countries(scale = "medium", returnclass = "sf")

# 2024 anomaly maps
df_anom_2024 <- as.data.frame(MCWD_anom_2024, xy = TRUE, na.rm = FALSE) %>%
  pivot_longer(
    cols = -c(x, y),
    names_to = "time",
    values_to = "value"
  )

p_anom <- ggplot() +
  geom_raster(
    data = df_anom_2024 %>%
      filter(month(paste0(time,"-01")) %in% c(2,5,8,11)) %>%
      filter(!is.na(value)),
    aes(x = x, y = y, fill = value)
  ) +
  geom_sf(data = world, fill = NA, color = "black", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "darkred", mid = "white", high = "darkblue",
    midpoint = 0, limits = c(-250, 250), oob = scales::squish
  ) +
  theme_bw() +
  facet_wrap(~ time) +
  coord_sf(xlim = c(-15, 65), ylim = c(-15, 10), expand = FALSE) +
  labs(fill = "MCWD anomaly", title = "2024 monthly MCWD anomaly vs climatology (excluding 2024)")

print(p_anom)

# Optional: standardized anomaly maps
df_z_2024 <- as.data.frame(MCWD_z_2024, xy = TRUE, na.rm = FALSE) %>%
  pivot_longer(
    cols = -c(x, y),
    names_to = "time",
    values_to = "value"
  )

p_z <- ggplot() +
  geom_raster(
    data = df_z_2024 %>%
      filter(month(paste0(time,"-01")) %in% c(2,5,8,11)) %>%
      filter(!is.na(value)),
    aes(x = x, y = y, fill = value)
  ) +
  geom_sf(data = world, fill = NA, color = "black", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "darkred", mid = "white", high = "darkblue",
    midpoint = 0, limits = c(-1, 1)*3, oob = scales::squish
  ) +
  theme_bw() +
  facet_wrap(~ time) +
  coord_sf(xlim = c(5, 35), ylim = c(-10, 5), expand = FALSE) +
  labs(fill = "z anomaly", title = "2024 standardized MCWD anomaly vs climatology (excluding 2024)")

print(p_z)

# ---------------------------
# 9. Time series plots
# ---------------------------

p_ts <- ggplot(df_2024_ts, aes(x = as.Date(paste0(time,"-01")))) +
  geom_line(aes(y = clim, group = 1), linewidth = 1) +
  geom_point(aes(y = clim), size = 2) +
  geom_line(aes(y = value_2024, group = 1), linewidth = 1, linetype = 2) +
  geom_point(aes(y = value_2024), size = 2) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Spatial mean MCWD",
    title = "2024 monthly MCWD vs climatology (excluding 2024)"
  )

print(p_ts)

p_anom_ts <- ggplot(df_2024_ts, aes(x = as.Date(paste0(time,"-01")), y = anomaly_2024)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line() +
  geom_point() +
  theme_bw() +
  labs(
    x = NULL,
    y = "Spatial mean MCWD anomaly",
    title = "2024 monthly MCWD anomaly vs climatology (excluding 2024)"
  )

print(p_anom_ts)

p_z_ts <- ggplot(df_2024_ts, aes(x = as.Date(paste0(time,"-01")), y = z_2024)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line() +
  geom_point() +
  theme_bw() +
  labs(
    x = NULL,
    y = "Spatial mean standardized anomaly",
    title = "2024 monthly standardized MCWD anomaly vs climatology (excluding 2024)"
  )

print(p_z_ts)

# ---------------------------
# 10. Optional annual comparison
# ---------------------------

MCWD_2024_annual <- mean(MCWD_2024[[year(time(MCWD_2024)) == 2024]], na.rm = TRUE)
MCWD_clim_annual <- mean(MCWD_clim_mean, na.rm = TRUE)
MCWD_2024_annual_anom <- MCWD_2024_annual - MCWD_clim_annual

print(global(MCWD_2024_annual, "mean", na.rm = TRUE))
print(global(MCWD_clim_annual, "mean", na.rm = TRUE))
print(global(MCWD_2024_annual_anom, "mean", na.rm = TRUE))


writeRaster(MCWD_mean,
            paste0("./outputs/MCWD.MEM.tif"),
            overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

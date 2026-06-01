rm(list = ls())

library(terra)
library(lubridate)
library(rnaturalearth)
library(sf)
library(dplyr)
library(ggplot2)
library(slider)

Extent <-  ext(-15,65,-15,10)

Radar.data <- rast("./outputs/Radar_all.years.tif")
Radar.data <- crop(Radar.data[[4:nlyr(Radar.data)]],
                   Extent)

# IF <- rast("./data/shapefiles/Peatlands_mask.tif")
# IF <- crop(rast("~/Documents/projects/CausalAI/data/Intact_forest_mask.tif"),
#            Extent)
# IF <- rasterize(vect(read_sf("~/Documents/projects/Drying.CB/data/shapefiles/Rainforests.shp")), Radar.data, field = "LC")
# IF <- crop(rast("./data/shapefiles/Intact_mask.tif"),
#            Extent)

IF <- crop(rast("./data/shapefiles/Intact_mask_mod.tif"),
           Extent)

# IF <- rasterize(vect(read_sf("~/Documents/projects/Drying.CB/data/shapefiles/CongoBasin.shp")), Radar.data, field = "BASIN_ID")
# IF <- rast("~/Documents/projects/CausalAI/data/Intact_forest_mask.tif")**0
# IF <- rasterize(as.polygons(ext(5, 35, -10, 5)), Radar.data[[1]], field = 1, background = 0)
# IF <- rasterize(as.polygons(ext(20, 25, 0, 2.5)), Radar.data[[1]], field = 1, background = 0)
# IF <- rasterize(as.polygons(ext(18,25,-2,4)), Radar.data[[1]], field = 1, background = 0)


plot(IF)


IF.crop <- crop(IF, ext(Radar.data))
IF.crop[is.na(IF.crop)] <- 0
IF.r <- resample(IF.crop, Radar.data, method = "bilinear")
IF.r <- round(IF.r)

plot(IF.crop)
plot(IF.r)

MCWD <- resample(crop(rast("~/Downloads/chirps_MCWD_rolling12m_GLEAMPET_all.months.tif"),
                      ext(Radar.data)),
                 Radar.data, method = "bilinear")
MCWD.IF <- mask(MCWD,IF.r,maskvalues = c(0,NA))
common.dates <- which(time(MCWD) %in% time(Radar.data))
MCWD.IF.common <- MCWD.IF[[common.dates]]

MCWD2plot <- MCWD.IF.common[[which(time(MCWD.IF.common) == as.Date("2024-08-01"))]]
plot(MCWD2plot)
hist(MCWD2plot)
summary(MCWD2plot)

Radar.data.IF <- mask(Radar.data,IF.r,maskvalues = c(0,NA))
pos <- which(time(Radar.data.IF) < as.Date("2024/09/01"))
Radar.data.IF <-  Radar.data.IF[[pos]]
Radar.data <- Radar.data[[pos]]

plot(Radar.data.IF[[length(Radar.data.IF)]])

mths <- as.integer(format(time(Radar.data.IF), "%m"))
clim <- tapp(Radar.data.IF, index = mths, fun = mean, na.rm = TRUE)
clim.sd <- tapp(Radar.data.IF, index = mths, fun = sd, na.rm = TRUE)

Radar.anom <- Radar.data.IF - clim[[mths]]
terra::time(Radar.anom) <- time(Radar.data.IF)
global(Radar.anom[[383]],mean,na.rm = TRUE)

Radar.anom.norm <- (Radar.data.IF - clim[[mths]]) / clim.sd[[mths]]
terra::time(Radar.anom.norm) <- time(Radar.data.IF)

Nsmooth = 9

ts.Radar <- data.frame(
  time = time(Radar.data.IF),
  MCWD = as.vector(global(MCWD.IF.common[[pos]],mean,na.rm = TRUE))[["mean"]],
  value = as.vector(global(Radar.data.IF,mean,na.rm = TRUE))[["mean"]],
  anomaly.av = as.vector(global(Radar.anom,mean,na.rm = TRUE))[["mean"]],
  anomaly.norm.av = as.vector(global(Radar.anom.norm,mean,na.rm = TRUE))[["mean"]]) %>%
  mutate(year = year(time),
         month = month(time)) %>%
  group_by(month) %>%
  mutate(clim_month = mean(value, na.rm = TRUE),
         sd_month = sd(value, na.rm = TRUE),
         anomaly = value - clim_month,
         anomaly.norm = (value - clim_month)/sd_month) %>%
  ungroup() %>%
  arrange(time) %>%
  mutate(
    MCWD_ma12 = slide_dbl(MCWD, mean, .before = (Nsmooth-1), .complete = TRUE, na.rm = TRUE),
    value_ma12 = slide_dbl(value, mean, .before = (Nsmooth-1), .complete = TRUE, na.rm = TRUE),
    anomaly_ma12 = slide_dbl(anomaly, mean, .before = (Nsmooth-1), .complete = TRUE, na.rm = TRUE),
    anomaly.norm_ma12 = slide_dbl(anomaly.norm, mean, .before = (Nsmooth-1), .complete = TRUE, na.rm = TRUE),
    anomaly.av_ma12 = slide_dbl(anomaly.av, mean, .before = (Nsmooth-1), .complete = TRUE, na.rm = TRUE),
    anomaly.norm.av_ma12 = slide_dbl(anomaly.norm.av, mean, .before = (Nsmooth-1), .complete = TRUE, na.rm = TRUE))

plot(ts.Radar %>%
       filter(year %in% c(2021:2024)) %>%
       pull(value),type = "l")

lines(ts.Radar %>%
       filter(year %in% c(2021:2024)) %>%
       pull(clim_month),type = "l",lty = 2)

abline(v = c(12,24,36) + 0.5, col = "red")

plot(ts.Radar$anomaly.av,type = 'l')


ggplot(data = ts.Radar,
       aes(x = time,
           y = MCWD)) +
  geom_rect(inherit.aes = FALSE,
            aes(xmin = as.Date("2005-01-01"),xmax = as.Date("2007-01-01"),
                ymin = -Inf,ymax = Inf),
            fill = "grey",
            alpha = 0.2) +
  geom_rect(inherit.aes = FALSE,
            aes(xmin = as.Date("2015-06-01"),xmax = as.Date("2017-06-01"),
                ymin = -Inf,ymax = Inf),
            fill = "grey",
            alpha = 0.2) +
  geom_rect(inherit.aes = FALSE,
            aes(xmin = as.Date("2023-01-01"),xmax = as.Date("2025-01-01"),
                ymin = -Inf,ymax = Inf),
            fill = "grey",
            alpha = 0.2) +
  geom_line(color = "darkgrey") +
  geom_line(aes(x = time,
                y = MCWD_ma12)) +
  stat_smooth(method = "lm", color = "black") +
  theme_bw()

ggplot(data = ts.Radar,
       aes(x = time,
           y = anomaly.av)) +
  geom_rect(inherit.aes = FALSE,
            aes(xmin = as.Date("2015-06-01"),xmax = as.Date("2017-06-01"),
                ymin = -Inf,ymax = Inf),
            fill = "grey",
            alpha = 0.2) +
  geom_rect(inherit.aes = FALSE,
            aes(xmin = as.Date("2024-01-01"),xmax = as.Date("2025-01-01"),
                ymin = -Inf,ymax = Inf),
            fill = "grey",
            alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(color = "darkgrey") +
  geom_line(aes(x = time,
                y = anomaly.av_ma12)) +
  stat_smooth(method = "lm", color = "black") +
  theme_bw() + labs(x = "",y = "C-band anomaly (dB)") +
  theme(text = element_text(size = 20))


summary(lm(data = ts.Radar,
           formula = anomaly.av_ma12 ~ year))

t <- time(Radar.anom)
t_num <- as.numeric(format(t, "%Y")) + (as.numeric(format(t, "%m")) - 1) / 12

# function to compute slope
fun_slope <- function(x) {
  if (all(is.na(x))) return(NA)
  coef(lm(x ~ t_num))[2]
}

Radar.slope <- app(Radar.anom, fun = fun_slope)
df.slope <- as.data.frame(Radar.slope, xy = TRUE, na.rm = TRUE) %>%
  rename(slope = lyr.1)

world <- ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  geom_raster(data = df.slope,
              aes(x = x, y = y, fill = slope)) +
  geom_sf(data = world, fill = NA, color = "black", linewidth = 0.3) +
  scale_fill_gradient2(low = "darkred",mid = "white",high = "darkblue",
                       midpoint = 0,
                       limits = c(-1,1)/100, oob = scales::squish) +
  theme_bw() +
  coord_sf(xlim = c(5, 35), ylim = c(-10, 5), expand = FALSE) +
  labs(
    fill = "Trend (per year)",
    x = "Longitude",
    y = "Latitude"
  )

hist(df.slope$slope)

mat <- values(Radar.anom, mat = TRUE)
xy <- as.data.frame(Radar.anom[[1]],xy = TRUE,na.rm = FALSE)

df.map <- data.frame(
  lon = rep(xy$x, times = nlyr(Radar.data)),
  lat = rep(xy$y, times = nlyr(Radar.data)),
  time = rep(time(Radar.data), each = nrow(mat)),
  value = as.vector(mat)
) %>%
  filter(!is.na(value))

ggplot() +
  geom_raster(data = df.map %>%
                filter(year(time) == 2024,
                       month(time) %in% 1:9),
              aes(x = lon, y = lat, fill = value)) +
  geom_sf(data = world, fill = NA, color = "black", linewidth = 0.3) +
  scale_fill_gradient2(low = "darkred",mid = "white",high = "darkblue",
                       midpoint = 0, limits = c(-0.5,0.5), oob = scales::squish) +
  theme_bw() +
  labs(x = "", y = "", fill = "Anomaly (dB)") +
  facet_wrap(~month(time)) +
  coord_sf(xlim = c(5, 35), ylim = c(-10, 5), expand = FALSE)


df_timing <- df.map %>%
  group_by(lon,lat) %>%
  summarise(time = year(time[which.min(value)]),
            .groups = "keep") %>%
  mutate(timing = case_when(time %in% c(2023:2024)~ "now",
                            TRUE ~ "then"))


ggplot() +
  geom_raster(data = df_timing,
              aes(x = lon, y = lat, fill = timing)) +
  geom_sf(data = world, fill = NA, color = "black", linewidth = 0.3) +
  theme_bw() +
  labs(x = "", y = "", fill = "Anomaly (dB)") +
  coord_sf(xlim = c(5, 35), ylim = c(-10, 5), expand = FALSE)

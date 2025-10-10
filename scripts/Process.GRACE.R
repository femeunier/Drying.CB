rm(list = ls())

library(terra)
library(ncdf4)
library(ggplot2)
library(lubridate)

# nc <- nc_open("~/Downloads/GRAVIS-3_GFZOP_0600_TWS_GRID_GFZ_0006.nc")
# times <- as.Date(ncvar_get(nc,"time"),
#                  origin = "2002/04/18")
# nc_close(nc)


R <- rast("~/Downloads/GRAVIS-3_GFZOP_0600_TWS_GRID_GFZ_0006.nc")
R2 <- rotate(R)  # moves extent to [-180, 180]

R2.crop <- crop(R2,
                ext(-180,180,-25,25))

plot(R2.crop[[1]])

variable_names <- c("tws","std_tws","leakage","model_atmosphere")
N <- nlyr(R2.crop)/4

for (i in seq(1,length(variable_names))){

  pos <- ((i-1)*N + 1) : ((i)*N)
  crast <- R2.crop[[pos]]
  cvar <- variable_names[i]

  print(cvar)

  writeRaster(crast,
              paste0("./data/GRACE/",cvar,"_all.years.tif"),
              overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))


}

B <- rast("./data/GRACE/tws_all.years.tif")
time(B)
plot(B[[235]])


B.cropped <- crop(B,ext(-25,65,-25,25))
ts_df <- global(B.cropped, fun = mean, na.rm = TRUE)
ts_df$date <- time(B.cropped)


ts.mod <- ts_df %>%
  mutate(year = year(date),
         month = month(date)) %>%
  mutate(period = case_when(date < "2019/08/01" ~ 1,
                            TRUE ~ 2)) %>%
  group_by(period) %>%
  mutate(m = mean(mean)) %>%
  mutate(corrected = mean - m)

ggplot(data = ts.mod) +
  geom_line(aes(x = year + (month - 1/2)/12,
                y = mean)) +
  # scale_x_continuous(limits = c(2019,2021)) +
  theme_bw()


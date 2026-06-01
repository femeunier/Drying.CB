rm(list = ls())

library(readr)
library(stringr)
library(ncdf4)
library(terra)
library(lubridate)
library(dplyr)
library(ggplot2)
library(sf)
library(zoo)


read_SWAT_file <- function(file,
                           header_line_num = 2,
                           data_start_line = 4,
                           skip = 3,
                           stringsAsFactors = FALSE) {

  lines <- readLines(file, warn = FALSE)

  header_line <- lines[header_line_num]

  # Column names from header
  col_names <- strsplit(trimws(header_line), "\\s+")[[1]]

  # Start positions of columns in the header
  starts <- gregexpr("\\S+", header_line)[[1]]

  # Use maximum line length, not header length
  max_width <- max(nchar(lines))

  # Column widths
  widths <- c(
    diff(starts),
    max_width - tail(starts, 1) + 1
  )

  # Read fixed-width file
  dat <- read.fwf(
    file,
    widths = widths,
    skip = skip,
    header = FALSE,
    stringsAsFactors = FALSE,
    strip.white = TRUE
  )

  names(dat) <- col_names

  return(dat)

}


file <- "~/Downloads/africa-congo 2/Scenarios/Default/TxtInOut/basin_pw_mon.txt"
plant.dat <- read_SWAT_file(file) %>%
  group_by(mon) %>%
  mutate(clim.tmn.m = mean(tmn,na.rm = TRUE),
         clim.tmn.sd = sd(tmn,na.rm = TRUE),

         clim.pr.m = mean(percn,na.rm = TRUE),
         clim.pr.sd = sd(percn,na.rm = TRUE)) %>%
  mutate(Zscore.tmn = (tmn - clim.tmn.m)/clim.tmn.sd,
         Zscore.pr = (percn - clim.pr.m)/clim.pr.sd)


ggplot(data = plant.dat) +
  geom_rect(aes(xmin = 2024, xmax = 2025,
                ymin = -Inf, ymax = Inf),
            fill = "grey80",
            alpha = 0.5,
            inherit.aes = FALSE) +
  geom_hline(yintercept = 0,linetype = 2) +
  geom_line(aes(x = yr + (mon - 1/2)/12,
                y = Zscore.pr)) +
  theme_bw()



file <- "~/Downloads/africa-congo 2/Scenarios/Default/TxtInOut/basin_wb_mon.txt"
wb.dat <- read_SWAT_file(file) %>%
  group_by(mon) %>%
  mutate(surq_runon.m = mean(surq_runon,na.rm = TRUE),
         surq_runon.sd = sd(surq_runon,na.rm = TRUE)) %>%
  mutate(Zscore.surq_runon = (surq_runon - surq_runon.m)/surq_runon.sd)


ggplot(data = wb.dat) +
  geom_rect(aes(xmin = 2024, xmax = 2025,
                ymin = -Inf, ymax = Inf),
            fill = "grey80",
            alpha = 0.5,
            inherit.aes = FALSE) +
  geom_hline(yintercept = 0,linetype = 2) +
  geom_line(aes(x = yr + (mon - 1/2)/12,
                y = Zscore.surq_runon)) +
  theme_bw()


###############################################################################

file <- "~/Downloads/africa-congo 2/Scenarios/Default/TxtInOut/hru_wb_mon.txt"
wb.dat <- read_SWAT_file(file)


###############################################################################

CB <- read_sf("./data/shapefiles/CongoBasin.shp")
CB <- read_sf("/Users/felicien/Downloads/africa-congo 2/Scenarios/Default/Results/deep_aquifers.shp")

f <- "/Users/felicien/Downloads/africa-congo 2/Scenarios/Default/TxtInOut/20crv3-era5.nc4"

# Load variables
tmin <- rast(f, subds = "tmin")
tmax <- rast(f, subds = "tmax")
pcp  <- rast(f, subds = "pcp")

# Extract dates
nc <- nc_open(f)
time <- ncvar_get(nc, "time")
dates <- as.Date(time, origin = "1900-01-01")
nc_close(nc)

terra::time(tmin) <- dates
terra::time(tmax) <- dates
terra::time(pcp)  <- dates

# Daily mean temperature
tas <- (tmin + tmax)/2
names(tas) <- format(dates, "%Y_%m_%d")
names(pcp) <- format(dates, "%Y_%m_%d")

tas <- (tmin + tmax) / 2
terra::time(tas) <- dates

# Monthly aggregation
month_id <- format(dates, "%Y-%m")

tas.mon <- tapp(tas, month_id, mean, na.rm = TRUE)
pcp.mon <- tapp(pcp, month_id, sum,  na.rm = TRUE)  # precip: monthly totals

# Add monthly dates
mon_dates <- as.Date(paste0(month_id[!duplicated(month_id)], "-15"))

terra::time(tas.mon) <- mon_dates
terra::time(pcp.mon) <- mon_dates

names(tas.mon) <- format(mon_dates, "%Y_%m")
names(pcp.mon) <- format(mon_dates, "%Y_%m")

tas.mon <- tas.mon[[year(mon_dates) >= 2001]]
pcp.mon <- pcp.mon[[year(mon_dates) >= 2001]]

zscore_monthly <- function(x) {

    mons <- month(terra::time(x))

    out <- x

    for(m in 1:12){

      idx <- which(mons == m)
      xm  <- x[[idx]]

      clim_mean <- app(xm, mean, na.rm = TRUE)
      clim_sd   <- app(xm, sd, na.rm = TRUE)

      out[[idx]] <- (xm - clim_mean) / clim_sd
    }

    out
}

CB2 <- vect(st_transform(CB, crs(tas.mon)))

tas.z <- mask(zscore_monthly(tas.mon),CB2)
pcp.z <- mask(zscore_monthly(pcp.mon),CB2)

plot(tas.z[[grepl("2024",names(tas.z))]])
plot(pcp.z[[grepl("2024",names(pcp.z))]])

ts.tas.z <- global(tas.z,mean,na.rm = TRUE)[["mean"]]
ts.pr.z <- global(pcp.z,mean,na.rm = TRUE)[["mean"]]

tas.rm12 <- rollmean(ts.tas.z, k = 12, fill = NA, align = "right")
pr.rs12 <- rollmean(ts.pr.z, k = 12, fill = NA, align = "right")


plot(time(tas.z),tas.rm12,type = "l",ylim = c(-3,3))
lines(time(tas.z),ts.tas.z,type = "l",lwd = 0.2)
abline(h = 0, lty = 2)

plot(time(pcp.z),pr.rs12,type = "l",ylim = c(-1,1)*2)
lines(time(tas.z),ts.pr.z,type = "l",lwd = 0.2)
abline(h = 0, lty = 2)


plant.file <- "~/Downloads/africa-congo 2/Scenarios/Default/TxtInOut/basin_pw_mon.txt"
plant.dat <- read_SWAT_file(plant.file)

water.file <- "~/Downloads/africa-congo 2/Scenarios/Default/TxtInOut/basin_wb_mon.txt"
water.dat <- read_SWAT_file(water.file)


plot(time(tas.z),global(mask(tas.mon,CB2),mean,na.rm = TRUE)[["mean"]],
     type = "l", ylim = c(22,27))
lines(as.Date(paste0(plant.dat$yr,"/",plant.dat$mon,"/01")),
      plant.dat$tmpav, col = "red")


plot(time(tas.z),global(mask(pcp.mon,CB2),mean,na.rm = TRUE)[["mean"]],
     type = "l",ylim = c(0,300))
lines(as.Date(paste0(water.dat$yr,"/",water.dat$mon,"/01")),
      water.dat$precip, col = "red")


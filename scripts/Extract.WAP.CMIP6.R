rm(list = ls())

library(terra)
library(lubridate)
library(tidyr)
library(slider)
library(dplyr)
library(ncdf4)
library(reshape2)
library(ncdf4.helpers)

CATL <- ext(-23,3,-2.4,5.5)
CC <- ext(22,29,-10,4)
variants <- c("r1i1p1f1")

main.dir <- '/data/gent/vo/000/gvo00074/felicien/CMIP6'

scenarios <- c("historical","ssp126",
               'ssp245',"ssp370","ssp585")

vars <- c("pr","wap","tos")

years <- 1850:2100

df.all <- data.frame()

r_template <- rast(xmin = -180, xmax = 180,
                   ymin = -90,  ymax = 90,
                   resolution = 2,
                   crs = "EPSG:4326")

all.rast <- list() ; raster.names <- c() ; compt <- 1

for (cvar in vars){
  for (cscenario in scenarios){

    cdir <- file.path(main.dir,cscenario,cvar)
    all.files <- list.files(cdir,
                            pattern = "*.nc$",
                            full.names = TRUE)

    if (cvar %in% c("wap","pr")){
      all.files <- all.files[grepl("Amon",
                                   all.files)]
    } else if (cvar == "tos"){
      all.files <- all.files[grepl("Omon",
                                   all.files)]
    }

    all.files <- all.files[grepl(paste0(variants,collapse = "|"),
                                 all.files)]

    # all.files <- all.files[grepl("MCM-UA-1-0",all.files)]

    if (length(all.files) == 0) {
      next
    }

    all.models <- sapply(strsplit(tools::file_path_sans_ext(basename(all.files)),
                                  "\\_"),"[[",3)
    models <- sort(unique(sapply(strsplit(tools::file_path_sans_ext(basename(all.files)),
                              "\\_"),"[[",3)))

    for (cmodel in models){

      cfiles <- all.files[all.models == cmodel]

      for (cfile in cfiles){
        csplit <- strsplit(tools::file_path_sans_ext(basename(cfile)),
                           "\\_")[[1]]
        cVar <- csplit[1]
        cModel <- cmodel

        cScenario <- csplit[4]
        cVariant <- csplit[5]
        cDates <- csplit[7]

        cyear.min <- as.numeric(substr(cDates,1,4))
        cyear.max <- as.numeric(substr(cDates,8,11))

        if (!(any(cyear.min:cyear.max %in% years))) next

        print(
          paste0(cScenario," - ",cVar," - ",cModel," - ",cDates))

        cR <- tryCatch(rast(cfile),
                       error = function(e) NULL)

        if (is.null(cR)){
          next
        }
        Extent <- as.numeric(as.vector(ext(cR)))

        if (Extent[1] == Extent[3]){

          nc  <- nc_open(cfile)
          values <- ncvar_get(nc, cVar)

          if (cVar == "wap" && length(dim(values)) == 3) {

            lon1d <- ncvar_get(nc, "longitude")  # length = i
            lat1d <- ncvar_get(nc, "latitude")   # length = i
            plev  <- ncvar_get(nc, "plev")       # Pa, length = 39
            dates <- nc.get.time.series(nc, v = cVar)

            fill_values <- ncatt_get(nc, cVar, "_FillValue")$value
            nc_close(nc)

            # choose nearest level to 50000 Pa (500 hPa)
            pos <- which.min(abs(plev - 50000))
            sel_plev <- plev[pos]
            if (sel_plev != 50000) {
              warning(paste("50000 Pa not available, using", sel_plev, "Pa instead"))
            }

            # reduce to [i, time]
            vals_it <- values[, pos, ]   # dims: i x time
            vals_it[vals_it >= fill_values/10] <- NA

            nt <- dim(vals_it)[2]
            mcoord <- is.finite(lon1d) & is.finite(lat1d)

            r_list <- vector("list", nt)
            for (t in seq_len(nt)) {
              v <- vals_it[, t]
              m <- mcoord & is.finite(v)

              df <- data.frame(
                lon = lon1d[m],
                lat = lat1d[m],
                val = v[m]
              )

              pts <- vect(df, geom = c("lon","lat"), crs = "EPSG:4326")
              r_list[[t]] <- rasterize(pts, r_template, field = "val", fun = "mean")
            }

            cR <- rast(r_list)
            names(cR) <- paste0("wap_", format(dates, "%Y%m%d"))
            crs(cR) <- "EPSG:4326"
            time(cR) <- as.Date(as.character(dates))

          } else if (cVar == "pr" && length(dim(values)) == 2){

            lon1d <- ncvar_get(nc, "longitude")  # length = i
            lat1d <- ncvar_get(nc, "latitude")   # length = i
            dates <- nc.get.time.series(nc, v = cVar)

            fill_values <- ncatt_get(nc, cVar, "_FillValue")$value
            nc_close(nc)

            # reduce to [i, time]
            vals_it <- values[, ]   # dims: i x time
            vals_it[vals_it >= fill_values/10] <- NA

            nt <- dim(vals_it)[2]
            mcoord <- is.finite(lon1d) & is.finite(lat1d)

            r_list <- vector("list", nt)
            for (t in seq_len(nt)) {
              v <- vals_it[, t]
              m <- mcoord & is.finite(v)

              df <- data.frame(
                lon = lon1d[m],
                lat = lat1d[m],
                val = v[m]
              )

              pts <- vect(df, geom = c("lon","lat"), crs = "EPSG:4326")
              r_list[[t]] <- rasterize(pts, r_template, field = "val", fun = "mean")
            }

            cR <- rast(r_list)
            names(cR) <- paste0("pr_", format(dates, "%Y%m%d"))
            crs(cR) <- "EPSG:4326"
            time(cR) <- as.Date(as.character(dates))


          } else if (cVar == "tos" & length(dim(values)) == 3){
            lon2d <- ncvar_get(nc,
                               intersect(c("lon", "longitude","nav_lon"), names(nc$var))[1])
            lon2d[lon2d>180] <-
              lon2d[lon2d>180] -360
            lat2d <- ncvar_get(nc,
                               intersect(c("lat", "latitude",
                                           "nav_lat"), names(nc$var))[1])


            fill_values <- ncatt_get(nc, cVar, "_FillValue")$value
            mv_lon <- ncatt_get(nc, intersect(c("lon", "longitude",
                                                "nav_lon"), names(nc$var))[1], "missing_value")$value
            mv_lat <- ncatt_get(nc, intersect(c("lat", "latitude",
                                                "nav_lat"), names(nc$var))[1], "missing_value")$value

            times <- ncvar_get(nc,"time")
            dates <- nc.get.time.series(nc, v = cVar)

            nc_close(nc)

            lon2d[lon2d == mv_lon] <- NA
            lat2d[lat2d == mv_lat] <- NA

            ni <- dim(lon2d)[1]
            nj <- dim(lon2d)[2]
            nt <- dim(values)[3]

            r_list <- vector("list", nt)

            for (t in seq_len(nt)) {
              vals <- values[ , , t]

              # mask invalid SSTs & coords
              vals[vals >= fill_values/10] <- NA
              m <- is.finite(vals) & is.finite(lon2d) & is.finite(lat2d)

              df <- data.frame(
                lon = as.vector(lon2d[m]),
                lat = as.vector(lat2d[m]),
                tos = as.vector(vals[m])
              )

              pts <- vect(df,
                          geom = c("lon", "lat"), crs = "EPSG:4326")

              # rasterize to regular grid
              r_t <- rasterize(pts, r_template, field = "tos", fun = "mean")
              r_list[[t]] <- r_t
            }

            # stack into one SpatRaster
            cR <- rast(r_list)
            names(cR) <- paste0(cVar, seq_len(nt))
            crs(cR)   <- "EPSG:4326"   # already in lon/lat

            time(cR) <- as.Date(as.character(dates))
          } else if (cVar == "tos") {

            lons  <- ncvar_get(nc,
                               intersect(c("lon", "longitude"), names(nc$var))[1])
            lats  <- ncvar_get(nc,
                               intersect(c("lat", "latitude"), names(nc$var))[1])
            times <- ncvar_get(nc,"time")
            dates <- nc.get.time.series(nc, v = cVar)

            nc_close(nc)

            n_time <- length(times)
            r_list <- vector("list", n_time)

            for (i in seq_len(n_time)) {
              # build a data.frame for this time slice
              df <- data.frame(lon = lons,
                               lat = lats,
                               val = values[, i])

              # convert to SpatVector
              pts <- vect(df, geom = c("lon", "lat"), crs = "EPSG:4326")

              # rasterize onto template (binning: mean of points per cell)
              r_i <- rasterize(pts, r_template, field = "val", fun = "mean")
              r_list[[i]] <- r_i
            }

            # stack into a single SpatRaster with nlyr = n_time
            cR <- rast(r_list)

            # name layers and attach time
            names(cR) <- paste0(cVar, "_", format(dates, "%Y%m%d"))
            time(cR)  <- as.Date(as.character(dates))    # terra time dimension

          }
        }

        if (max(as.vector(ext(cR))) >= 350){
          cR <- rotate(cR)
        }

        if (max(as.vector(ext(cR))) >= 200){
          nc  <- nc_open(cfile)
          lon2d <- ncvar_get(nc, "longitude")  # [i,j] = [320,384]
          lat2d <- ncvar_get(nc, "latitude")   # [i,j]
          nc_close(nc)

          lon1d <- apply(lon2d, 1, mean, na.rm = TRUE)   # length 320
          lat1d <- apply(lat2d, 2, mean, na.rm = TRUE)   # length 384

          ext(cR) <- ext(
            min(lon1d, na.rm = TRUE),
            max(lon1d, na.rm = TRUE),
            min(lat1d, na.rm = TRUE),
            max(lat1d, na.rm = TRUE)
          )

          crs(cR) <- "EPSG:4326"   # lon/lat WGS84

          if (max(as.vector(ext(cR))) >= 350){
            cR <- rotate(cR)
          }

          time(cR) <- as.Date(as.character(dates))

        }

        cR.CATL <- crop(cR,CATL)
        cR.CC <- crop(cR,CC)

        if (cvar == "pr"){
          raw <- bind_rows(data.frame(file = basename(cfile),
                                      time = time(cR),
                                      var = cVar,
                                      model = cModel,
                                      scenario = cScenario,
                                      variant = cVariant,
                                      region = "CC",
                                      value = as.numeric(
                                        global(cR.CC,mean,na.rm = TRUE)[["mean"]])))
        } else if (cvar == "tos"){
          raw <- bind_rows(data.frame(file = basename(cfile),
                                      time = time(cR),
                                      var = cVar,
                                      model = cModel,
                                      scenario = cScenario,
                                      variant = cVariant,
                                      region = "CATL",
                                      value = as.numeric(
                                        global(cR.CATL,mean,na.rm = TRUE)[["mean"]])))
        } else { # wap

          if (!is.null(depth(cR))) {
            depths <- depth(cR)
            pos <- which(depths == 50000)

            if (length(pos) == 0){
              uni.depths <- unique(depths)
              cdepth <- uni.depths[which.min(abs(uni.depths - 50000))]
              warning(paste("50000 not available, using", cdepth, "instead"))
              pos <- which(depths == cdepth)
            }

            cR <- cR[[pos]]
          }

          cR.CATL <- crop(cR,CATL)
          cR.CC   <- crop(cR,CC)

          raw <- bind_rows(
            data.frame(file = basename(cfile),
                       time = time(cR),
                       var = cVar,
                       model = cModel,
                       scenario = cScenario,
                       variant = cVariant,
                       region = "CATL",
                       value = as.numeric(global(cR.CATL, mean, na.rm = TRUE)[["mean"]])),
            data.frame(file = basename(cfile),
                       time = time(cR),
                       var = cVar,
                       model = cModel,
                       scenario = cScenario,
                       variant = cVariant,
                       region = "CC",
                       value = as.numeric(global(cR.CC, mean, na.rm = TRUE)[["mean"]])))

        }

        if (all(is.na(raw$time))){
          new.times <- as.Date(seq(as.Date(paste0(substr(cDates,1,4),"/",substr(cDates,5,6),"/01")),
              as.Date(paste0(substr(cDates,8,11),"/",substr(cDates,12,13),"/01")),
              by = "month"))
          if (length(new.times) != nrow(raw)){
            new.times <- as.Date(seq(as.Date(paste0(substr(cDates,1,4),"/",substr(cDates,5,6),"/01")),
                                     31+as.Date(paste0(substr(cDates,8,11),"/",substr(cDates,12,13),"/01")),
                                     by = "month"))
          }

          raw$time <- new.times

        }


        df.all <- bind_rows(df.all,
                            raw %>%
                              filter(year(time) %in% years))

        if (cfile == cfiles[1]){

           all.rast[[compt]] <- cR[[1]]
           raster.names[compt] <- paste0(cmodel,".",
                                         cScenario,".",
                                         cVariant,".",
                                         cVar,".",
                                         cDates)
           compt <- compt + 1

        }
      }
    }
  }
}

saveRDS(df.all,
        "./outputs/CMIP6.summary.RDS")

Scollection <- sprc(all.rast)
names(Scollection) <- raster.names
saveRDS(Scollection,
        "./outputs/CMIP6.checks.RDS")

# scp /Users/felicien/Documents/projects/Drying.CB/scripts/Extract.WAP.CMIP6.R hpc:/data/gent/vo/000/gvo00074/felicien/R/


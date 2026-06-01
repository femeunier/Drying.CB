rm(list = ls())

library(terra)
library(lubridate)
library(tidyr)
library(slider)
library(dplyr)
library(ncdf4)
library(reshape2)
library(sf)
library(ncdf4.helpers)


regions <- list(CATL = ext(-23,3,-2.4,5.5),
                CC = ext(22,29,-10,4),
                CB = read_sf("./data/CongoBasin.shp"),
                rainforest = read_sf("./data/Rainforests.shp"))

variants <- c("r1i1p1f1")

main.dir <- '/data/gent/vo/000/gvo00074/felicien/CMIP6'

scenarios <- c("historical","ssp126",
               'ssp245',"ssp370","ssp585")

vars <- c("pr")

years <- 1850:2100

df.all <- data.frame()


r_template <- rast(xmin = -180, xmax = 180,
                   ymin = -90,  ymax = 90,
                   resolution = 2,
                   crs = "EPSG:4326")

for (cvar in vars){
  for (cscenario in scenarios){

    cdir <- file.path(main.dir,cscenario,cvar)
    all.files <- list.files(cdir,
                            pattern = "*.nc$",
                            full.names = TRUE)
    all.files <- all.files[grepl("Amon",
                                 all.files)]
    all.files <- all.files[grepl(paste0(variants,collapse = "|"),
                                 all.files)]

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

          if(cVar == "pr" && length(dim(values)) == 2){

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


        for (ishape in seq(1,length(regions))){
          raw <- bind_rows(data.frame(file = basename(cfile),
                                      time = time(cR),
                                      var = cVar,
                                      model = cModel,
                                      scenario = cScenario,
                                      variant = cVariant,
                                      region = names(regions)[ishape],
                                      value = as.numeric(
                                        global(crop(cR,regions[[ishape]]),
                                               mean,
                                               na.rm = TRUE)[["mean"]])))

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


        }

      }
    }
  }
}

saveRDS(df.all,
        "./outputs/CMIP6.precip.RDS")


# scp /Users/felicien/Documents/projects/Drying.CB/scripts/Extract.precip.CMIP6.R hpc:/data/gent/vo/000/gvo00074/felicien/R/


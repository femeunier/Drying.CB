rm(list = ls())

library(dplyr)
library(terra)
library(lubridate)

files <- c("./outputs/All.climatevars.CA.RDS",
           "./outputs/All.CMIP6.vars.CA.RDS",
           "./outputs/All.GLEAM.CA.RDS",
           "./outputs/All.GPP.CA.RDS",
           "./outputs/All.GRACE.CA.RDS",
           "./outputs/All.Radar.CA.RDS",
           "./outputs/All.MODIS.CA.RDS")

types <- c("Observational",
           "CMIP6",
           "Observational",
           "Observational",
           "Observational",
           "Observational",
           "Observational")

all.ts <- data.frame()
for (ifile in seq(1,length(files))){

  cfile <- files[ifile]
  ctype <- types[ifile]

  print(cfile)

  cdata <- readRDS(cfile)

  if (!("scenario" %in% colnames(cdata))){
    cdata <- cdata %>%
      mutate(scenario = "historical")

  }

  if ("product" %in% colnames(cdata)){
    cdata <- cdata %>%
      rename(model = product)

  }

  all.ts <- bind_rows(all.ts,
                      cdata %>%
                        mutate(type = ctype))

}

saveRDS(all.ts,
        "./outputs/All.timeseries.CA.RDS")

###################################################################

All.trends.files <-
  c(list.files("/data/gent/vo/000/gvo00074/felicien/R/outputs/Drying.CB/",
               "*_trends.tif",
               full.names = TRUE),
    list.files("/data/gent/vo/000/gvo00074/felicien/R/outputs/Drying.CB/CMIP6/",
               "*_trends.tif",
               full.names = TRUE))

all.slopes <- data.frame()
for (ifile in seq(1,length(All.trends.files))){

  cfile <- All.trends.files[ifile]

  print(cfile)

  csplit <- strsplit(tools::file_path_sans_ext(
    basename(cfile)),"\\_")[[1]]

  if (length(csplit) == 3){
    cmodel <- csplit[1]
    cvar <- csplit[2]
    cperiod <- "historical"
    ctype <- "Observational"
  } else {
    cmodel <- csplit[1]
    cperiod <- csplit[2]
    cvar <- csplit[3]
    ctype <- "CMIP6"
  }

  cr <- rast(cfile)
  cintercept <- cr[[1]]
  cslope <- cr[[2]]

  cslope.df <- as.data.frame(cslope,xy = TRUE) %>%
    dplyr::rename(lon = x,
                  lat = y)

  cintercept.df <- as.data.frame(cintercept,xy = TRUE) %>%
    dplyr::rename(lon = x,
                  lat = y)


  all.slopes <- bind_rows(
    all.slopes,
    cslope.df %>%
      left_join(cintercept.df,
                by = c("lon","lat")) %>%
      mutate(model = cmodel,
             period = cperiod,
             var = cvar,
             type = ctype,
             rel_slope = slope_per_year/abs(intercept_t0)))

}

saveRDS(all.slopes,
        "./outputs/All.slopes.CA.RDS")

###################################################################


All.anomalies.files <-
  list.files("/data/gent/vo/000/gvo00074/felicien/R/outputs/Drying.CB/",
               "*_anomalies.tif",
               full.names = TRUE)

all.anomalies <- data.frame()
for (ifile in seq(1,length(All.anomalies.files))){

  cfile <- All.anomalies.files[ifile]

  print(cfile)

  csplit <- strsplit(tools::file_path_sans_ext(
    basename(cfile)),"\\_")[[1]]

  if (length(csplit) == 3){
    cmodel <- csplit[1]
    cvar <- csplit[2]
    cperiod <- "historical"
    ctype <- "Observational"
  } else {
    cmodel <- csplit[1]
    cperiod <- csplit[2]
    cvar <- csplit[3]
    ctype <- "CMIP6"
  }

  cr <- rast(cfile)
  pos <- which(year(time(cr)) == 2024)

  if (length(pos) == 0) next()

  canomalies <- mean(cr[[pos]],
                     na.rm = TRUE)

  canomalies.df <- as.data.frame(canomalies,xy = TRUE) %>%
    dplyr::rename(lon = x,
                  lat = y)

  all.anomalies <- bind_rows(
    all.anomalies,
    canomalies.df %>%
      mutate(model = cmodel,
             period = cperiod,
             var = cvar,
             type = ctype))

}

saveRDS(all.anomalies,
        "./outputs/All.anomalies.CA.RDS")


###################################################################

All.Zanomalies.files <-
  list.files("/data/gent/vo/000/gvo00074/felicien/R/outputs/Drying.CB/",
               "*_Zanomalies.tif",
               full.names = TRUE)

all.Zanomalies <- data.frame()
for (ifile in seq(1,length(All.Zanomalies.files))){

  cfile <- All.Zanomalies.files[ifile]

  print(cfile)

  csplit <- strsplit(tools::file_path_sans_ext(
    basename(cfile)),"\\_")[[1]]

  if (length(csplit) == 3){
    cmodel <- csplit[1]
    cvar <- csplit[2]
    cperiod <- "historical"
    ctype <- "Observational"
  } else {
    cmodel <- csplit[1]
    cperiod <- csplit[2]
    cvar <- csplit[3]
    ctype <- "CMIP6"
  }

  cr <- rast(cfile)
  pos <- which(year(time(cr)) == 2024)

  if (length(pos) == 0) next()

  canomalies <- mean(cr[[pos]],
                     na.rm = TRUE)

  canomalies.df <- as.data.frame(canomalies,xy = TRUE) %>%
    dplyr::rename(lon = x,
                  lat = y)

  all.Zanomalies <- bind_rows(
    all.Zanomalies,
    canomalies.df %>%
      mutate(model = cmodel,
             period = cperiod,
             var = cvar,
             type = ctype))

}

saveRDS(all.Zanomalies,
        "./outputs/All.Zanomalies.CA.RDS")


###################################################################

All.trendsanomalies.files <-
  c(list.files("/data/gent/vo/000/gvo00074/felicien/R/outputs/Drying.CB/",
               "*_trendsanomalies.tif",
               full.names = TRUE),
    list.files("/data/gent/vo/000/gvo00074/felicien/R/outputs/Drying.CB/CMIP6/",
               "*_trendsanomalies.tif",
               full.names = TRUE))

all.slopes.anomalies <- data.frame()
for (ifile in seq(1,length(All.trendsanomalies.files))){

  cfile <- All.trendsanomalies.files[ifile]

  print(cfile)

  csplit <- strsplit(tools::file_path_sans_ext(
    basename(cfile)),"\\_")[[1]]

  if (length(csplit) == 3){
    cmodel <- csplit[1]
    cvar <- csplit[2]
    cperiod <- "historical"
    ctype <- "Observational"
  } else {
    cmodel <- csplit[1]
    cperiod <- csplit[2]
    cvar <- csplit[3]
    ctype <- "CMIP6"
  }

  cr <- rast(cfile)
  cslope <- cr[[2]]

  cslope.df <- as.data.frame(cslope,xy = TRUE) %>%
    dplyr::rename(lon = x,
                  lat = y)

  all.slopes.anomalies <- bind_rows(
    all.slopes.anomalies,
    cslope.df %>%
      left_join(cintercept.df,
                by = c("lon","lat")) %>%
      mutate(model = cmodel,
             period = cperiod,
             var = cvar,
             type = ctype))

}

saveRDS(all.slopes.anomalies,
        "./outputs/All.slopesanomalies.CA.RDS")


###################################################################

All.trendsZanomalies.files <-
  c(list.files("/data/gent/vo/000/gvo00074/felicien/R/outputs/Drying.CB/",
               "*_trendsZanomalies.tif",
               full.names = TRUE),
    list.files("/data/gent/vo/000/gvo00074/felicien/R/outputs/Drying.CB/CMIP6/",
               "*_trendsZanomalies.tif",
               full.names = TRUE))

all.slopes.Zanomalies <- data.frame()
for (ifile in seq(1,length(All.trendsZanomalies.files))){

  cfile <- All.trendsZanomalies.files[ifile]

  print(cfile)

  csplit <- strsplit(tools::file_path_sans_ext(
    basename(cfile)),"\\_")[[1]]

  if (length(csplit) == 3){
    cmodel <- csplit[1]
    cvar <- csplit[2]
    cperiod <- "historical"
    ctype <- "Observational"
  } else {
    cmodel <- csplit[1]
    cperiod <- csplit[2]
    cvar <- csplit[3]
    ctype <- "CMIP6"
  }

  cr <- rast(cfile)
  cslope <- cr[[2]]

  cslope.df <- as.data.frame(cslope,xy = TRUE) %>%
    dplyr::rename(lon = x,
                  lat = y)

  all.slopes.Zanomalies <- bind_rows(
    all.slopes.Zanomalies,
    cslope.df %>%
      left_join(cintercept.df,
                by = c("lon","lat")) %>%
      mutate(model = cmodel,
             period = cperiod,
             var = cvar,
             type = ctype))

}

saveRDS(all.slopes.Zanomalies,
        "./outputs/All.slopesZanomalies.CA.RDS")



# scp /home/femeunier/Documents/projects/Drying.CB/scripts/Compile.all.outputs.R hpc:/kyukon/data/gent/vo/000/gvo00074/felicien/R/


rm(list = ls())

library(dplyr)
library(terra)
library(lubridate)

All.anomalies.files <-
  list.files("/data/gent/vo/000/gvo00074/felicien/R/outputs/Drying.CB/",
             "*_anomalies_Amazon.tif",
             full.names = TRUE)

all.anomalies <- data.frame()
for (ifile in seq(1,length(All.anomalies.files))){

  cfile <- All.anomalies.files[ifile]

  print(cfile)

  csplit <- strsplit(tools::file_path_sans_ext(
    basename(cfile)),"\\_")[[1]]

  if (length(csplit) == 4){
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
  pos <- which(year(time(cr)) == 2023)

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
        "./outputs/All.anomalies.Amazon.RDS")

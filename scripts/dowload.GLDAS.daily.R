rm(list = ls())

# Download daily GLDAS
library(ncdf4)
library(rvest)
library(reshape2)
library(dplyr)
library(Drying.CB)


years <- 2000:2025
all.days <- 1:366
dest.dir <- "/data/gent/vo/000/gvo00074/ED_common_data/met/GLDAS/daily/"

max.try = 100

# Function
retry.func <- function (expr, isError = function(x) inherits(x, "try-error"),
                        maxErrors = 5, sleep = 0) {
  attempts = 0
  retval = try(eval(expr))
  while (isError(retval)) {
    attempts = attempts + 1
    if (attempts >= maxErrors) {
      msg = sprintf("retry: too many retries [[%s]]",
                    utils::capture.output(utils::str(retval)))
      # warning(msg)
      return(0)
    }
    else {
      msg = sprintf("retry: error in attempt %i/%i [[%s]]",
                    attempts, maxErrors, utils::capture.output(utils::str(retval)))
      # warning(msg)
    }
    if (sleep > 0)
      Sys.sleep(sleep)
    retval = try(eval(expr))
  }
  return(1)
}

options(timeout=6000)

all.df.month <- data.frame()
for (cyear in years){

  op.file <- paste0(dest.dir,"GLDAS.year.tas.",cyear,".RDS")

  if (file.exists(op.file)){

    cdf.year <- readRDS(op.file)
    all.df.month <- bind_rows(all.df.month,
                              cdf.year)
    next()
  }

  all.df <- data.frame()

  for (cday in all.days){

    print(paste0(cyear," - ",cday))

    curl <- paste0("https://hydro1.gesdisc.eosdis.nasa.gov/data/GLDAS/GLDAS_NOAH025_3H.2.1/",
                  cyear,"/",sprintf("%03d",cday))


    # Read the HTML content of the directory listing
    page <- tryCatch(read_html(curl),
                     error = function(e) NULL)

    if (is.null(page)){
      closeAllConnections()
      next()
    }

    # Extract links that end with .nc4
    nc_files <- page %>%
      html_elements("a") %>%
      html_attr("href") %>%
      grep("\\.nc4$", ., value = TRUE)

    # Create full URLs
    full_urls <- unique(paste0(curl,"/",nc_files))

    if (length(full_urls) == 0){
      next()
    }

    all.attributes <- strsplit(basename(full_urls),
                               split = "\\.")

    dates <- as.character(as.vector(unlist(purrr::map_dfr(1:length(all.attributes),
                                                            function(i){
                                                              data.frame(var = all.attributes[[2]][2])}))))

    cyears <- as.numeric(substr(dates,2,5))
    cmonths <- as.numeric(substr(dates,6,7))
    cdays <- as.numeric(substr(dates,8,9))


    for (ifile in seq(1,length(full_urls))){

      cfile <- full_urls[ifile]

      cdest.file <- file.path(dest.dir,basename(cfile))
      does.file.exist <- file.exists(cdest.file)
      ctry <- 1

     while(!does.file.exist & ctry < max.try){

        dumb <- retry.func(system2("wget",
                                   c(cfile,
                                     "-P", dest.dir)),
                           maxErrors = 10,
                           sleep = 0)
        does.file.exist <- file.exists(cdest.file)
        ctry <- ctry + 1

        if (!does.file.exist){
          print("Pausing")
          Sys.sleep(10)
        }
      }

      if (!file.exists(cdest.file)){
        warning(paste0("Couldn't download:",cdest.file))
        next()
      }

      nc <- nc_open(cdest.file)

      all.lons <- ncvar_get(nc,"lon")
      all.lats <- ncvar_get(nc,"lat")

      pos.lats <- which(abs(all.lats) <= 30)
      lats <- all.lats[pos.lats]

      tmp <- ncvar_get(nc,"Tair_f_inst",
                       start = c(1,min(pos.lats),1),
                       count = c(-1,length(pos.lats),-1))

      Qair <- ncvar_get(nc,"Qair_f_inst",
                        start = c(1,min(pos.lats),1),
                        count = c(-1,length(pos.lats),-1))


      cdf <- melt(tmp) %>%
        mutate(lon = all.lons[Var1],
               lat = lats[Var2],
               year = cyears[ifile],
               month = cmonths[ifile],
               day = cdays[ifile]) %>%
        dplyr::select(lon,lat,year,month,day,value) %>%
        left_join(melt(Qair) %>%
                    mutate(lon = all.lons[Var1],
                           lat = lats[Var2],
                           year = cyears[ifile],
                           month = cmonths[ifile],
                           day = cdays[ifile]) %>%
                    rename(Qair = value),
                  by = c("lon","lat","year","month","day")) %>%
        mutate(vpd = vpd_from_T_q(value,Qair,101.325))

      nc_close(nc)

      all.df <- bind_rows(all.df,
                          cdf %>%
                            na.omit())

    }

    closeAllConnections()
  }

  system2("rm",
          c("-rf",
            paste0(dest.dir,"*",cyear,"*.nc4")))

  cdf.year <-  all.df %>%
    group_by(lon,lat,year,month,day) %>%
    summarise(tas.max = max(value,
                            na.rm = TRUE),
              tas.min = min(value,
                            na.rm = TRUE),
              tas.m = mean(value,
                           na.rm = TRUE),
              N = n(),
              .groups = "keep") %>%
    group_by(lon,lat,year,month) %>%
    summarise(tasmax = mean(tas.max,
                            na.rm = TRUE),
              tasmin = mean(tas.min,
                            na.rm = TRUE),
              N = sum(N),
              tas = mean(tas.m,
                         na.rm = TRUE),
              .groups = "keep")

  saveRDS(cdf.year,
          op.file)

  all.df.month <- bind_rows(all.df.month,
                            cdf.year)
}


# scp /home/femeunier/Documents/projects/Drying.CB/scripts/dowload.GLDAS.daily.R hpc:/kyukon/data/gent/vo/000/gvo00074/felicien/R

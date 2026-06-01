library(reticulate)
library(future)
library(purrr)
library(furrr)
library(lubridate)

setwd("./outputs/test.ERA5land/")

plan(multicore)

years <- 2024:2025
cmonths <- 1:12

c(years) %>%
  future_map(function(year) {

    # you need to have an account for downloaing the files
    # Read the documantion for how to setup your account and settings before trying this
    # https://confluence.ecmwf.int/display/CKB/How+to+download+ERA5#HowtodownloadERA5-3-DownloadERA5datathroughtheCDSAPI
    cdsapi <-import("cdsapi")
    c <- cdsapi$Client()

    c$retrieve(
      'reanalysis-era5-land',
      list(
        'product_type' = 'reanalysis',
        'format' = 'netcdf',
        'download_format' = 'unarchived',
        'day' = list('01','02','03',
                     '04','05','06',
                     '07','08','09',
                     '10','11','12',
                     '13','14','15',
                     '16','17','18',
                     '19','20','21',
                     '22','23','24',
                     '25','26','27',
                     '28','29','30',
                     '31'),
        'time' = list('00:00','03:00','06:00',
                      '09:00','12:00','15:00',
                      '18:00','21:00'),
        'month' = cmonths,
        'year' = as.character(year),
        'area' = "1/24/0/25", # max lat/min lon/min lat/max lon
        'variable' = list("2m_temperature")
      ),
      paste0('ERA5land',year,'.nc')
    )
  })

# scp /Users/felicien/Documents/projects/Drying.CB/scripts/download.ERA5land.R hpc:/data/gent/vo/000/gvo00074/felicien/R




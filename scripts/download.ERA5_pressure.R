library(reticulate)
library(future)
library(purrr)
library(furrr)
library(lubridate)

# setwd("/home/femeunier/Documents/projects/YGB/outputs/") # change this to your own working directory
setwd("/data/gent/vo/000/gvo00074/ED_common_data/met/CB/ERA5")

plan(multicore)

files.already.dow <- tools::file_path_sans_ext(list.files(
  getwd(),pattern = "*pressure*"))
yrs <- as.numeric(sub(".*\\_","",files.already.dow))
years <- 1940:2024
cmonths <- 1:12
years <- years[!(years %in% yrs)]

c(years) %>%
  future_map(function(year) {

    # you need to have an account for downloaing the files
    # Read the documantion for how to setup your account and settings before trying this
    # https://confluence.ecmwf.int/display/CKB/How+to+download+ERA5#HowtodownloadERA5-3-DownloadERA5datathroughtheCDSAPI
    cdsapi <-import("cdsapi")
    c <- cdsapi$Client()

    c$retrieve(
      'reanalysis-era5-pressure-levels',
      list(
        'product_type' = 'reanalysis',
        'format' = 'netcdf',
        'pressure_level' = list('500'),
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
        'area' = "25/-20/-25/60", #"5/-53/5.5/-52.5",
        'grid' = "0.5/0.5",
        'variable' = list("vertical_velocity")
      ),
      paste0('ERA5_CB_pressure_',year,'.nc')
    )
  })

# scp /Users/felicien/Documents/projects/Drying.CB/scripts/download.ERA5_pressure.R hpc:/data/gent/vo/000/gvo00074/felicien/R




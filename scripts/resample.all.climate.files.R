rm(list = ls())

library(raster)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(lubridate)
library(TrENDY.analyses)
library(tidyr)
library(raster)

# system2("rsync",
#         c("-avz",
#           "hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/*climate.RDS",
#           "./outputs/"))

files <- list.files("./outputs/",
                    pattern = "*.climate.RDS",
                    full.names = TRUE)
files <- files[grepl("3IMERG|Berk|CAMS|chirps|chirpsv3|CRU|CRUJRA3Q|ERA5|GLDAS|GPCC|MSWEP|NCEP",
                     files)]
files <- files[grepl("NCEP|GLDAS",
                     files)]


coord <- expand.grid(lon = seq(-179.75,179.75,0.5),
                     lat = seq(-30.25,30.25,0.5)) %>%
  mutate(value = 1)
craster <- raster::rasterFromXYZ(coord)

all.df <- data.frame()

overwrite = TRUE

for (ifile in seq(1,length(files))){

  cfile <- files[ifile]
  cdataset <- strsplit(basename(cfile),
                       "\\.")[[1]][2]
  print(cdataset)

  OP.file <- paste0(tools::file_path_sans_ext(cfile),".rspld.RDS")

  if (file.exists(OP.file) & !overwrite){
    next()
  }

  cdf <- readRDS(files[ifile])

  colnames(cdf)[colnames(cdf) %in% c("pr","MAP","Pmm",'precip',"prate")] <- "pre"
  colnames(cdf)[colnames(cdf) %in% c("tmn","tmin","T2MMIN")] <- "tasmin"
  colnames(cdf)[colnames(cdf) %in% c("tmx","tmax","T2MMAX")] <- "tasmax"
  colnames(cdf)[colnames(cdf) %in% c("tmp","T2MMEAN")] <- "tas"
  cn <- colnames(cdf)

  cdf <- cdf %>%
    dplyr::select(any_of(c("year","lon","lat","month","pre","tasmin","tasmax","tas")))

  cdf.rspld <- data.frame()

  for (cyear in unique(cdf$year)){

    print(paste("-",cyear))
    ccdf.rspld <- resample.df.all.col(bigdf = cdf %>%
                                        filter(year == cyear),
                                     raster2resample = craster,
                                     var.names = cn[which(cn %in% c('pre',"tas","tasmax","tasmin"))],
                                     res = NULL,
                                     verbose = FALSE)
    cdf.rspld <- bind_rows(cdf.rspld,
                           ccdf.rspld)
  }

  saveRDS(cdf.rspld,
         OP.file)

}

# scp /home/femeunier/Documents/projects/Drying.CB/scripts/resample.all.climate.files.R hpc:/data/gent/vo/000/gvo00074/felicien/R

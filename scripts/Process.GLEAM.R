rm(list = ls())

library(terra)

# To download:
# scp -P 2225 -r gleamuser@aether.ugent.be:/data/v4.3a/monthly/SMrz .
# password: GLEAM4#h-cel_111

all.vars <- list()
# c("SMs","SMrz","S","Ep")

for (cvar in c("SMs","SMrz","S","Ep")){

  all.vars[[cvar]] <- list()

  all.dates <- c()

  for (cyear in 1980:2025){

    print(paste0(cvar,"-",cyear))

    cfile <- paste0("/data/gent/vo/000/gvo00074/felicien/GLEAM/",
                    cvar,"_",cyear,"_GLEAM_v4.2a_MO.nc")
    cr <- rast(cfile)

    all.vars[[cvar]][[as.character(cyear)]] <- cr

    all.dates <- c(all.dates,
                   as.Date(paste0(cyear,"/",1:12,"/01")))
  }

  ctemp <-  crop(aggregate(rast(all.vars[[cvar]]),
                           fact = 5,
                           na.rm = TRUE),
                 ext(-180,180,-25,25))
  time(ctemp) <- all.dates

  all.vars[[cvar]] <- ctemp

  writeRaster(ctemp,
              paste0("/data/gent/vo/000/gvo00074/felicien/GLEAM/",
                     cvar,"_all_years.tif"),
              overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

}

# scp /Users/felicien/Documents/projects/Drying.CB/scripts/Process.GLEAM.R hpc:/data/gent/vo/000/gvo00074/felicien/R/

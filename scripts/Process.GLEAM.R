rm(list = ls())


library(terra)

all.vars <- list()
# c("SMs","SMrz","S","Ep")

for (cvar in c("SMs","SMrz","S","Ep")){

  all.vars[[cvar]] <- list()

  all.dates <- c()

  for (cyear in 1980:2024){

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

# scp /home/femeunier/Documents/projects/Drying.CB/scripts/Process.GLEAM.R hpc:/data/gent/vo/000/gvo00074/felicien/R/


A <- rast("~/Documents/data/GLEAM/Ep_all_years.tif")

Comp <- rast("~/Documents/data/GLEAM/Ep_2023_GLEAM_v4.2a_MO.nc")

plot(A[[526]])
plot(crop(Comp[[10]],
          ext(-180,180,-25,25)))

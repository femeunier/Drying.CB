rm(list = ls())

library(ncdf4)
library(terra)

# To download
# gdown --folder "https://drive.google.com/drive/folders/1Qr6-dQLhdbS-V459wPigPate8fuENAdB?usp=sharing" -O . --remaining-ok


files <- list.files("/data/gent/vo/000/gvo00074/felicien/MODIS",
                    pattern = "*.hdf$",full.names = TRUE)

vars <- c("\"CMG 0.05 Deg Monthly NDVI\"",
          "\"CMG 0.05 Deg Monthly EVI\"" ,
          "\"CMG 0.05 Deg Monthly pixel reliability\"")
nice.name <- c("NDVI",
                "EVI",
                "Pixel.Reliability")


for (ivar in seq(1,length(vars))){

  cvar <- vars[ivar]

  clist <- list()

  times <- as.Date(sub(".*A(\\d{7}).*", "\\1", basename(files)),
                     format = "%Y%j")

  for (ifile in seq(1,length(files))){

    print(paste0(cvar,"-",ifile))

    cfile <- files[ifile]

    A <- rast(cfile)
    A.crop <- crop(A,
                   ext(-180,180,-25,25))

    pos <- which(names(A.crop) == cvar)
    clist[[ifile]] <- A.crop[[pos]]
  }

  crast <- rast(clist)
  time(crast) <- times

  writeRaster(crast,
              paste0("./outputs/",nice.name[ivar],"_all.years.tif"),
              overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

}

# scp /home/femeunier/Documents/projects/Drying.CB/scripts/process.MODIS.data.R hpc:/data/gent/vo/000/gvo00074/felicien/R/



# A <- rast("~/Downloads/MOD13C2.A2000183.061.2020051165038.hdf")
# names(A)
# plot(A.crop[[2]])

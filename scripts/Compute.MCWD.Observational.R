àrm(list = ls())

library(raster)
library(terra)
library(lubridate)
library(SPEI)

mcwd_calc <- function(v, yy, uy, reset_jan = FALSE) {
  m    <- length(v)
  cur  <- 0                    # running cumulative deficit (≤ 0)
  out  <- rep(NA_real_, length(uy))  # one value per year (most negative)

  for (i in seq_len(m)) {
    # reset at Jan if using calendar-year flavour
    if (reset_jan && (i == 1L || yy[i] != yy[i-1L])) cur <- 0

    vi <- v[i]
    if (is.finite(vi)) {
      # update cumulative deficit; never allow > 0
      cur <- min(0, cur + vi)
      yi  <- match(yy[i], uy)
      # store the most negative value reached in that year
      out[yi] <- if (is.na(out[yi])) cur else min(out[yi], cur)
    }
  }
  out
}


get_dates <- function(x) {
  tt <- time(x)
  if (!is.null(tt)) return(as.Date(tt))
  # fallback: try to parse YYYY-MM or YYYY-MM-DD from layer names
  as.Date(sub(".*(\\d{4}-\\d{2}-\\d{2}|\\d{4}-\\d{2}).*$", "\\1", names(x)))
}

Ra_for_J <- function(lat_r, J) {
  dr    <- 1 + 0.033 * cos(2*pi * J/365)
  delta <- 0.409 * sin(2*pi * J/365 - 1.39)

  # x = -tan(phi) * tan(delta), clamp to [-1, 1] BEFORE acos
  x  <- -tan(lat_r) * tan(delta)        # SpatRaster
  x  <- clamp(x, -1, 1)
  ws <- acos(x)

  ((24*60)/pi) * Gsc * dr *
    (ws * sin(lat_r) * sin(delta) + cos(lat_r) * cos(delta) * sin(ws))
}

dir <- "/data/gent/vo/000/gvo00074/felicien/R/outputs/all.climate"
files <- list.files(dir,
                    pattern = "*.tif",
                    full.names = TRUE)

files.split <- strsplit(basename(files),"\\_")
models <- sapply(files.split,"[[",1)
models.uni <- sort(unique(models))
vars <- sapply(files.split,"[[",2)

vars.uni <- sort(unique(vars))

necessary.vars <- c("tasmin","tasmax","tas","pre")

Gsc <- 0.0820  # MJ m^-2 min^-1

for (cmodel in models.uni){

  cfiles <- files[which(models == cmodel)]
  cvars <- vars[which(models == cmodel)]

  missing.vars <- necessary.vars[!(necessary.vars %in% cvars)]

  if (length(missing.vars) > 0){
    for (cmissing.var in missing.vars){
      cfiles <- c(cfiles,
                  files[which(models == "MEM" & vars == cmissing.var)])
    }
  }

  tasmin <- rast(cfiles[grepl("(^|_)tasmin(_|\\.)", basename(cfiles))])
  tasmax <- rast(cfiles[grepl("(^|_)tasmax(_|\\.)", basename(cfiles))])
  tas    <- rast(cfiles[grepl("(^|_)tas(_|\\.)",    basename(cfiles))])  # excludes tasmin/max
  pre <- rast(cfiles[grepl("pre",cfiles)])

  R <- list(tasmin = tasmin, tasmax = tasmax, tas = tas, pre = pre)

  tz <- lapply(R, get_dates)
  common_dates <- as.Date(Reduce(intersect, tz))

  R_common <- Map(function(x, tvec) {
    idx <- match(common_dates, tvec)
    x[[idx]]
  }, R, tz)

  tmin <- R_common$tasmin
  tmax  <- R_common$tasmax
  tmean <- R_common$tas
  pre     <- R_common$pre

  tt      <- as.Date(time(tmean))
  ndays   <- as.vector(days_in_month(tt))
  middate <- floor_date(tt, "month") + days(round((ndays - 1) / 2))
  Jmid    <- yday(middate)

  lat_deg <- init(tmin[[1]], "y")           # one-layer raster of latitude (°)

  # 5) Helper to compute Ra vector (MJ m^-2 day^-1) for a given latitude
  Ra_vec <- function(lat_deg_one, J) {
    if (is.na(lat_deg_one)) return(rep(NA_real_, length(J)))
    lat <- lat_deg_one * pi/180
    dr    <- 1 + 0.033 * cos(2*pi * J/365)
    delta <- 0.409 * sin(2*pi * J/365 - 1.39)
    x     <- -tan(lat) * tan(delta)
    x     <- pmin(1, pmax(-1, x))
    ws    <- acos(x)
    Gsc   <- 0.0820  # MJ m^-2 min^-1
    ((24*60)/pi) * Gsc * dr * (ws * sin(lat) * sin(delta) + cos(lat) * cos(delta) * sin(ws))
  }

  X <- c(tmin, tmax, lat_deg)

  start_ym <- c(lubridate::year(tt[1]), lubridate::month(tt[1]))

  pet_mm_month <- app(X, fun = function(v) {
    mm <- (length(v) - 1L) / 2L
    if (mm <= 0 || mm != floor(mm)) return(NA_real_)  # layer count mismatch
    mm <- as.integer(mm)

    tn  <- v[1:mm]
    tx  <- v[(mm+1):(2*mm)]
    lat <- v[length(v)]

    # If latitude is missing, no way to compute Ra → all NA for this cell
    if (!is.finite(lat)) return(rep(NA_real_, mm))

    # Fix occasional inversions: swap where tx < tn
    swap <- is.finite(tn) & is.finite(tx) & (tx < tn)
    if (any(swap)) {
      tmp    <- tn[swap]
      tn[swap] <- tx[swap]
      tx[swap] <- tmp
    }

    # Build Ra for the full series; must be length mm and finite
    Ra <- Ra_vec(lat, Jmid[seq_len(mm)])
    if (length(Ra) != mm) return(rep(NA_real_, mm))
    Ra[!is.finite(Ra)] <- NA_real_

    # Use ts so months/years are preserved; incomplete years are fine
    tn_ts <- ts(tn, frequency = 12, start = start_ym)
    tx_ts <- ts(tx, frequency = 12, start = start_ym)

    pet_mon <- suppressWarnings(
      SPEI::hargreaves(Tmin = tn_ts, Tmax = tx_ts, Ra = Ra, na.rm = TRUE)
    )

    pet <- as.numeric(pet_mon)
    pet[!is.finite(pet)] <- NA_real_
    pmax(0, pet)
  })

  time(pet_mm_month) <- time(tmax)

  writeRaster(pet_mm_month,
              paste0("/data/gent/vo/000/gvo00074/felicien/R/outputs/all.climate/",
                     cmodel,"_PET_all.years.tif"),
              overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))


  wb <- pre - pet_mm_month   # SpatRaster, mm/month

  tt <- as.Date(time(wb))
  yy <- year(tt)
  uy <- sort(unique(yy))
  ny <- length(uy)

  mcwd_year_calendar   <- app(wb, fun = function(v)
    mcwd_calc(v, yy, uy, reset_jan = TRUE))

  time(mcwd_year_calendar)    <- as.Date(paste0(uy, "-07-01"))

  writeRaster(mcwd_year_calendar,
              paste0("/data/gent/vo/000/gvo00074/felicien/R/outputs/all.climate/",
                     cmodel,"_MCWD_all.years.tif"),
              overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

}

# scp /home/femeunier/Documents/projects/Drying.CB/scripts/Compute.MCWD.Observational.R hpc:/data/gent/vo/000/gvo00074/felicien/R/


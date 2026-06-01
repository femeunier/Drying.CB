rm(list = ls())

library(terra)
library(lubridate)

###############################################################################
# Functions
###############################################################################

get_dates <- function(x) {
  tt <- time(x)

  if (!is.null(tt)) {
    tt <- as.Date(tt, origin = "1970-01-01")
    if (all(!is.na(tt))) return(tt)
  }

  nms <- names(x)

  # YYYY_M
  out1 <- suppressWarnings(as.Date(paste0(
    sub("^(\\d{4})_(\\d{1,2})$", "\\1", nms), "-",
    sprintf("%02d", as.integer(sub("^(\\d{4})_(\\d{1,2})$", "\\2", nms))), "-01"
  )))
  if (all(!is.na(out1))) return(out1)

  # YYYY-MM-DD somewhere in the name
  out2 <- suppressWarnings(as.Date(sub(".*(\\d{4}-\\d{2}-\\d{2}).*$", "\\1", nms)))
  if (all(!is.na(out2))) return(out2)

  # YYYY-MM somewhere in the name
  out3 <- suppressWarnings(as.Date(paste0(
    sub(".*(\\d{4}-\\d{2}).*$", "\\1", nms), "-01"
  )))
  if (all(!is.na(out3))) return(out3)

  stop("Could not parse dates from raster time or layer names.")
}

mcwd_one_window <- function(v) {
  cur <- 0
  out <- NA_real_

  for (i in seq_along(v)) {
    vi <- v[i]
    if (is.finite(vi)) {
      cur <- min(0, cur + vi)
      out <- if (is.na(out)) cur else min(out, cur)
    }
  }

  out
}

mcwd_rolling_12m <- function(v, k = 12) {
  n <- length(v)
  out <- rep(NA_real_, n)

  if (n < k) return(out)

  for (i in k:n) {
    win <- v[(i - k + 1):i]
    out[i] <- mcwd_one_window(win)
  }

  out
}

###############################################################################
# Paths
###############################################################################

precip_dir <- "/data/gent/vo/000/gvo00074/felicien/R/outputs/all.climate"
pet_file   <- "/data/gent/vo/000/gvo00074/felicien/GLEAM/Ep_all_years.tif"
out_dir    <- "/data/gent/vo/000/gvo00074/felicien/R/outputs/all.climate"

###############################################################################
# PET
###############################################################################

PET <- rast(pet_file)
pet_dates <- get_dates(PET)
time(PET) <- pet_dates
print(summary(pet_dates))

###############################################################################
# Precip files
###############################################################################

files <- list.files(
  precip_dir,
  pattern = "_pre_.*\\.tif$",
  full.names = TRUE
)

# product names = part before "_pre_"
prod_names <- sub("_pre_.*$", "", basename(files))

precip_df <- data.frame(
  product = prod_names,
  file = files,
  stringsAsFactors = FALSE
)


############
# SPEI param

spei_scale <- 6   # e.g. 1, 3, 6, 12 months

###############################################################################
# Loop over precipitation products
###############################################################################

for (i in seq_len(nrow(precip_df))) {

  cprod <- precip_df$product[i]
  cfile <- precip_df$file[i]

  message("Processing precipitation product: ", cprod)

  pre <- rast(cfile)
  pre_dates <- get_dates(pre)
  print(summary(pre_dates))
  time(pre) <- pre_dates

  PET_sub_full <- PET
  pet_dates <- get_dates(PET_sub_full)
  time(PET_sub_full) <- pet_dates

  common_dates <- intersect(pre_dates, pet_dates)
  common_dates <- sort(as.Date(common_dates, origin = "1970-01-01"))

  if (length(common_dates) < 12) {
    message("Skipping ", cprod, ": fewer than 12 common months.")
    next
  }

  idx_pre <- match(common_dates, pre_dates)
  idx_pet <- match(common_dates, pet_dates)

  pre_sub <- pre[[idx_pre]]
  PET_sub <- PET_sub_full[[idx_pet]]

  time(pre_sub) <- common_dates
  time(PET_sub) <- common_dates

  time(pre_sub) <- common_dates
  time(PET_sub) <- common_dates

  # harmonize CRS label if needed
  crs(pre_sub) <- crs(PET_sub)

  # crop to PET extent
  pre_crop <- crop(pre_sub, ext(PET_sub))

  # aggregate precipitation to PET resolution if needed
  res_pre <- res(pre_crop)
  res_pet <- res(PET_sub)

  if (all(abs(res_pre - res_pet) < 1e-10)) {
    pre_on_pet <- resample(pre_crop, PET_sub, method = "near")
  } else {
    fact_x <- res_pet[1] / res_pre[1]
    fact_y <- res_pet[2] / res_pre[2]

    if (abs(fact_x - round(fact_x)) < 1e-10 &&
        abs(fact_y - round(fact_y)) < 1e-10 &&
        round(fact_x) == round(fact_y)) {

      pre_agg <- aggregate(
        pre_crop,
        fact = round(fact_x),
        fun = mean,
        na.rm = TRUE
      )

      pre_on_pet <- resample(pre_agg, PET_sub, method = "near")
    } else {
      pre_on_pet <- resample(pre_crop, PET_sub, method = "bilinear")
    }
  }

  compareGeom(pre_on_pet, PET_sub, stopOnError = TRUE)

  # water balance
  wb <- pre_on_pet - PET_sub
  time(wb) <- common_dates

  # water balance
  wb0 <- pre_on_pet - 100*(PET_sub**0)
  time(wb0) <- common_dates


  writeRaster(
    wb,
    file.path(
      out_dir,
      paste0(cprod, "_wb.tif")
    ),
    overwrite = TRUE,
    gdal = c("COMPRESS=NONE", "TFW=YES")
  )

  writeRaster(
    wb0,
    file.path(
      out_dir,
      paste0(cprod, "_wb0.tif")
    ),
    overwrite = TRUE,
    gdal = c("COMPRESS=NONE", "TFW=YES")
  )


  # rolling 12-month MCWD
  mcwd_monthly_12m <- app(wb, fun = function(v) {
    mcwd_rolling_12m(v, k = 12)
  })

  # rolling 12-month MCWD with 100mm PET
  mcwd_monthly_12m0 <- app(wb0, fun = function(v) {
    mcwd_rolling_12m(v, k = 12)
  })

  time(mcwd_monthly_12m) <- common_dates
  time(mcwd_monthly_12m0) <- common_dates

  # write output
  out_file <- file.path(
    out_dir,
    paste0(cprod, "_MCWD_rolling12m_GLEAMPET_all.months.tif")
  )

  # write output
  out_file0 <- file.path(
    out_dir,
    paste0(cprod, "_MCWD_rolling12m_PET100mm_all.months.tif")
  )

  writeRaster(
    mcwd_monthly_12m,
    out_file,
    overwrite = TRUE,
    gdal = c("COMPRESS=NONE", "TFW=YES")
  )

  writeRaster(
    mcwd_monthly_12m0,
    out_file0,
    overwrite = TRUE,
    gdal = c("COMPRESS=NONE", "TFW=YES")
  )

  ########################################################################
  # ---- SPEI ----

  spei_rast <- app(wb, fun = function(v) {

    if (sum(is.finite(v)) < 24) {
      return(rep(NA_real_, length(v)))
    }

    wb_ts <- ts(v,
                frequency = 12,
                start = c(year(common_dates[1]), month(common_dates[1])))

    spei_fit <- tryCatch(
      SPEI::spei(wb_ts,
                 scale = spei_scale,
                 na.rm = TRUE,
                 verbose = FALSE),
      error = function(e) NULL
    )

    if (is.null(spei_fit)) {
      return(rep(NA_real_, length(v)))
    }

    out <- as.numeric(spei_fit$fitted)
    out[!is.finite(out)] <- NA_real_
    out
  })

  time(spei_rast) <- common_dates

  writeRaster(spei_rast,
              paste0("/data/gent/vo/000/gvo00074/felicien/R/outputs/all.climate/",
                     cprod, "_SPEI", spei_scale, "_all.years.tif"),
              overwrite = TRUE,
              gdal = c("COMPRESS=NONE", "TFW=YES"))


  spei_rast0 <- app(wb0, fun = function(v) {

    if (sum(is.finite(v)) < 24) {
      return(rep(NA_real_, length(v)))
    }

    wb_ts <- ts(v,
                frequency = 12,
                start = c(year(common_dates[1]), month(common_dates[1])))

    spei_fit <- tryCatch(
      SPEI::spei(wb_ts,
                 scale = spei_scale,
                 na.rm = TRUE,
                 verbose = FALSE),
      error = function(e) NULL
    )

    if (is.null(spei_fit)) {
      return(rep(NA_real_, length(v)))
    }

    out <- as.numeric(spei_fit$fitted)
    out[!is.finite(out)] <- NA_real_
    out
  })

  time(spei_rast0) <- common_dates

  writeRaster(spei_rast,
              paste0("/data/gent/vo/000/gvo00074/felicien/R/outputs/all.climate/",
                     cprod, "_SPEI0", spei_scale, "_all.years.tif"),
              overwrite = TRUE,
              gdal = c("COMPRESS=NONE", "TFW=YES"))

}

# scp /Users/felicien/Documents/projects/Drying.CB/scripts/Compute.monthly.MCWD.Observational.R hpc:/data/gent/vo/000/gvo00074/felicien/R/
# scp /Users/felicien/Documents/projects/Drying.CB/scripts/Compute.monthly.MCWD.Observational.R hydra:/user/gent/425/vsc42558/R/

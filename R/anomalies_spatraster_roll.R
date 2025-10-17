anomalies_spatraster_roll <- function(
    input,
    baseline_start = as.Date("1981-01-01"),
    baseline_end   = as.Date("2010-12-31"),
    detrend = FALSE,
    roll_window  = 12L,
    roll_align   = c("right","center","left"),
    roll_min_obs = NULL,
    float32 = TRUE,
    show_progress = FALSE,
    tmpdir = NULL,
    todisk = TRUE
) {
  stopifnot(inherits(input, "SpatRaster"))
  roll_align <- match.arg(roll_align)
  if (is.null(roll_min_obs)) roll_min_obs <- ceiling(roll_window/2)

  # ---- terra options (socket-safe) ----
  old <- terra::terraOptions()
  on.exit(terra::terraOptions(progress = old$progress,
                              tempdir  = old$tempdir,
                              todisk   = old$todisk,
                              datatype = old$datatype,
                              memfrac  = old$memfrac), add = TRUE)
  if (!show_progress) terra::terraOptions(progress = 0)
  if (!is.null(tmpdir)) terra::terraOptions(tempdir = tmpdir)
  if (isTRUE(todisk))  terra::terraOptions(todisk  = TRUE)
  if (float32)         terra::terraOptions(datatype = "FLT4S")
  terra::terraOptions(memfrac = 0.8)

  # ---- helpers ----
  get_year  <- function(d) as.POSIXlt(d)$year + 1900L
  get_month <- function(d) as.POSIXlt(d)$mon  + 1L

  # Streaming OLS across layers (intercept at x=0, slope per x-unit)
  ols_trend_stream <- function(Y, x_centered) {
    stopifnot(nlyr(Y) == length(x_centered))
    proto <- Y[[1]]; zero <- proto * 0
    sumy <- sumx <- sumxx <- sumxy <- nvalid <- zero
    for (i in seq_len(nlyr(Y))) {
      yi  <- Y[[i]]
      vi  <- !is.na(yi)
      yi0 <- ifel(vi, yi, 0)
      xi  <- x_centered[i]; xi2 <- xi * xi
      nvalid <- nvalid + ifel(vi, 1, 0)
      sumy   <- sumy   + yi0
      sumx   <- sumx   + ifel(vi, xi, 0)
      sumxx  <- sumxx  + ifel(vi, xi2, 0)
      sumxy  <- sumxy  + (yi0 * xi)
    }
    den   <- nvalid * sumxx - sumx * sumx
    slope <- (nvalid * sumxy - sumx * sumy) / den
    slope <- ifel(den == 0, NA, slope)
    intercept <- (sumy - slope * sumx) / nvalid
    intercept <- ifel(nvalid == 0, NA, intercept)
    names(intercept) <- "intercept_t0"; names(slope) <- "slope_per_year"
    c(intercept, slope)
  }

  # Rolling mean via terra::roll(); enforce min_obs with a rolled count
  roll_mean_terra <- function(x, k, align, min_obs) {
    n <- nlyr(x); if (k < 1L || k > n) return(rast())
    tt <- time(x); if (!inherits(tt, "Date")) tt <- as.Date(tt); time(x) <- tt

    type <- switch(align,
                   right  = "from",   # window ends at current layer
                   left   = "to",     # window starts at current layer
                   center = "around") # centered window

    out <- roll(x, n = k, fun = mean, type = type, na.rm = TRUE)
    cnt <- roll(ifel(is.na(x), 0, 1), n = k, fun = sum,  type = type, na.rm = FALSE)
    out <- ifel(cnt >= min_obs, out, NA)

    names(out) <- paste0("rollmean_", format(tt, "%Y-%m"))
    time(out)  <- tt
    out
  }

  # ---- time ----
  tt <- time(input)
  if (is.null(tt)) stop("`input` must have a time vector (terra::time(input)).")
  if (!inherits(tt, "Date")) tt <- as.Date(tt)
  time(input) <- tt

  # ---- baseline subset ----
  i_base <- which(!is.na(tt) & tt >= baseline_start & tt <= baseline_end)
  if (!length(i_base)) stop("No layers fall within the chosen baseline period.")
  rb <- input[[i_base]]; tb <- tt[i_base]

  # continuous year coordinate (month midpoint), centered
  tn_full <- get_year(tt) + (get_month(tt) - 0.5)/12
  tn_base <- get_year(tb) + (get_month(tb) - 0.5)/12
  t0_full <- mean(tn_full, na.rm = TRUE)
  t0_base <- mean(tn_base, na.rm = TRUE)
  tn_c_full <- tn_full - t0_full
  tn_c_base <- tn_base - t0_base

  # ---- baseline trend on RAW (streaming OLS) ----
  trend_raw <- ols_trend_stream(input[[i_base]], tn_c_base)
  names(trend_raw) <- c("intercept_t0","slope_per_year")

  # optional detrend for anomaly computation
  if (isTRUE(detrend)) {
    input <- input - (trend_raw[[1]] + trend_raw[[2]] * (tn_full - t0_base))
  }

  # ---- monthly climatology & SD on baseline (compiled sums) ----
  m_all  <- get_month(tt)
  m_base <- get_month(tb)

  ones_b <- ifel(is.na(rb), 0, 1)
  cnt_m  <- tapp(ones_b,      index = m_base, fun = sum, na.rm = TRUE)
  sum_m  <- tapp(rb,          index = m_base, fun = sum, na.rm = TRUE)
  ssq_m  <- tapp(rb * rb,     index = m_base, fun = sum, na.rm = TRUE)

  mean_m <- sum_m / cnt_m
  var_m  <- (ssq_m - (sum_m * sum_m) / cnt_m) / (cnt_m - 1)
  var_m  <- ifel(cnt_m <= 1, NA, var_m)
  sd_m   <- sqrt(ifel(var_m < 0, 0, var_m))

  present <- sort(unique(m_base))
  na_r <- rb[[1]] * NA_real_
  make_full12 <- function(stack_idx) {
    lst <- lapply(1:12, function(m) { pos <- match(m, present); if (is.na(pos)) na_r else stack_idx[[pos]] })
    out <- rast(lst); names(out) <- month.abb; out
  }
  clim12 <- make_full12(mean_m)
  sd12   <- make_full12(sd_m)

  # ---- anomalies & z ----
  clim_exp <- clim12[[m_all]]
  sd_exp   <- sd12[[m_all]]

  anom  <- input - clim_exp; names(anom)  <- paste0("anom_", names(input))
  z_anom <- (input - clim_exp) / sd_exp
  z_anom <- ifel(sd_exp == 0, NA, z_anom); names(z_anom) <- paste0("z_", names(input))

  # ---- trends (streaming OLS) ----
  trend_anom   <- ols_trend_stream(anom,   tn_c_full); names(trend_anom)   <- c("anom_intercept_t0","anom_slope_per_year")
  trend_z_anom <- ols_trend_stream(z_anom, tn_c_full); names(trend_z_anom) <- c("z_anom_intercept_t0","z_anom_slope_per_year")

  # ---- rolling means via terra::roll() and their trends ----
  r_input <- roll_mean_terra(input,  k = roll_window, align = roll_align, min_obs = roll_min_obs)
  r_anom  <- roll_mean_terra(anom,   k = roll_window, align = roll_align, min_obs = roll_min_obs)
  r_zanom <- roll_mean_terra(z_anom, k = roll_window, align = roll_align, min_obs = roll_min_obs)

  if (nlyr(r_input) > 0) {
    tt_r  <- time(r_input)
    tnr   <- get_year(tt_r) + (get_month(tt_r) - 0.5)/12
    t0r   <- mean(tnr, na.rm = TRUE)
    tnr_c <- tnr - t0r
    trend_r_input <- ols_trend_stream(r_input,  tnr_c); names(trend_r_input) <- c("roll_input_intercept_t0","roll_input_slope_per_year")
    trend_r_anom  <- ols_trend_stream(r_anom,   tnr_c); names(trend_r_anom)  <- c("roll_anom_intercept_t0","roll_anom_slope_per_year")
    trend_r_zanom <- ols_trend_stream(r_zanom,  tnr_c); names(trend_r_zanom) <- c("roll_z_anom_intercept_t0","roll_z_anom_slope_per_year")
  } else {
    trend_r_input <- trend_r_anom <- trend_r_zanom <- NULL
  }

  # warn if baseline misses months
  missing_months <- month.abb[!(1:12 %in% present)]
  if (length(missing_months)) warning("Baseline missing months: ", paste(missing_months, collapse = ", "),
                                      ". Corresponding climatology/SD months are NA.")

  list(
    anom               = anom,
    z_anom             = z_anom,
    clim12             = clim12,
    sd12               = sd12,
    trend              = trend_raw,
    trend_anom         = trend_anom,
    trend_z_anom       = trend_z_anom,
    roll_mean_input    = r_input,
    roll_mean_anom     = r_anom,
    roll_mean_z_anom   = r_zanom,
    trend_roll_input   = trend_r_input,
    trend_roll_anom    = trend_r_anom,
    trend_roll_z_anom  = trend_r_zanom,
    trend_t0           = as.Date(round(mean(as.numeric(time(input)), na.rm = TRUE)), origin = "1970-01-01"),
    roll_times         = if (nlyr(r_input) > 0) time(r_input) else as.Date(character())
  )
}

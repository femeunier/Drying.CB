anomalies_spatraster <- function(input,
                                 baseline_start = as.Date("1900-01-01"),
                                 baseline_end   = as.Date("1900-12-31"),
                                 detrend = FALSE) {

  # Time vector
  tt <- time(input)
  if (is.null(tt)) stop("`input` must have a time vector (terra::time(input)).")
  if (!inherits(tt, "Date")) tt <- as.Date(tt)

  # Keep a copy of the raw series for trend estimation
  input_raw <- input

  # Baseline subset
  i_base <- which(!is.na(tt) & tt >= baseline_start & tt <= baseline_end)
  if (length(i_base) == 0) stop("No layers fall within the chosen baseline period.")
  rb <- input[[i_base]]
  tb <- tt[i_base]

  # numeric time in years for regressions
  tn_years <- year(tb) + (month(tb) - 1/2)/12
  t0_years <- mean(tn_years, na.rm = TRUE)


  # --- linear trend on RAW data (two-layer raster: intercept_at_t0, slope_per_year)
  trend_stack <- app(input_raw[[i_base]],
                     fun = trend_fit,
                     tn = tn_years, t0 = t0_years)


  names(trend_stack) <- c("intercept_t0", "slope_per_year")

  # --- optional detrend for anomaly computation
  if (isTRUE(detrend)) {
    input <- app(input, fun = detrend_linear, tn = tn_years)
  }



  # Month indices without format()
  m_all  <- as.POSIXlt(tt)$mon + 1L
  m_base <- as.POSIXlt(tb)$mon + 1L

  # Monthly climatology & SD via tapp
  clim_idx <- tapp(rb, index = m_base, fun = function(x) mean(x, na.rm = TRUE))
  sd_idx   <- tapp(rb, index = m_base, fun = function(x) stats::sd(x, na.rm = TRUE))

  # Rebuild full Jan..Dec stacks (fill missing months with NA rasters)
  present <- sort(unique(m_base))
  make_full12 <- function(stack_idx) {
    lst <- lapply(1:12, function(m) {
      pos <- match(m, present)
      if (is.na(pos)) app(rb[[1]], fun = function(x) NA_real_) else stack_idx[[pos]]
    })
    out <- rast(lst); names(out) <- month.abb
    out
  }
  clim12 <- make_full12(clim_idx)
  sd12   <- make_full12(sd_idx)

  # Expand to full timeline and compute anomalies
  clim_expanded <- clim12[[m_all]]
  sd_expanded   <- sd12[[m_all]]

  anom <- input - clim_expanded
  names(anom) <- paste0("anom_", names(input))

  z_anom <- (input - clim_expanded) / sd_expanded
  z_anom <- ifel(sd_expanded == 0, NA, z_anom)
  names(z_anom) <- paste0("z_", names(input))

  # warn if baseline misses months
  missing_months <- month.abb[!(1:12 %in% present)]
  if (length(missing_months)) {
    warning("Baseline missing months: ", paste(missing_months, collapse = ", "),
            ". Corresponding climatology/SD months are NA.")
  }

  list(
    anom        = anom,
    z_anom      = z_anom,
    clim12      = clim12,
    sd12        = sd12,
    trend       = trend_stack,   # two layers: intercept_t0, slope_per_year
    trend_t0    = as.Date( round(mean(as.numeric(tt), na.rm=TRUE)) , origin = "1970-01-01")
  )
}

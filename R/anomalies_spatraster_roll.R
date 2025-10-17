anomalies_spatraster <- function(input,
                                 baseline_start = as.Date("1981-01-01"),
                                 baseline_end   = as.Date("2010-12-31"),
                                 detrend = FALSE,
                                 roll_window  = 12L,
                                 roll_align   = c("right","center","left"),
                                 roll_min_obs = NULL) {

  roll_align <- match.arg(roll_align)
  if (is.null(roll_min_obs)) roll_min_obs <- ceiling(roll_window/2)

  # --- helpers ---
  get_year  <- function(d) as.POSIXlt(d)$year + 1900L
  get_month <- function(d) as.POSIXlt(d)$mon  + 1L

  # map our align arg to terra::roll "type"
  roll_type <- switch(roll_align,
                      right  = "from",    # window ends at current layer
                      left   = "to",      # window starts at current layer
                      center = "around")  # centered

  if (!exists("trend_fit", mode = "function")) {
    trend_fit <- function(v, tn, t0) {
      if (all(!is.finite(v))) return(c(NA_real_, NA_real_))
      ok <- is.finite(v) & is.finite(tn)
      if (sum(ok) < 2) return(c(NA_real_, NA_real_))
      X <- cbind(1, tn[ok] - t0)
      beta <- tryCatch(qr.solve(X, v[ok]),
                       error = function(e) c(NA_real_, NA_real_))
      if (any(!is.finite(beta))) return(c(NA_real_, NA_real_))
      beta  # c(intercept_at_t0, slope_per_year)
    }
  }

  # --- time vector ---
  tt <- time(input)
  if (is.null(tt)) stop("`input` must have a time vector (terra::time(input)).")
  if (!inherits(tt, "Date")) tt <- as.Date(tt)
  time(input) <- tt  # ensure Dates are written back (avoids format() issues)

  nlay <- nlyr(input)
  input_raw <- input

  # --- baseline ---
  i_base <- which(!is.na(tt) & tt >= baseline_start & tt <= baseline_end)
  if (length(i_base) == 0) stop("No layers fall within the chosen baseline period.")
  rb <- input[[i_base]]
  tb <- tt[i_base]

  tn        <- get_year(tt) + (get_month(tt) - 0.5)/12
  tn_years  <- get_year(tb) + (get_month(tb) - 0.5)/12
  t0_years  <- mean(tn_years, na.rm = TRUE)

  # raw trend over baseline
  trend_stack <- app(input_raw[[i_base]], fun = trend_fit, tn = tn_years, t0 = t0_years)
  names(trend_stack) <- c("intercept_t0", "slope_per_year")

  # optional detrend for anomaly calc
  if (isTRUE(detrend)) {
    input <- input - (trend_stack[[1]] + trend_stack[[2]] * tn)
  }

  # --- climatology & SD ---
  m_all  <- get_month(tt)
  m_base <- get_month(tb)

  clim_idx <- tapp(rb, index = m_base, fun = function(x) mean(x, na.rm = TRUE))
  sd_idx   <- tapp(rb, index = m_base, fun = function(x) stats::sd(x, na.rm = TRUE))

  present <- sort(unique(m_base))
  make_full12 <- function(stack_idx) {
    lst <- lapply(1:12, function(m) {
      pos <- match(m, present)
      if (is.na(pos)) app(rb[[1]], fun = function(x) NA_real_) else stack_idx[[pos]]
    })
    out <- rast(lst); names(out) <- month.abb; out
  }
  clim12 <- make_full12(clim_idx)
  sd12   <- make_full12(sd_idx)

  # --- anomalies & z-anomalies ---
  clim_expanded <- clim12[[m_all]]
  sd_expanded   <- sd12[[m_all]]

  anom <- input - clim_expanded
  names(anom) <- paste0("anom_", names(input))

  z_anom <- (input - clim_expanded) / sd_expanded
  z_anom <- ifel(sd_expanded == 0, NA, z_anom)
  names(z_anom) <- paste0("z_", names(input))

  # --- trends over FULL series ---
  trend_anom <- app(anom,  fun = trend_fit, tn = tn, t0 = mean(tn, na.rm = TRUE))
  names(trend_anom) <- c("anom_intercept_t0", "anom_slope_per_year")

  trend_z_anom <- app(z_anom, fun = trend_fit, tn = tn, t0 = mean(tn, na.rm = TRUE))
  names(trend_z_anom) <- c("z_anom_intercept_t0", "z_anom_slope_per_year")

  # --- rolling means via terra::roll() & their trends ---------------------
  # roll() keeps the same number of layers as input; we enforce min_obs by
  # rolling a count (sum of non-NA indicators) and masking when < roll_min_obs.

  # helper to roll a mean with min-obs mask
  roll_mean <- function(x) {
    # count of valid obs in the window at each layer
    cnt <- roll(ifel(is.na(x), 0, 1), n = roll_window, fun = sum, type = roll_type, na.rm = FALSE)
    # mean with na.rm=TRUE, then mask by min obs
    rm  <- roll(x, n = roll_window, fun = mean, type = roll_type, na.rm = TRUE)
    rm  <- ifel(cnt >= roll_min_obs, rm, NA)
    # tidy names & time
    names(rm) <- paste0("rollmean_", format(tt, "%Y-%m"))
    time(rm)  <- tt
    rm
  }

  roll_mean_input  <- roll_mean(input)
  roll_mean_anom   <- roll_mean(anom)
  roll_mean_z_anom <- roll_mean(z_anom)

  # time for rolled series (same length as tt)
  tt_roll <- tt
  tn_roll <- get_year(tt_roll) + (get_month(tt_roll) - 0.5)/12

  # trends on rolling means (reusing your trend_fit)
  trend_roll_input <- app(roll_mean_input, fun = trend_fit,
                          tn = tn_roll, t0 = mean(tn_roll, na.rm = TRUE))
  names(trend_roll_input) <- c("roll_input_intercept_t0", "roll_input_slope_per_year")

  trend_roll_anom <- app(roll_mean_anom,  fun = trend_fit,
                         tn = tn_roll, t0 = mean(tn_roll, na.rm = TRUE))
  names(trend_roll_anom) <- c("roll_anom_intercept_t0", "roll_anom_slope_per_year")

  trend_roll_z_anom <- app(roll_mean_z_anom, fun = trend_fit,
                           tn = tn_roll, t0 = mean(tn_roll, na.rm = TRUE))
  names(trend_roll_z_anom) <- c("roll_z_anom_intercept_t0", "roll_z_anom_slope_per_year")

  # warn if baseline misses months
  missing_months <- month.abb[!(1:12 %in% present)]
  if (length(missing_months)) {
    warning("Baseline missing months: ", paste(missing_months, collapse = ", "),
            ". Corresponding climatology/SD months are NA.")
  }

  list(
    anom               = anom,
    z_anom             = z_anom,

    clim12             = clim12,
    sd12               = sd12,

    trend              = trend_stack,      # raw trend over baseline
    trend_anom         = trend_anom,       # trend on anomalies (full series)
    trend_z_anom       = trend_z_anom,     # trend on z-anomalies (full series)

    trend_roll_input   = trend_roll_input, # trend on rolling-mean input
    trend_roll_anom    = trend_roll_anom,  # trend on rolling-mean anomalies
    trend_roll_z_anom  = trend_roll_z_anom,

    roll_mean_input    = roll_mean_input,  # rolling mean (input) via terra::roll()
    roll_mean_anom     = roll_mean_anom,   # rolling mean (anomalies)
    roll_mean_z_anom   = roll_mean_z_anom, # rolling mean (z-anomalies)

    trend_t0           = as.Date(round(mean(as.numeric(tt), na.rm=TRUE)),
                                 origin = "1970-01-01"),
    roll_times         = tt_roll
  )
}

anomalies_spatraster_fast <- function(input,
                                 baseline_start = as.Date("1981-01-01"),
                                 baseline_end   = as.Date("2010-12-31"),
                                 detrend = FALSE,
                                 # rolling options
                                 roll_window  = 12L,
                                 roll_align   = c("right","center","left"),
                                 roll_min_obs = NULL,
                                 # performance options
                                 cores = max(1L, parallel::detectCores() - 1L),
                                 float32 = TRUE) {

  stopifnot(inherits(input, "SpatRaster"))
  roll_align <- match.arg(roll_align)
  if (is.null(roll_min_obs)) roll_min_obs <- ceiling(roll_window/2)

  # ----- perf knobs (optional) -----
  if (float32) {
    terraOptions(datatype = "FLT4S")     # halve I/O vs FLT8S
  }
  terraOptions(progress = 1, memfrac = 0.8)

  # ----- helpers -----
  .get_year  <- function(d) as.POSIXlt(d)$year + 1900L
  .get_month <- function(d) as.POSIXlt(d)$mon  + 1L

  # Streaming OLS: for Y(t) and centered time x(t)
  # Returns c(intercept_at_x0, slope_per_unit_x)
  .ols_trend_stream <- function(Y, x_centered) {
    stopifnot(nlyr(Y) == length(x_centered))
    proto <- Y[[1]]
    zero  <- proto * 0

    sumy   <- zero
    sumx   <- zero
    sumxx  <- zero
    sumxy  <- zero
    nvalid <- zero

    for (i in seq_len(nlyr(Y))) {
      yi <- Y[[i]]
      vi <- !is.na(yi)               # valid mask
      yi0 <- ifel(vi, yi, 0)         # NA -> 0 for sums
      xi  <- x_centered[i]
      xi2 <- xi * xi

      # accumulate (all rasterized & compiled)
      nvalid <- nvalid + ifel(vi, 1, 0)
      sumy   <- sumy   + yi0
      sumx   <- sumx   + ifel(vi, xi, 0)     # only where valid
      sumxx  <- sumxx  + ifel(vi, xi2, 0)
      sumxy  <- sumxy  + (yi0 * xi)
    }

    den   <- nvalid * sumxx - sumx * sumx
    slope <- (nvalid * sumxy - sumx * sumy) / den
    slope <- ifel(den == 0, NA, slope)

    intercept <- (sumy - slope * sumx) / nvalid
    intercept <- ifel(nvalid == 0, NA, intercept)

    names(intercept) <- "intercept_t0"
    names(slope)     <- "slope_per_year"
    c(intercept, slope)
  }

  # Fast rolling mean using cumulative sums (NA-aware)
  .fast_roll_mean <- function(x, k, align = "right", min_obs = ceiling(k/2)) {
    n <- nlyr(x)
    if (k < 1L || k > n) return(rast())

    ones <- ifel(is.na(x), 0, 1)
    x0   <- ifel(is.na(x), 0, x)

    cs_x <- cumsum(x0)
    cs_n <- cumsum(ones)

    zero <- cs_x[[1]] * 0
    diff_layers <- function(cs, s, e) { (if (s > 1L) cs[[s-1L]] else zero); cs[[e]] - (if (s > 1L) cs[[s-1L]] else zero) }

    if (align == "right") {
      e <- k:n; s <- 1:(n-k+1L); stamp <- e
    } else if (align == "left") {
      s <- 1:(n-k+1L); e <- k:n; stamp <- s
    } else {
      hl <- floor((k-1L)/2L); hr <- k - hl - 1L
      s <- 1:(n-k+1L); e <- k:n; stamp <- (1+hl):(n-hr)
    }

    sum_list <- vector("list", length(e))
    n_list   <- vector("list", length(e))
    for (i in seq_along(e)) {
      sum_list[[i]] <- cs_x[[e[i]]] - (if (s[i] > 1L) cs_x[[s[i]-1L]] else zero)
      n_list[[i]]   <- cs_n[[e[i]]] - (if (s[i] > 1L) cs_n[[s[i]-1L]] else zero)
    }
    sum_k <- rast(sum_list)
    n_k   <- rast(n_list)

    m <- sum_k / n_k
    m <- ifel(n_k < min_obs, NA, m)

    tt <- time(x)
    names(m) <- paste0("rollmean_", format(tt[stamp], "%Y-%m"))
    time(m)  <- tt[stamp]
    m
  }

  # ----- time & baseline -----
  tt <- time(input)
  if (is.null(tt)) stop("`input` must have a time vector (terra::time(input)).")
  if (!inherits(tt, "Date")) tt <- as.Date(tt)
  input_raw <- input

  i_base <- which(!is.na(tt) & tt >= baseline_start & tt <= baseline_end)
  if (length(i_base) == 0) stop("No layers fall within the chosen baseline period.")
  rb <- input[[i_base]]
  tb <- tt[i_base]

  tn_full <- .get_year(tt) + (.get_month(tt) - 0.5)/12
  tn_base <- .get_year(tb) + (.get_month(tb) - 0.5)/12

  t0_base <- mean(tn_base, na.rm = TRUE)
  t0_full <- mean(tn_full, na.rm = TRUE)

  tn_c_base <- tn_base - t0_base
  tn_c_full <- tn_full - t0_full

  # ----- baseline trend on RAW (streaming OLS) -----
  trend_stack <- .ols_trend_stream(input_raw[[i_base]], tn_c_base)
  names(trend_stack) <- c("intercept_t0", "slope_per_year")

  # optional detrend
  if (isTRUE(detrend)) {
    input <- input - (trend_stack[[1]] + trend_stack[[2]] * (tn_full - t0_base))
  }

  # ----- monthly climatology & SD on baseline -----
  m_all  <- .get_month(tt)
  m_base <- .get_month(tb)

  clim_idx <- tapp(rb, index = m_base, fun = mean, na.rm = TRUE, cores = cores)
  sd_idx   <- tapp(rb, index = m_base, fun = function(x) stats::sd(x, na.rm = TRUE), cores = cores)

  present <- sort(unique(m_base))
  na_r <- rb[[1]] * NA_real_
  make_full12 <- function(stack_idx) {
    lst <- lapply(1:12, function(m) {
      pos <- match(m, present)
      if (is.na(pos)) na_r else stack_idx[[pos]]
    })
    out <- rast(lst); names(out) <- month.abb; out
  }
  clim12 <- make_full12(clim_idx)
  sd12   <- make_full12(sd_idx)

  # ----- anomalies & z-anomalies -----
  clim_expanded <- clim12[[m_all]]
  sd_expanded   <- sd12[[m_all]]

  anom  <- input - clim_expanded
  names(anom) <- paste0("anom_", names(input))

  z_anom <- (input - clim_expanded) / sd_expanded
  z_anom <- ifel(sd_expanded == 0, NA, z_anom)
  names(z_anom) <- paste0("z_", names(input))

  # ----- trends over FULL series (streaming OLS) -----
  trend_anom   <- .ols_trend_stream(anom,   tn_c_full)
  names(trend_anom) <- c("anom_intercept_t0", "anom_slope_per_year")

  trend_z_anom <- .ols_trend_stream(z_anom, tn_c_full)
  names(trend_z_anom) <- c("z_anom_intercept_t0", "z_anom_slope_per_year")

  # ----- rolling means & their trends -----
  roll_mean_input  <- .fast_roll_mean(input,  k = roll_window, align = roll_align, min_obs = roll_min_obs)
  roll_mean_anom   <- .fast_roll_mean(anom,   k = roll_window, align = roll_align, min_obs = roll_min_obs)
  roll_mean_z_anom <- .fast_roll_mean(z_anom, k = roll_window, align = roll_align, min_obs = roll_min_obs)

  if (nlyr(roll_mean_input) > 0) {
    tt_roll   <- time(roll_mean_input)
    tn_roll   <- .get_year(tt_roll) + (.get_month(tt_roll) - 0.5)/12
    t0_roll   <- mean(tn_roll, na.rm = TRUE)
    tn_c_roll <- tn_roll - t0_roll

    trend_roll_input  <- .ols_trend_stream(roll_mean_input,  tn_c_roll)
    names(trend_roll_input) <- c("roll_input_intercept_t0", "roll_input_slope_per_year")

    trend_roll_anom   <- .ols_trend_stream(roll_mean_anom,   tn_c_roll)
    names(trend_roll_anom)  <- c("roll_anom_intercept_t0",  "roll_anom_slope_per_year")

    trend_roll_z_anom <- .ols_trend_stream(roll_mean_z_anom, tn_c_roll)
    names(trend_roll_z_anom) <- c("roll_z_anom_intercept_t0", "roll_z_anom_slope_per_year")
  } else {
    trend_roll_input  <- NULL
    trend_roll_anom   <- NULL
    trend_roll_z_anom <- NULL
  }

  # ----- warning on missing months -----
  missing_months <- month.abb[!(1:12 %in% present)]
  if (length(missing_months)) {
    warning("Baseline missing months: ", paste(missing_months, collapse = ", "),
            ". Corresponding climatology/SD months are NA.")
  }

  # ----- return -----
  list(
    anom               = anom,
    z_anom             = z_anom,
    clim12             = clim12,
    sd12               = sd12,
    trend              = trend_stack,        # raw trend over baseline
    trend_anom         = trend_anom,         # trend on anomalies (full series)
    trend_z_anom       = trend_z_anom,       # trend on z-anomalies (full series)
    roll_mean_input    = roll_mean_input,
    trend_roll_input   = trend_roll_input,
    roll_mean_anom     = roll_mean_anom,
    trend_roll_anom    = trend_roll_anom,
    roll_mean_z_anom   = roll_mean_z_anom,
    trend_roll_z_anom  = trend_roll_z_anom,
    trend_t0           = as.Date(round(mean(as.numeric(tt), na.rm = TRUE)), origin = "1970-01-01"),
    roll_times         = if (nlyr(roll_mean_input) > 0) time(roll_mean_input) else as.Date(character())
  )
}

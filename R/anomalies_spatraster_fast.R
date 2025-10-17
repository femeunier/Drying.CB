anomalies_spatraster_fast <- function(
    input,
    baseline_start = as.Date("1981-01-01"),
    baseline_end   = as.Date("2010-12-31"),
    detrend = FALSE,
    roll_window  = 12L,
    roll_align   = c("right","center","left"),  # "right" recommended
    roll_min_obs = NULL,
    float32 = TRUE,
    tmpdir = NULL,            # optionally point to a fast SSD for temp files
    todisk = TRUE             # let terra spill to disk for large ops
) {
  stopifnot(inherits(input, "SpatRaster"))
  roll_align <- match.arg(roll_align)
  if (is.null(roll_min_obs)) roll_min_obs <- ceiling(roll_window/2)

  # --- Terra options (no parallel, no sockets) ---
  if (!is.null(tmpdir)) terraOptions(tempdir = tmpdir)
  if (todisk) terraOptions(todisk = TRUE)
  if (float32) terraOptions(datatype = "FLT4S")
  terraOptions(progress = 1, memfrac = 0.8)

  # --- Helpers ---
  .get_year  <- function(d) as.POSIXlt(d)$year + 1900L
  .get_month <- function(d) as.POSIXlt(d)$mon  + 1L

  # Streaming OLS across layers (intercept at x=0, slope per x-unit)
  .ols_trend_stream <- function(Y, x_centered) {
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

  # True streaming rolling mean (sliding window; O(1) memory w.r.t. layers)
  .roll_mean_stream <- function(x, k, align = "right", min_obs = ceiling(k/2)) {
    n <- nlyr(x); if (k < 1L || k > n) return(rast())
    tt <- time(x); if (!inherits(tt, "Date")) tt <- as.Date(tt)

    # Precompute per-layer valid-count rasters on the fly
    get_ones <- function(r) ifel(is.na(r), 0, 1)
    get_vals <- function(r) ifel(is.na(r), 0, r)

    # initialize window [1..k]
    sum_k <- x[[1]] * 0
    n_k   <- x[[1]] * 0
    for (i in 1:k) {
      sum_k <- sum_k + get_vals(x[[i]])
      n_k   <- n_k   + get_ones(x[[i]])
    }
    out_list <- list(sum_k / n_k)
    # stream remaining windows
    for (end in (k+1):n) {
      start <- end - k
      sum_k <- sum_k + get_vals(x[[end]])  - get_vals(x[[start]])
      n_k   <- n_k   + get_ones(x[[end]])  - get_ones(x[[start]])
      out_list[[length(out_list)+1L]] <- sum_k / n_k
    }
    out <- rast(out_list)

    # enforce minimum observations
    # rebuild n_k per output layer (cheap: reuse streaming with ones only)
    # We re-stream counts quickly:
    n_out <- list()
    n_w   <- x[[1]] * 0
    for (i in 1:k) n_w <- n_w + get_ones(x[[i]])
    n_out[[1]] <- n_w
    for (end in (k+1):n) {
      start <- end - k
      n_w <- n_w + get_ones(x[[end]]) - get_ones(x[[start]])
      n_out[[length(n_out)+1L]] <- n_w
    }
    n_out <- rast(n_out)
    out   <- ifel(n_out < min_obs, NA, out)

    # timestamps & names
    if (align == "right") {
      stamp <- k:n
    } else if (align == "left") {
      stamp <- 1:(n - k + 1L)
    } else { # center
      hl <- floor((k - 1L)/2L); hr <- k - hl - 1L
      stamp <- (1 + hl):(n - hr)
    }
    names(out) <- paste0("rollmean_", format(tt[stamp], "%Y-%m"))
    time(out)  <- tt[stamp]
    out
  }

  # --- Time & baseline prep ---
  tt <- time(input); if (is.null(tt)) stop("`input` must have a time vector.")
  if (!inherits(tt, "Date")) tt <- as.Date(tt)
  time(input) <- tt

  i_base <- which(!is.na(tt) & tt >= baseline_start & tt <= baseline_end)
  if (!length(i_base)) stop("No layers fall within the chosen baseline period.")
  rb <- input[[i_base]]; tb <- tt[i_base]

  tn_full <- .get_year(tt) + (.get_month(tt) - 0.5)/12
  tn_base <- .get_year(tb) + (.get_month(tb) - 0.5)/12
  t0_full <- mean(tn_full, na.rm = TRUE)
  t0_base <- mean(tn_base, na.rm = TRUE)
  tn_c_full <- tn_full - t0_full
  tn_c_base <- tn_base - t0_base

  # --- Baseline trend on RAW (streaming OLS) ---
  trend_raw <- .ols_trend_stream(input[[i_base]], tn_c_base)
  names(trend_raw) <- c("intercept_t0","slope_per_year")

  # optional detrend for anomaly calc
  if (isTRUE(detrend)) {
    input <- input - (trend_raw[[1]] + trend_raw[[2]] * (tn_full - t0_base))
  }

  # --- Monthly climatology & SD on baseline (compiled sums) ---
  m_all  <- .get_month(tt)
  m_base <- .get_month(tb)

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

  # --- Anomalies & z-anomalies ---
  clim_exp <- clim12[[m_all]]
  sd_exp   <- sd12[[m_all]]

  anom  <- input - clim_exp; names(anom)  <- paste0("anom_", names(input))
  z_anom <- (input - clim_exp) / sd_exp
  z_anom <- ifel(sd_exp == 0, NA, z_anom); names(z_anom) <- paste0("z_", names(input))

  # --- Trends on anomalies & z (streaming OLS) ---
  trend_anom   <- .ols_trend_stream(anom,   tn_c_full); names(trend_anom)   <- c("anom_intercept_t0","anom_slope_per_year")
  trend_z_anom <- .ols_trend_stream(z_anom, tn_c_full); names(trend_z_anom) <- c("z_anom_intercept_t0","z_anom_slope_per_year")

  # --- Rolling means (streaming) & their trends ---
  r_input <- .roll_mean_stream(input,  k = roll_window, align = roll_align, min_obs = roll_min_obs)
  r_anom  <- .roll_mean_stream(anom,   k = roll_window, align = roll_align, min_obs = roll_min_obs)
  r_zanom <- .roll_mean_stream(z_anom, k = roll_window, align = roll_align, min_obs = roll_min_obs)

  if (nlyr(r_input) > 0) {
    tt_r   <- time(r_input)
    tnr    <- .get_year(tt_r) + (.get_month(tt_r) - 0.5)/12
    t0r    <- mean(tnr, na.rm = TRUE)
    tnr_c  <- tnr - t0r

    trend_r_input <- .ols_trend_stream(r_input,  tnr_c); names(trend_r_input) <- c("roll_input_intercept_t0","roll_input_slope_per_year")
    trend_r_anom  <- .ols_trend_stream(r_anom,   tnr_c); names(trend_r_anom)  <- c("roll_anom_intercept_t0","roll_anom_slope_per_year")
    trend_r_zanom <- .ols_trend_stream(r_zanom,  tnr_c); names(trend_r_zanom) <- c("roll_z_anom_intercept_t0","roll_z_anom_slope_per_year")
  } else {
    trend_r_input <- trend_r_anom <- trend_r_zanom <- NULL
  }

  # --- Warning if baseline misses months ---
  missing_months <- month.abb[!(1:12 %in% present)]
  if (length(missing_months)) warning("Baseline missing months: ", paste(missing_months, collapse = ", "), ". Corresponding climatology/SD months are NA.")

  # --- Return ---
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

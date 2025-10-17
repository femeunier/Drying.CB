anomalies_spatraster <- function(input,
                                 baseline_start = as.Date("1981-01-01"),
                                 baseline_end   = as.Date("2010-12-31"),
                                 detrend = FALSE,
                                 roll_window  = 12L,
                                 roll_align   = c("right","center","left"),
                                 roll_min_obs = NULL) {

  roll_align <- match.arg(roll_align)
  if (is.null(roll_min_obs)) roll_min_obs <- ceiling(roll_window/2)

  # helpers ---------------------------------------------------------------
  get_year  <- function(d) as.POSIXlt(d)$year + 1900L
  get_month <- function(d) as.POSIXlt(d)$mon  + 1L

  # If trend_fit is not in scope, define a safe default
  if (!exists("trend_fit", mode = "function")) {
    trend_fit <- function(v, tn, t0) {
      if (all(!is.finite(v))) return(c(NA_real_, NA_real_))
      ok <- is.finite(v) & is.finite(tn)
      if (sum(ok) < 2) return(c(NA_real_, NA_real_))
      X <- cbind(1, tn[ok] - t0)
      beta <- tryCatch(qr.solve(X, v[ok]), error = function(e) c(NA_real_, NA_real_))
      if (any(!is.finite(beta))) return(c(NA_real_, NA_real_))
      beta  # c(intercept_at_t0, slope_per_year)
    }
  }

  # Build rolling index list for terra::tapp over layers (time)
  build_roll_indices <- function(n, k, align=c("right","center","left")) {
    align <- match.arg(align)
    if (k < 1L || k > n) return(list())
    if (align == "right") {
      lapply(seq_len(n - k + 1L), function(i) i:(i + k - 1L))
    } else if (align == "left") {
      lapply(seq_len(n - k + 1L), function(i) i:(i + k - 1L))  # same windows;
      # we'll align timestamps differently below
    } else { # center
      # use odd/even handling by centering floor(k/2) to the left, remainder to the right
      half_left  <- floor((k - 1L)/2L)
      half_right <- k - half_left - 1L
      starts <- (1L + 0L):(n - k + 1L)
      lapply(starts, function(s) s:(s + k - 1L))
    }
  }

  # Rolling timestamps aligned to the chosen side (we’ll pick the "index time" per window)
  roll_times <- function(tt, k, align=c("right","center","left")) {
    align <- match.arg(align)
    n <- length(tt)
    if (k < 1L || k > n) return(as.Date(character()))
    if (align == "right") {
      tt[k:n]
    } else if (align == "left") {
      tt[1:(n - k + 1L)]
    } else { # center
      half_left  <- floor((k - 1L)/2L)
      idx <- (1L + half_left):(n - (k - half_left - 1L))
      tt[idx]
    }
  }

  # Time vector -----------------------------------------------------------
  tt <- time(input)
  if (is.null(tt)) stop("`input` must have a time vector (terra::time(input)).")
  if (!inherits(tt, "Date")) tt <- as.Date(tt)
  nlay <- nlyr(input)

  # Keep a copy of the raw series for trend estimation on raw and rolling
  input_raw <- input

  # Baseline subset -------------------------------------------------------
  i_base <- which(!is.na(tt) & tt >= baseline_start & tt <= baseline_end)
  if (length(i_base) == 0) stop("No layers fall within the chosen baseline period.")
  rb <- input[[i_base]]
  tb <- tt[i_base]

  # numeric time in years for regressions (full & baseline)
  tn        <- get_year(tt) + (get_month(tt) - 0.5)/12
  tn_years  <- get_year(tb) + (get_month(tb) - 0.5)/12
  t0_years  <- mean(tn_years, na.rm = TRUE)

  # --- linear trend on RAW data (two-layer raster: intercept_at_t0, slope_per_year)
  trend_stack <- app(input_raw[[i_base]],
                     fun = trend_fit,
                     tn = tn_years, t0 = t0_years)
  names(trend_stack) <- c("intercept_t0", "slope_per_year")

  # --- optional detrend for anomaly computation --------------------------
  if (isTRUE(detrend)) {
    # broadcast tn across layers via recycling inside terra
    input <- input - (trend_stack[[1]] + trend_stack[[2]] * tn)
  }

  # Month indices without format()
  m_all  <- get_month(tt)
  m_base <- get_month(tb)

  # Monthly climatology & SD via tapp -------------------------------------
  clim_idx <- tapp(rb,
                   index = m_base,
                   fun = function(x) mean(x, na.rm = TRUE))
  sd_idx   <- tapp(rb,
                   index = m_base,
                   fun = function(x) stats::sd(x, na.rm = TRUE))

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

  # Expand to full timeline and compute anomalies -------------------------
  clim_expanded <- clim12[[m_all]]
  sd_expanded   <- sd12[[m_all]]

  anom <- input - clim_expanded
  names(anom) <- paste0("anom_", names(input))

  z_anom <- (input - clim_expanded) / sd_expanded
  z_anom <- ifel(sd_expanded == 0, NA, z_anom)
  names(z_anom) <- paste0("z_", names(input))

  # Trend on anomalies (over the full series) -----------------------------
  trend_anom <- app(anom,
                    fun = trend_fit,
                    tn  = tn,
                    t0  = mean(tn, na.rm = TRUE))
  names(trend_anom) <- c("anom_intercept_t0", "anom_slope_per_year")

  # Rolling means (input & anomalies) via tapp ----------------------------
  roll_idx_list <- build_roll_indices(nlay, roll_window, roll_align)
  tt_roll <- roll_times(tt, roll_window, roll_align)
  tn_roll <- get_year(tt_roll) + (get_month(tt_roll) - 0.5)/12

  if (length(roll_idx_list)) {
    roll_mean_input <- tapp(input, index = roll_idx_list, fun = function(x) {
      ok <- is.finite(x)
      if (sum(ok) < roll_min_obs) return(NA_real_)
      mean(x[ok])
    })
    names(roll_mean_input) <- paste0("rollmean_", format(tt_roll, "%Y-%m"))

    roll_mean_anom <- tapp(anom, index = roll_idx_list, fun = function(x) {
      ok <- is.finite(x)
      if (sum(ok) < roll_min_obs) return(NA_real_)
      mean(x[ok])
    })
    names(roll_mean_anom) <- paste0("rollmean_anom_", format(tt_roll, "%Y-%m"))

    # Trend on rolling-mean anomalies -------------------------------------
    trend_roll_anom <- app(roll_mean_anom,
                           fun = trend_fit,
                           tn  = tn_roll,
                           t0  = mean(tn_roll, na.rm = TRUE))
    names(trend_roll_anom) <- c("roll_anom_intercept_t0", "roll_anom_slope_per_year")
  } else {
    roll_mean_input <- NULL
    roll_mean_anom  <- NULL
    trend_roll_anom <- NULL
  }

  # warn if baseline misses months
  missing_months <- month.abb[!(1:12 %in% present)]
  if (length(missing_months)) {
    warning("Baseline missing months: ", paste(missing_months, collapse = ", "),
            ". Corresponding climatology/SD months are NA.")
  }

  list(
    anom         = anom,
    z_anom       = z_anom,

    clim12       = clim12,
    sd12         = sd12,

    roll_mean_input = roll_mean_input,  # rolling mean of input
    roll_mean_anom  = roll_mean_anom,   # rolling mean of anomalies

    trend        = trend_stack,         # raw trend over baseline (intercept_t0, slope_per_year)
    trend_anom   = trend_anom,          # trend on anomalies over full series
    trend_roll_anom = trend_roll_anom,  # trend on rolling-mean anomalies

    trend_t0     = as.Date(round(mean(as.numeric(tt), na.rm=TRUE)),
                           origin = "1970-01-01"),
    roll_times   = tt_roll            # dates associated with rolling outputs
  )
}

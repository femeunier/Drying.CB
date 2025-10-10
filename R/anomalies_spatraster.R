anomalies_spatraster <- function(input,
                                 baseline_start = as.Date("1981-01-01"),
                                 baseline_end   = as.Date("2010-12-31"),
                                 detrend = FALSE) {

  # require time stamps
  t <- time(input)
  if (is.null(t)) stop("`input` must have a time vector (use terra::time(input)).")

  # --- optional linear de-trend per cell (remove slope, keep mean level) ---
  if (isTRUE(detrend)) {
    # numeric time in years (better conditioning than raw days)
    tn <- as.numeric(t) / 365.2425
    input <- app(input, fun = detrend_linear, tn = tn)
  }

  # select baseline window
  i_base <- which(!is.na(t) & t >= baseline_start & t <= baseline_end)
  if (length(i_base) == 0) stop("No layers fall within the chosen baseline period.")
  rb <- input[[i_base]]

  # month indices
  m_all  <- as.integer(format(t, "%m"))
  m_base <- as.integer(format(time(rb), "%m"))

  # monthly climatology (Jan..Dec order; robust to missing months)
  clim_list <- lapply(1:12, function(m) {
    lyr <- which(m_base == m)
    if (length(lyr) == 0) {
      app(rb[[1]], fun = function(x) NA_real_)
    } else {
      mean(rb[[lyr]], na.rm = TRUE)
    }
  })
  clim12 <- rast(clim_list); names(clim12) <- month.abb

  # expand climatology to full timeline
  clim_expanded <- clim12[[m_all]]

  # absolute anomalies
  anom <- input - clim_expanded
  names(anom) <- paste0("anom_", names(input))

  # monthly SD for z-scores
  sd_list <- lapply(1:12, function(m) {
    lyr <- which(m_base == m)
    if (length(lyr) == 0) {
      app(rb[[1]], fun = function(x) NA_real_)
    } else {
      stdev(rb[[lyr]], na.rm = TRUE)
    }
  })
  sd12 <- rast(sd_list); names(sd12) <- month.abb
  sd_expanded <- sd12[[m_all]]

  # z-anomalies, protect against sd = 0
  z_anom <- (input - clim_expanded) / sd_expanded
  z_anom <- ifel(sd_expanded == 0, NA, z_anom)
  names(z_anom) <- paste0("z_", names(input))

  # warn if some months are entirely missing in baseline
  missing_months <- month.abb[which(sapply(1:12, function(k) sum(m_base == k)) == 0)]
  if (length(missing_months)) {
    warning("Baseline missing months: ", paste(missing_months, collapse = ", "),
            ". Corresponding anomalies will be NA for those months.")
  }

  list(
    anom   = anom,
    z_anom = z_anom,
    clim12 = clim12,
    sd12   = sd12
  )
}

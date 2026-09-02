# Blockwise replacement for Drying.CB::anomalies_spatraster_roll()
#
# The function reads a few raster rows at a time, performs all time-series
# calculations in ordinary matrices, and writes every requested result directly
# to a file. Peak RAM therefore depends on block_mem_mb rather than on the full
# number of raster cells.

anomalies_spatraster_roll_blockwise <- function(
    input,
    baseline_start = as.Date("1981-01-01"),
    baseline_end   = as.Date("2010-12-31"),
    detrend = FALSE,
    roll_window  = 12L,
    roll_align   = c("right", "center", "left"),
    roll_min_obs = NULL,
    SD_threshold = NULL,
    Z_anom_threshold = NULL,
    output_dir = tempdir(),
    prefix = "anomalies",
    outputs = c(
      "anom", "z_anom", "clim12", "sd12", "trend",
      "trend_anom", "trend_z_anom",
      "trend_roll_input", "trend_roll_anom", "trend_roll_z_anom",
      "roll_mean_input", "roll_mean_anom", "roll_mean_z_anom"
    ),
    filenames = NULL,
    overwrite = FALSE,
    block_mem_mb = 512,
    block_nrows = NULL,
    datatype = "FLT4S",
    gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES", "BIGTIFF=YES"),
    finalize_retries = 30L,
    finalize_wait_seconds = 2,
    progress = TRUE) {

  if (!inherits(input, "SpatRaster")) {
    stop("`input` must be a terra SpatRaster.")
  }

  roll_align <- match.arg(roll_align)
  roll_window <- as.integer(roll_window)

  if (length(roll_window) != 1L || is.na(roll_window) || roll_window < 2L) {
    stop("`roll_window` must be one integer >= 2.")
  }

  if (is.null(roll_min_obs)) {
    roll_min_obs <- ceiling(roll_window / 2)
  }
  roll_min_obs <- as.integer(roll_min_obs)

  if (length(roll_min_obs) != 1L || is.na(roll_min_obs) ||
      roll_min_obs < 1L || roll_min_obs > roll_window) {
    stop("`roll_min_obs` must be between 1 and `roll_window`.")
  }

  if (!is.null(SD_threshold)) {
    if (length(SD_threshold) != 1L || !is.finite(SD_threshold) ||
        SD_threshold < 0) {
      stop("`SD_threshold` must be NULL or one non-negative number.")
    }
  }

  if (!is.null(Z_anom_threshold)) {
    if (length(Z_anom_threshold) != 1L || !is.finite(Z_anom_threshold) ||
        Z_anom_threshold <= 0) {
      stop("`Z_anom_threshold` must be NULL or one positive number.")
    }
  }

  finalize_retries <- as.integer(finalize_retries)
  if (length(finalize_retries) != 1L || is.na(finalize_retries) ||
      finalize_retries < 1L) {
    stop("`finalize_retries` must be a positive integer.")
  }
  if (length(finalize_wait_seconds) != 1L ||
      !is.finite(finalize_wait_seconds) || finalize_wait_seconds < 0) {
    stop("`finalize_wait_seconds` must be one non-negative number.")
  }

  all_outputs <- c(
    "anom", "z_anom", "clim12", "sd12", "trend",
    "trend_anom", "trend_z_anom",
    "trend_roll_input", "trend_roll_anom", "trend_roll_z_anom",
    "roll_mean_input", "roll_mean_anom", "roll_mean_z_anom"
  )

  if (length(outputs) == 1L && identical(outputs, "all")) {
    outputs <- all_outputs
  }
  outputs <- unique(as.character(outputs))

  unknown_outputs <- setdiff(outputs, all_outputs)
  if (length(unknown_outputs)) {
    stop("Unknown `outputs`: ", paste(unknown_outputs, collapse = ", "))
  }
  if (!length(outputs)) stop("At least one output must be requested.")

  tt <- terra::time(input)
  if (is.null(tt) || length(tt) != terra::nlyr(input)) {
    stop("`input` must have one time value per raster layer.")
  }
  if (!inherits(tt, "Date")) tt <- as.Date(tt)
  if (anyNA(tt)) stop("The input time vector contains missing dates.")

  n_time <- terra::nlyr(input)
  if (roll_window > n_time) {
    stop("`roll_window` cannot exceed the number of input layers.")
  }

  get_year  <- function(d) as.POSIXlt(d)$year + 1900L
  get_month <- function(d) as.POSIXlt(d)$mon + 1L

  years  <- get_year(tt)
  months <- get_month(tt)
  tn <- years + (months - 0.5) / 12

  i_base <- which(tt >= baseline_start & tt <= baseline_end)
  if (!length(i_base)) {
    stop("No layers fall within the selected baseline period.")
  }

  months_base <- months[i_base]
  tn_base <- tn[i_base]
  t0_base <- mean(tn_base)
  t0_full <- mean(tn)
  present_months <- sort(unique(months_base))

  missing_months <- month.abb[!(seq_len(12) %in% present_months)]
  if (length(missing_months)) {
    warning(
      "Baseline missing months: ", paste(missing_months, collapse = ", "),
      ". Corresponding climatology and SD layers will be NA."
    )
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) stop("Could not create `output_dir`.")

  # A named `filenames` vector/list can override any automatically generated
  # path. This is useful when the caller needs an established naming scheme.
  if (is.null(filenames)) filenames <- character()
  filenames <- unlist(filenames, use.names = TRUE)
  if (length(filenames) && (is.null(names(filenames)) || any(!nzchar(names(filenames))))) {
    stop("`filenames` must be named with output identifiers such as `anom`.")
  }

  output_paths <- setNames(character(length(outputs)), outputs)
  for (key in outputs) {
    if (key %in% names(filenames)) {
      output_paths[[key]] <- filenames[[key]]
    } else {
      output_paths[[key]] <- file.path(output_dir, paste0(prefix, "_", key, ".tif"))
    }
    dir.create(dirname(output_paths[[key]]), recursive = TRUE, showWarnings = FALSE)
  }

  existing <- output_paths[file.exists(output_paths)]
  if (length(existing) && !isTRUE(overwrite)) {
    stop(
      "Output file(s) already exist and `overwrite = FALSE`: ",
      paste(existing, collapse = ", ")
    )
  }

  layer_spec <- list(
    anom = list(
      n = n_time,
      names = paste0("anom_", names(input)),
      time = tt
    ),
    z_anom = list(
      n = n_time,
      names = paste0("z_", names(input)),
      time = tt
    ),
    clim12 = list(n = 12L, names = month.abb, time = NULL),
    sd12   = list(n = 12L, names = month.abb, time = NULL),
    trend = list(
      n = 2L,
      names = c("intercept_t0", "slope_per_year"),
      time = NULL
    ),
    trend_anom = list(
      n = 2L,
      names = c("anom_intercept_t0", "anom_slope_per_year"),
      time = NULL
    ),
    trend_z_anom = list(
      n = 2L,
      names = c("z_anom_intercept_t0", "z_anom_slope_per_year"),
      time = NULL
    ),
    trend_roll_input = list(
      n = 2L,
      names = c("roll_input_intercept_t0", "roll_input_slope_per_year"),
      time = NULL
    ),
    trend_roll_anom = list(
      n = 2L,
      names = c("roll_anom_intercept_t0", "roll_anom_slope_per_year"),
      time = NULL
    ),
    trend_roll_z_anom = list(
      n = 2L,
      names = c("roll_z_anom_intercept_t0", "roll_z_anom_slope_per_year"),
      time = NULL
    ),
    roll_mean_input = list(
      n = n_time,
      names = paste0("rollmean_", format(tt, "%Y-%m")),
      time = tt
    ),
    roll_mean_anom = list(
      n = n_time,
      names = paste0("rollmean_", format(tt, "%Y-%m")),
      time = tt
    ),
    roll_mean_z_anom = list(
      n = n_time,
      names = paste0("rollmean_", format(tt, "%Y-%m")),
      time = tt
    )
  )

  make_template <- function(n_layers, layer_names, layer_time = NULL) {
    r <- terra::rast(
      nrows = terra::nrow(input),
      ncols = terra::ncol(input),
      nlyrs = n_layers,
      xmin = terra::xmin(input),
      xmax = terra::xmax(input),
      ymin = terra::ymin(input),
      ymax = terra::ymax(input),
      crs = terra::crs(input)
    )
    names(r) <- layer_names
    if (!is.null(layer_time)) terra::time(r) <- layer_time
    r
  }

  wopt <- list(datatype = datatype, gdal = gdal)
  writers <- setNames(vector("list", length(outputs)), outputs)

  for (key in outputs) {
    spec <- layer_spec[[key]]
    r <- make_template(spec$n, spec$names, spec$time)
    terra::writeStart(
      r,
      filename = output_paths[[key]],
      overwrite = overwrite,
      wopt = wopt
    )
    writers[[key]] <- list(raster = r, open = TRUE)
  }

  input_open <- FALSE
  completed <- FALSE

  on.exit({
    if (input_open) try(terra::readStop(input), silent = TRUE)
    if (!completed && exists("writers", inherits = FALSE)) {
      for (key in names(writers)) {
        if (isTRUE(writers[[key]]$open)) {
          try(terra::writeStop(writers[[key]]$raster), silent = TRUE)
        }
      }
    }
  }, add = TRUE)

  terra::readStart(input)
  input_open <- TRUE

  # Roughly eight full block matrices can coexist at the local peak. The
  # estimate deliberately errs on the conservative side.
  if (is.null(block_nrows)) {
    bytes_available <- block_mem_mb * 1024^2
    cells_per_block <- floor(bytes_available / (8 * n_time * 8))
    block_nrows <- floor(cells_per_block / terra::ncol(input))
    block_nrows <- max(1L, block_nrows)
  }
  block_nrows <- min(as.integer(block_nrows), terra::nrow(input))
  if (is.na(block_nrows) || block_nrows < 1L) {
    stop("`block_nrows` must be NULL or a positive integer.")
  }

  starts <- seq.int(1L, terra::nrow(input), by = block_nrows)
  n_blocks <- length(starts)

  write_block <- function(key, value, row, nrows) {
    if (!(key %in% names(writers))) return(invisible(NULL))
    if (!is.matrix(value)) value <- as.matrix(value)
    terra::writeValues(writers[[key]]$raster, value, row, nrows)
    invisible(NULL)
  }

  # Row-wise mean and sample SD, ignoring all non-finite observations. Climate
  # rasters should not contain Inf; treating it as missing prevents an Inf from
  # contaminating every derived layer for that pixel.
  row_mean_sd <- function(x) {
    if (!is.matrix(x)) x <- as.matrix(x)
    ok <- is.finite(x)
    n <- rowSums(ok)
    x[!ok] <- 0
    sx <- rowSums(x)
    sx2 <- rowSums(x * x)

    avg <- sx / n
    avg[n == 0] <- NA_real_

    variance <- (sx2 - (sx * sx) / n) / (n - 1)
    variance[n < 2] <- NA_real_
    variance[is.finite(variance) & variance < 0] <- 0

    list(mean = avg, sd = sqrt(variance))
  }

  # Vectorized OLS for every matrix row. The returned intercept is evaluated at
  # t0, matching the original trend_fit() parameterization.
  fit_trend_matrix <- function(y, x, t0) {
    if (!is.matrix(y)) y <- as.matrix(y)
    ok <- is.finite(y)
    y[!ok] <- 0

    n <- rowSums(ok)
    sx <- as.vector(ok %*% x)
    sxx <- as.vector(ok %*% (x * x))
    sy <- rowSums(y)
    sxy <- as.vector(y %*% x)

    denominator <- n * sxx - sx * sx
    tolerance <- .Machine$double.eps * pmax(1, abs(n * sxx), abs(sx * sx))
    valid <- n >= 2 & is.finite(denominator) & abs(denominator) > tolerance

    slope <- rep(NA_real_, nrow(y))
    intercept <- rep(NA_real_, nrow(y))

    slope[valid] <- (
      n[valid] * sxy[valid] - sx[valid] * sy[valid]
    ) / denominator[valid]

    intercept[valid] <- (
      sy[valid] - slope[valid] * (sx[valid] - n[valid] * t0)
    ) / n[valid]

    cbind(intercept, slope)
  }

  make_anomaly <- function(x, climatology) {
    out <- x
    for (m in seq_len(12)) {
      j <- which(months == m)
      if (length(j)) {
        out[, j] <- sweep(x[, j, drop = FALSE], 1L, climatology[, m], "-")
      }
    }
    out
  }

  make_z_anomaly <- function(anomaly, standard_deviation) {
    out <- anomaly

    for (m in seq_len(12)) {
      j <- which(months == m)
      if (!length(j)) next

      denom <- standard_deviation[, m]
      valid <- is.finite(denom) & denom != 0
      if (!is.null(SD_threshold)) valid <- valid & denom >= SD_threshold

      out[, j] <- sweep(anomaly[, j, drop = FALSE], 1L, denom, "/")
      if (any(!valid)) out[!valid, j] <- NA_real_
    }

    if (!is.null(Z_anom_threshold)) {
      too_large <- is.finite(out) & abs(out) > Z_anom_threshold
      out[too_large] <- NA_real_
    }

    out
  }

  # Full windows are required at the time-series boundaries, as in
  # terra::roll(..., circular = FALSE). `right` means a trailing window,
  # `left` a forward window. For an even centered window, the extra position is
  # placed after the focal layer.
  roll_matrix <- function(x) {
    out <- matrix(NA_real_, nrow = nrow(x), ncol = ncol(x))

    if (roll_align == "right") {
      before <- roll_window - 1L
      after <- 0L
    } else if (roll_align == "left") {
      before <- 0L
      after <- roll_window - 1L
    } else {
      before <- floor((roll_window - 1L) / 2L)
      after <- (roll_window - 1L) - before
    }

    first_i <- 1L + before
    last_i <- ncol(x) - after

    if (first_i > last_i) return(out)

    for (i in seq.int(first_i, last_i)) {
      j <- seq.int(i - before, i + after)
      window_values <- x[, j, drop = FALSE]
      ok <- !is.na(window_values)
      count <- rowSums(ok)
      window_values[!ok] <- 0
      value <- rowSums(window_values) / count
      value[count < roll_min_obs | count == 0] <- NA_real_
      out[, i] <- value
    }

    out
  }

  need_roll_input <- any(c("roll_mean_input", "trend_roll_input") %in% outputs)
  need_roll_anom <- any(c("roll_mean_anom", "trend_roll_anom") %in% outputs)
  need_roll_z <- any(c("roll_mean_z_anom", "trend_roll_z_anom") %in% outputs)

  if (isTRUE(progress)) {
    message(
      "Processing ", n_blocks, " spatial blocks (", block_nrows,
      " raster row(s) per block)."
    )
  }

  for (iblock in seq_along(starts)) {
    row_start <- starts[iblock]
    n_rows <- min(block_nrows, terra::nrow(input) - row_start + 1L)

    if (isTRUE(progress)) {
      message("Block ", iblock, "/", n_blocks, ": rows ", row_start,
              "-", row_start + n_rows - 1L)
    }

    raw <- terra::readValues(
      input,
      row = row_start,
      nrows = n_rows,
      mat = TRUE
    )
    if (!is.matrix(raw)) raw <- matrix(raw, ncol = n_time)

    # Raw baseline trend.
    trend_raw <- fit_trend_matrix(raw[, i_base, drop = FALSE], tn_base, t0_base)
    write_block("trend", trend_raw, row_start, n_rows)

    # Correct detrending: fitted values are intercept_t0 + slope * (t - t0).
    # The original function omitted `- t0`, which only changed the level but
    # made the returned detrended rolling-input intercept difficult to interpret.
    work <- raw
    if (isTRUE(detrend)) {
      for (j in seq_len(n_time)) {
        fitted_j <- trend_raw[, 1L] + trend_raw[, 2L] * (tn[j] - t0_base)
        work[, j] <- raw[, j] - fitted_j
      }
    }

    # With detrending enabled, calculate baseline statistics from the detrended
    # series. With detrending disabled this is identical to the original code.
    climatology <- matrix(NA_real_, nrow = nrow(work), ncol = 12L)
    standard_deviation <- matrix(NA_real_, nrow = nrow(work), ncol = 12L)

    for (m in present_months) {
      j <- i_base[months_base == m]
      stats_m <- row_mean_sd(work[, j, drop = FALSE])
      climatology[, m] <- stats_m$mean
      standard_deviation[, m] <- stats_m$sd
    }

    write_block("clim12", climatology, row_start, n_rows)
    write_block("sd12", standard_deviation, row_start, n_rows)

    rm(raw)

    if (need_roll_input) {
      rolled <- roll_matrix(work)
      write_block("roll_mean_input", rolled, row_start, n_rows)
      if ("trend_roll_input" %in% outputs) {
        write_block(
          "trend_roll_input",
          fit_trend_matrix(rolled, tn, t0_full),
          row_start,
          n_rows
        )
      }
      rm(rolled)
    }

    anomaly <- make_anomaly(work, climatology)
    write_block("anom", anomaly, row_start, n_rows)

    if ("trend_anom" %in% outputs) {
      write_block(
        "trend_anom",
        fit_trend_matrix(anomaly, tn, t0_full),
        row_start,
        n_rows
      )
    }

    if (need_roll_anom) {
      rolled <- roll_matrix(anomaly)
      write_block("roll_mean_anom", rolled, row_start, n_rows)
      if ("trend_roll_anom" %in% outputs) {
        write_block(
          "trend_roll_anom",
          fit_trend_matrix(rolled, tn, t0_full),
          row_start,
          n_rows
        )
      }
      rm(rolled)
    }

    z_anomaly <- make_z_anomaly(anomaly, standard_deviation)
    write_block("z_anom", z_anomaly, row_start, n_rows)

    if ("trend_z_anom" %in% outputs) {
      write_block(
        "trend_z_anom",
        fit_trend_matrix(z_anomaly, tn, t0_full),
        row_start,
        n_rows
      )
    }

    rm(anomaly)

    if (need_roll_z) {
      rolled <- roll_matrix(z_anomaly)
      write_block("roll_mean_z_anom", rolled, row_start, n_rows)
      if ("trend_roll_z_anom" %in% outputs) {
        write_block(
          "trend_roll_z_anom",
          fit_trend_matrix(rolled, tn, t0_full),
          row_start,
          n_rows
        )
      }
      rm(rolled)
    }

    rm(work, climatology, standard_deviation, trend_raw, z_anomaly)
    gc(verbose = FALSE)
  }

  terra::readStop(input)
  input_open <- FALSE

  # On parallel filesystems, writeStop() can finish the GDAL write but fail
  # while immediately reopening the new GeoTIFF because its directory entry is
  # not visible yet. Close every writer first and defer reopening until all
  # handles have been finalized.
  close_errors <- setNames(vector("list", length(writers)), names(writers))

  for (key in names(writers)) {
    closed <- try(terra::writeStop(writers[[key]]$raster), silent = TRUE)
    writers[[key]]$open <- FALSE

    if (inherits(closed, "try-error")) {
      close_errors[[key]] <- conditionMessage(attr(closed, "condition"))
    }
  }

  # Release the writer objects and flush filesystem/GDAL state before testing
  # whether the completed files can be opened.
  rm(writers)
  gc(verbose = FALSE)

  reopen_with_retry <- function(path) {
    last_error <- NULL

    for (attempt in seq_len(finalize_retries)) {
      candidate <- try(terra::rast(path), silent = TRUE)

      if (!inherits(candidate, "try-error")) {
        return(candidate)
      }

      last_error <- conditionMessage(attr(candidate, "condition"))

      if (attempt < finalize_retries && finalize_wait_seconds > 0) {
        Sys.sleep(finalize_wait_seconds)
      }
    }

    stop(
      "Output could not be opened after ", finalize_retries,
      " attempts: ", path,
      if (!is.null(last_error)) paste0("\nLast GDAL/terra error: ", last_error)
    )
  }

  result <- setNames(vector("list", length(all_outputs)), all_outputs)
  for (key in names(output_paths)) {
    result[[key]] <- reopen_with_retry(output_paths[[key]])

    if (!is.null(close_errors[[key]]) && isTRUE(progress)) {
      message(
        "Output finalized after a transient writeStop() error: ",
        basename(output_paths[[key]])
      )
    }
  }

  completed <- TRUE

  result$trend_t0 <- as.Date(
    round(mean(as.numeric(tt))),
    origin = "1970-01-01"
  )
  result$roll_times <- tt
  result$files <- output_paths
  result
}



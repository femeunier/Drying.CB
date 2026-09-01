#!/usr/bin/env Rscript

rm(list = ls())

suppressPackageStartupMessages({
  library(ncdf4)
  library(data.table)
  library(matrixStats)
})

# ---------------------------------------------------------------------------
# User settings
# ---------------------------------------------------------------------------

years <- 1970:2026

input_dir <- paste0(
  "/data/gent/vo/000/gvo00074/",
  "ED_common_data/met/Tropics"
)

output_file <- "./outputs/df.ERA5.Tropics.climate.RDS"

temperature_pattern <- "ERA5_tropics_hourly_temperature_%04d.nc"
precipitation_pattern <- "ERA5_tropics_hourly_precipitation_%04d.nc"

dir.create(
  dirname(output_file),
  recursive = TRUE,
  showWarnings = FALSE
)

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

find_netcdf_name <- function(nc, candidates, description, path) {

  available_names <- unique(
    c(
      names(nc$dim),
      names(nc$var)
    )
  )

  selected_name <- candidates[candidates %in% available_names][1]

  if (
    length(selected_name) == 0L ||
    is.na(selected_name)
  ) {
    stop(
      "Could not identify the ",
      description,
      " in ",
      path,
      ". Tried: ",
      paste(candidates, collapse = ", ")
    )
  }

  selected_name
}


require_variable <- function(nc, variable, path) {

  if (!variable %in% names(nc$var)) {
    stop(
      "Variable '",
      variable,
      "' is absent from ",
      path,
      ".\nAvailable data variables: ",
      paste(names(nc$var), collapse = ", ")
    )
  }
}


get_grid_info <- function(nc, path) {

  longitude_name <- find_netcdf_name(
    nc,
    candidates = c("longitude", "lon"),
    description = "longitude coordinate",
    path = path
  )

  latitude_name <- find_netcdf_name(
    nc,
    candidates = c("latitude", "lat"),
    description = "latitude coordinate",
    path = path
  )

  time_name <- find_netcdf_name(
    nc,
    candidates = c("valid_time", "time"),
    description = "time coordinate",
    path = path
  )

  list(
    longitude_name = longitude_name,
    latitude_name = latitude_name,
    time_name = time_name,
    longitude = ncvar_get(nc, longitude_name),
    latitude = ncvar_get(nc, latitude_name)
  )
}


parse_netcdf_time <- function(nc, time_name, path) {

  time_values <- ncvar_get(nc, time_name)

  units_attribute <- ncatt_get(
    nc,
    time_name,
    "units"
  )

  time_units <- if (isTRUE(units_attribute$hasatt)) {
    units_attribute$value
  } else if (
    time_name %in% names(nc$dim) &&
    !is.null(nc$dim[[time_name]]$units)
  ) {
    nc$dim[[time_name]]$units
  } else {
    NA_character_
  }

  if (
    length(time_units) != 1L ||
    is.na(time_units) ||
    !nzchar(time_units)
  ) {
    stop(
      "Time units are missing from '",
      time_name,
      "' in ",
      path
    )
  }

  unit_match <- regexec(
    paste0(
      "^(seconds?|minutes?|hours?|days?)\\s+since\\s+",
      "(.+)$"
    ),
    time_units,
    ignore.case = TRUE
  )

  unit_parts <- regmatches(
    time_units,
    unit_match
  )[[1]]

  if (length(unit_parts) != 3L) {
    stop(
      "Unsupported NetCDF time units in ",
      path,
      ": ",
      time_units
    )
  }

  time_unit <- tolower(unit_parts[2])
  origin_text <- trimws(unit_parts[3])
  origin_text <- sub("\\s+UTC$", "", origin_text, ignore.case = TRUE)
  origin_text <- gsub("T", " ", origin_text, fixed = TRUE)

  origin_time <- as.POSIXct(
    origin_text,
    format = "%Y-%m-%d %H:%M:%OS",
    tz = "UTC"
  )

  if (is.na(origin_time)) {
    origin_time <- as.POSIXct(
      origin_text,
      format = "%Y-%m-%d",
      tz = "UTC"
    )
  }

  if (is.na(origin_time)) {
    stop(
      "Could not parse the time origin in ",
      path,
      ": ",
      origin_text
    )
  }

  seconds_per_unit <- switch(
    sub("s$", "", time_unit),
    second = 1,
    minute = 60,
    hour = 3600,
    day = 86400,
    stop(
      "Unsupported time unit in ",
      path,
      ": ",
      time_unit
    )
  )

  as.POSIXct(
    as.numeric(origin_time) +
      as.numeric(time_values) * seconds_per_unit,
    origin = "1970-01-01",
    tz = "UTC"
  )
}


assert_contiguous <- function(indices, description) {

  if (length(indices) == 0L) {
    return(invisible(NULL))
  }

  expected <- seq.int(
    min(indices),
    max(indices)
  )

  if (!identical(as.integer(indices), as.integer(expected))) {
    stop(description, " indices are not contiguous.")
  }

  invisible(NULL)
}


read_variable_block <- function(
    nc,
    variable,
    grid_info,
    latitude_indices,
    time_indices,
    path
) {

  require_variable(
    nc,
    variable,
    path
  )

  assert_contiguous(
    latitude_indices,
    "Latitude"
  )

  assert_contiguous(
    time_indices,
    "Time"
  )

  variable_dimensions <- nc$var[[variable]]$dim

  dimension_names <- vapply(
    variable_dimensions,
    function(x) x$name,
    character(1)
  )

  longitude_position <- match(
    grid_info$longitude_name,
    dimension_names
  )

  latitude_position <- match(
    grid_info$latitude_name,
    dimension_names
  )

  time_position <- match(
    grid_info$time_name,
    dimension_names
  )

  if (
    anyNA(
      c(
        longitude_position,
        latitude_position,
        time_position
      )
    )
  ) {
    stop(
      "Variable '",
      variable,
      "' in ",
      path,
      " does not use the expected longitude, latitude and time dimensions. ",
      "Its dimensions are: ",
      paste(dimension_names, collapse = ", ")
    )
  }

  dimension_lengths <- as.integer(
    vapply(
      variable_dimensions,
      function(x) x$len,
      numeric(1)
    )
  )

  start <- rep(1L, length(dimension_names))
  count <- dimension_lengths

  start[latitude_position] <- min(latitude_indices)
  count[latitude_position] <- length(latitude_indices)

  start[time_position] <- min(time_indices)
  count[time_position] <- length(time_indices)

  values <- ncvar_get(
    nc,
    variable,
    start = start,
    count = count,
    collapse_degen = FALSE
  )

  principal_positions <- c(
    longitude_position,
    latitude_position,
    time_position
  )

  extra_positions <- setdiff(
    seq_along(dimension_names),
    principal_positions
  )

  values <- aperm(
    values,
    c(
      principal_positions,
      extra_positions
    )
  )

  nlon <- length(grid_info$longitude)
  nlat <- length(latitude_indices)
  ntime <- length(time_indices)
  ncells <- nlon * nlat

  if (
    !identical(
      as.integer(dim(values)[1:3]),
      as.integer(c(nlon, nlat, ntime))
    )
  ) {
    stop(
      "Unexpected array dimensions when reading '",
      variable,
      "' from ",
      path,
      "."
    )
  }

  # ERA5 files can occasionally include an additional expver dimension.
  # Collapse any such dimensions by taking the first non-missing value.
  extra_size <- if (length(extra_positions) == 0L) {
    1L
  } else {
    as.integer(
      prod(dim(values)[-(1:3)])
    )
  }

  dim(values) <- c(
    ncells * ntime,
    extra_size
  )

  merged_values <- values[, 1L]

  if (extra_size > 1L) {
    for (icolumn in 2:extra_size) {
      replace_indices <- which(
        is.na(merged_values) &
          !is.na(values[, icolumn])
      )

      if (length(replace_indices) > 0L) {
        merged_values[replace_indices] <- values[
          replace_indices,
          icolumn
        ]
      }
    }
  }

  dim(merged_values) <- c(
    ncells,
    ntime
  )

  merged_values
}


grids_are_equal <- function(first, second) {

  isTRUE(
    all.equal(
      as.numeric(first$longitude),
      as.numeric(second$longitude),
      tolerance = 1e-10,
      check.attributes = FALSE
    )
  ) &&
    isTRUE(
      all.equal(
        as.numeric(first$latitude),
        as.numeric(second$latitude),
        tolerance = 1e-10,
        check.attributes = FALSE
      )
    )
}


times_are_equal <- function(first, second) {

  length(first) == length(second) &&
    all(
      abs(
        as.numeric(first) -
          as.numeric(second)
      ) < 1
    )
}


process_year <- function(current_year, iyear, number_of_years) {

  message(
    sprintf(
      "Processing %d (%d/%d)",
      current_year,
      iyear,
      number_of_years
    )
  )

  temperature_file <- file.path(
    input_dir,
    sprintf(
      temperature_pattern,
      current_year
    )
  )

  precipitation_file <- file.path(
    input_dir,
    sprintf(
      precipitation_pattern,
      current_year
    )
  )

  missing_files <- c(
    temperature_file,
    precipitation_file
  )[
    !file.exists(
      c(
        temperature_file,
        precipitation_file
      )
    )
  ]

  if (length(missing_files) > 0L) {
    message(
      "  Skipping year; missing file(s):\n    ",
      paste(missing_files, collapse = "\n    ")
    )

    return(NULL)
  }

  temperature_nc <- nc_open(temperature_file)

  on.exit(
    nc_close(temperature_nc),
    add = TRUE
  )

  precipitation_nc <- nc_open(precipitation_file)

  on.exit(
    nc_close(precipitation_nc),
    add = TRUE
  )

  require_variable(
    temperature_nc,
    "t2m",
    temperature_file
  )

  has_d2m <- "d2m" %in% names(temperature_nc$var)

  if (!has_d2m) {
    message(
      "  Variable 'd2m' is absent; continuing with temperature and ",
      "precipitation only."
    )
  }

  require_variable(
    precipitation_nc,
    "tp",
    precipitation_file
  )

  temperature_grid <- get_grid_info(
    temperature_nc,
    temperature_file
  )

  precipitation_grid <- get_grid_info(
    precipitation_nc,
    precipitation_file
  )

  if (
    !grids_are_equal(
      temperature_grid,
      precipitation_grid
    )
  ) {
    stop(
      "Temperature and precipitation grids differ for ",
      current_year,
      "."
    )
  }

  temperature_times <- parse_netcdf_time(
    temperature_nc,
    temperature_grid$time_name,
    temperature_file
  )

  precipitation_times <- parse_netcdf_time(
    precipitation_nc,
    precipitation_grid$time_name,
    precipitation_file
  )

  if (
    !times_are_equal(
      temperature_times,
      precipitation_times
    )
  ) {
    stop(
      "Temperature and precipitation times differ for ",
      current_year,
      "."
    )
  }

  time_dates <- as.Date(
    temperature_times,
    tz = "UTC"
  )

  time_year <- as.integer(
    format(
      temperature_times,
      "%Y",
      tz = "UTC"
    )
  )

  time_month <- as.integer(
    format(
      temperature_times,
      "%m",
      tz = "UTC"
    )
  )

  latitude_indices <- which(
    abs(temperature_grid$latitude) <= 30 + 1e-10
  )

  if (length(latitude_indices) == 0L) {
    stop(
      "No grid cells occur between 30 degrees south and 30 degrees north in ",
      temperature_file,
      "."
    )
  }

  assert_contiguous(
    latitude_indices,
    "Tropical latitude"
  )

  longitudes <- temperature_grid$longitude
  latitudes <- temperature_grid$latitude[latitude_indices]

  nlon <- length(longitudes)
  nlat <- length(latitudes)
  ncells <- nlon * nlat

  # Coordinates match the linearised longitude x latitude arrays.
  cell_longitude <- rep(
    longitudes,
    times = nlat
  )

  cell_latitude <- rep(
    latitudes,
    each = nlon
  )

  monthly_results <- vector(
    "list",
    12L
  )

  for (imonth in 1:12) {

    time_indices <- which(
      time_year == current_year &
        time_month == imonth
    )

    if (length(time_indices) == 0L) {
      message(
        sprintf(
          "  Month %02d is missing",
          imonth
        )
      )

      next
    }

    message(
      sprintf(
        "  Month %02d",
        imonth
      )
    )

    assert_contiguous(
      time_indices,
      sprintf("Time indices for %04d-%02d", current_year, imonth)
    )

    month_dates <- time_dates[time_indices]

    # -----------------------------------------------------------------------
    # Hourly precipitation accumulation (m) to monthly total (mm)
    # -----------------------------------------------------------------------

    total_precipitation <- read_variable_block(
      nc = precipitation_nc,
      variable = "tp",
      grid_info = precipitation_grid,
      latitude_indices = latitude_indices,
      time_indices = time_indices,
      path = precipitation_file
    )

    monthly_precipitation <- rowSums2(
      total_precipitation,
      na.rm = FALSE
    ) * 1000

    rm(total_precipitation)

    # -----------------------------------------------------------------------
    # Hourly air temperature
    # -----------------------------------------------------------------------

    air_temperature <- read_variable_block(
      nc = temperature_nc,
      variable = "t2m",
      grid_info = temperature_grid,
      latitude_indices = latitude_indices,
      time_indices = time_indices,
      path = temperature_file
    )

    days <- unique(month_dates)
    ndays <- length(days)

    monthly_tmin_sum <- numeric(ncells)
    monthly_tmax_sum <- numeric(ncells)
    monthly_tmean_sum <- numeric(ncells)

    for (iday in seq_along(days)) {

      day_columns <- which(
        month_dates == days[iday]
      )

      daily_tmin <- rowMins(
        air_temperature,
        cols = day_columns,
        na.rm = FALSE
      )

      daily_tmax <- rowMaxs(
        air_temperature,
        cols = day_columns,
        na.rm = FALSE
      )

      daily_tmean <- rowMeans2(
        air_temperature,
        cols = day_columns,
        na.rm = FALSE
      )

      monthly_tmin_sum <- monthly_tmin_sum + daily_tmin
      monthly_tmax_sum <- monthly_tmax_sum + daily_tmax
      monthly_tmean_sum <- monthly_tmean_sum + daily_tmean
    }

    monthly_tmin <- monthly_tmin_sum / ndays
    monthly_tmax <- monthly_tmax_sum / ndays
    monthly_tmean <- monthly_tmean_sum / ndays

    rm(
      daily_tmin,
      daily_tmax,
      daily_tmean,
      monthly_tmin_sum,
      monthly_tmax_sum,
      monthly_tmean_sum
    )

    # -----------------------------------------------------------------------
    # Final monthly table
    # -----------------------------------------------------------------------

    monthly_result <- data.table(
      year = current_year,
      month = imonth,
      lat = cell_latitude,
      lon = cell_longitude,
      pre = monthly_precipitation,
      tmin = monthly_tmin,
      tmax = monthly_tmax,
      tmp = monthly_tmean
    )

    if (has_d2m) {

      # ---------------------------------------------------------------------
      # Hourly dew-point temperature
      # ---------------------------------------------------------------------

      dewpoint_temperature <- read_variable_block(
        nc = temperature_nc,
        variable = "d2m",
        grid_info = temperature_grid,
        latitude_indices = latitude_indices,
        time_indices = time_indices,
        path = temperature_file
      )

      monthly_dewpoint <- rowMeans2(
        dewpoint_temperature,
        na.rm = TRUE
      )

      # ---------------------------------------------------------------------
      # Hourly vapour pressure deficit averaged over the month
      # ---------------------------------------------------------------------

      vpd_sum <- numeric(ncells)
      vpd_n <- integer(ncells)

      block_size <- 24L

      blocks <- split(
        seq_len(ncol(air_temperature)),
        ceiling(
          seq_len(ncol(air_temperature)) /
            block_size
        )
      )

      for (columns in blocks) {

        temperature_c <- (
          air_temperature[, columns, drop = FALSE] -
            273.15
        )

        dewpoint_c <- (
          dewpoint_temperature[, columns, drop = FALSE] -
            273.15
        )

        saturation_vapour_pressure <- 0.6108 * exp(
          17.27 * temperature_c /
            (temperature_c + 237.3)
        )

        actual_vapour_pressure <- 0.6108 * exp(
          17.27 * dewpoint_c /
            (dewpoint_c + 237.3)
        )

        vpd <- pmax(
          saturation_vapour_pressure -
            actual_vapour_pressure,
          0
        )

        valid <- !is.na(vpd)

        vpd_sum <- vpd_sum + rowSums2(
          vpd,
          na.rm = TRUE
        )

        vpd_n <- vpd_n + rowSums2(valid)

        rm(
          temperature_c,
          dewpoint_c,
          saturation_vapour_pressure,
          actual_vapour_pressure,
          vpd,
          valid
        )
      }

      monthly_vpd <- vpd_sum / vpd_n
      monthly_vpd[vpd_n == 0L] <- NaN

      monthly_result[
        ,
        `:=`(
          d2m = monthly_dewpoint,
          vpd = monthly_vpd
        )
      ]

      rm(
        dewpoint_temperature,
        monthly_dewpoint,
        monthly_vpd,
        vpd_sum,
        vpd_n,
        blocks
      )
    }

    monthly_results[[imonth]] <- monthly_result

    rm(
      air_temperature,
      monthly_precipitation,
      monthly_tmin,
      monthly_tmax,
      monthly_tmean,
      monthly_result
    )

    gc()
  }

  if (all(vapply(monthly_results, is.null, logical(1)))) {
    return(NULL)
  }

  rbindlist(
    monthly_results,
    use.names = TRUE,
    fill = TRUE
  )
}

# ---------------------------------------------------------------------------
# Process all years
# ---------------------------------------------------------------------------

year_results <- vector(
  "list",
  length(years)
)

for (iyear in seq_along(years)) {
  year_results[[iyear]] <- process_year(
    current_year = years[iyear],
    iyear = iyear,
    number_of_years = length(years)
  )

  gc()
}

# ---------------------------------------------------------------------------
# Combine all years and save
# ---------------------------------------------------------------------------

if (all(vapply(year_results, is.null, logical(1)))) {
  stop("No complete annual input pairs were found.")
}

df.ERA5 <- rbindlist(
  year_results,
  use.names = TRUE,
  fill = TRUE
)

saveRDS(
  df.ERA5,
  output_file,
  compress = FALSE
)

message(
  "Saved ",
  nrow(df.ERA5),
  " rows to ",
  output_file
)

# scp /Users/felicien/Documents/projects/Drying.CB/scripts/summarise.climate.ERA5_fast.R hpc:/data/gent/vo/000/gvo00074/felicien/R

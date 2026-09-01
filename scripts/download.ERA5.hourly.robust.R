#!/usr/bin/env Rscript

rm(list = ls())

# Download hourly ERA5 into two separate annual files:
#
#   ERA5_tropics_hourly_temperature_YYYY.nc
#       t2m
#
#   ERA5_tropics_hourly_precipitation_YYYY.nc
#       tp
#
# Complete historical years and the complete months of the current year are
# processed in parallel. Each worker creates its own Python/CDS
# client and downloads temperature followed by precipitation for one year.
# Temperature and precipitation are never merged.
#
# Worker selection, in decreasing priority:
#   1. First command-line argument, for example:
#        Rscript download_ERA5_hourly_parallel.R 4
#   2. SLURM_CPUS_PER_TASK when submitted through Slurm.
#   3. Two workers otherwise.
#
# Requirements:
#   module load CDO
#   R package: reticulate
#   Python package: cdsapi
#   CDS credentials in ~/.cdsapirc

suppressPackageStartupMessages({
  library(parallel)
  library(reticulate)
})

# ---------------------------------------------------------------------------
# User settings
# ---------------------------------------------------------------------------

current_year <- as.integer(format(Sys.Date(), "%Y"))
# years <- 1970:current_year
years <- 1940:1969

output_dir <- paste0(
  "/data/gent/vo/000/gvo00074/",
  "ED_common_data/met/Tropics"
)

times <- sprintf("%02d:00", 0:23)

# CDS area order: north, west, south, east.
area <- list(30, -180, -30, 180)

grid <- list(0.5, 0.5)

# The active selection downloads t2m only.
temperature_variables <- "2m_temperature"
temperature_short_names <- "t2m"

precipitation_variables <- "total_precipitation"
precipitation_short_names <- "tp"

max_attempts <- 4L
default_workers <- 2L

# ---------------------------------------------------------------------------
# Determine the number of parallel workers
# ---------------------------------------------------------------------------

get_worker_count <- function() {

  command_args <- commandArgs(trailingOnly = TRUE)

  if (length(command_args) >= 1L) {

    requested_workers <- suppressWarnings(
      as.integer(command_args[1])
    )

    if (
      is.na(requested_workers) ||
      requested_workers < 1L
    ) {
      stop(
        "The first command-line argument must be a positive integer: ",
        command_args[1]
      )
    }

    return(requested_workers)
  }

  slurm_workers <- Sys.getenv(
    "SLURM_CPUS_PER_TASK",
    unset = ""
  )

  if (nzchar(slurm_workers)) {

    requested_workers <- suppressWarnings(
      as.integer(slurm_workers)
    )

    if (
      !is.na(requested_workers) &&
      requested_workers >= 1L
    ) {
      return(requested_workers)
    }
  }

  default_workers
}

n_workers <- get_worker_count()
n_workers <- min(n_workers, length(years))

# ---------------------------------------------------------------------------
# Initial checks
# ---------------------------------------------------------------------------

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cdo <- unname(Sys.which("cdo"))

if (!nzchar(cdo)) {
  stop(
    "CDO was not found. Load it before starting R:\n",
    "  module load CDO"
  )
}

if (!py_module_available("cdsapi")) {
  stop(
    "Python module 'cdsapi' is unavailable in the Python ",
    "environment used by reticulate."
  )
}

# Every worker will explicitly use this same Python installation.
python_bin <- py_config()$python

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

cdo_values <- function(operator, path) {

  output <- system2(
    cdo,
    args = c(
      "-s",
      operator,
      shQuote(path)
    ),
    stdout = TRUE,
    stderr = FALSE
  )

  status <- attr(output, "status")

  if (!is.null(status) && status != 0L) {
    return(character())
  }

  output_text <- trimws(
    paste(output, collapse = " ")
  )

  if (!nzchar(output_text)) {
    return(character())
  }

  strsplit(output_text, "\\s+")[[1]]
}


months_to_download <- function(year) {

  current_year <- as.integer(format(Sys.Date(), "%Y"))
  current_month <- as.integer(format(Sys.Date(), "%m"))

  if (year < current_year) {
    return(1:12)
  }

  if (year == current_year) {

    # Download complete calendar months only. This avoids treating the
    # incomplete current month as though it were a complete annual record.
    if (current_month == 1L) {
      return(integer())
    }

    return(seq_len(current_month - 1L))
  }

  stop("Cannot download a future year: ", year)
}


expected_timesteps_for_months <- function(
    year,
    selected_months
) {

  if (length(selected_months) == 0L) {
    return(0L)
  }

  first_dates <- as.Date(
    sprintf(
      "%04d-%02d-01",
      year,
      selected_months
    )
  )

  first_dates_next_month <- as.Date(
    format(
      first_dates + 32,
      "%Y-%m-01"
    )
  )

  days_per_month <- as.integer(
    first_dates_next_month - first_dates
  )

  # One timestep for every requested hour of every selected day.
  as.integer(sum(days_per_month) * length(times))
}


netcdf_is_valid <- function(
    path,
    expected_variables,
    expected_timesteps
) {

  if (
    !file.exists(path) ||
    is.na(file.info(path)$size) ||
    file.info(path)$size <= 0
  ) {
    return(FALSE)
  }

  variable_names <- tryCatch(
    cdo_values("showname", path),
    error = function(e) character()
  )

  if (!all(expected_variables %in% variable_names)) {
    return(FALSE)
  }

  ntime <- tryCatch(
    as.integer(cdo_values("ntime", path)[1]),
    error = function(e) NA_integer_
  )

  if (
    is.na(ntime) ||
    ntime != expected_timesteps
  ) {
    return(FALSE)
  }

  TRUE
}


retrieve_with_retry <- function(
    client,
    request,
    target,
    expected_variables,
    expected_timesteps,
    year,
    data_type
) {

  if (
    netcdf_is_valid(
      target,
      expected_variables,
      expected_timesteps
    )
  ) {
    message(
      "[", year, "][", data_type, "] Reusing completed file: ",
      target
    )

    return(invisible(target))
  }

  partial_file <- paste0(target, ".part")

  for (attempt in seq_len(max_attempts)) {

    if (file.exists(partial_file)) {
      unlink(partial_file)
    }

    message(
      "[", year, "][", data_type, "] Downloading ",
      basename(target),
      " (attempt ",
      attempt,
      "/",
      max_attempts,
      ")"
    )

    download_error <- tryCatch(
      {
        client$retrieve(
          "reanalysis-era5-single-levels",
          request,
          partial_file
        )

        NULL
      },
      error = identity
    )

    if (
      inherits(download_error, "error") &&
      grepl(
        "403|cost limits exceeded|request.*too large",
        conditionMessage(download_error),
        ignore.case = TRUE
      )
    ) {
      stop(
        "CDS rejected the annual ",
        data_type,
        " request for ",
        year,
        " as too large:\n",
        conditionMessage(download_error),
        "\n\nParallel processing cannot reduce the size of an ",
        "individual request. This request must be split into smaller periods."
      )
    }

    if (
      is.null(download_error) &&
      netcdf_is_valid(
        partial_file,
        expected_variables,
        expected_timesteps
      )
    ) {

      # Remove an old invalid target only after the replacement has been
      # completely downloaded and validated.
      if (file.exists(target)) {
        unlink(target)
      }

      if (!file.rename(partial_file, target)) {
        stop(
          "Could not move completed download to: ",
          target
        )
      }

      message(
        "[", year, "][", data_type, "] Completed ",
        basename(target)
      )

      return(invisible(target))
    }

    if (attempt == max_attempts) {

      if (inherits(download_error, "error")) {
        stop(
          "CDS download failed for ",
          year,
          " (",
          data_type,
          "): ",
          conditionMessage(download_error)
        )
      }

      stop(
        "CDS returned an invalid or incomplete file for ",
        year,
        " (",
        data_type,
        "): ",
        partial_file
      )
    }

    wait_seconds <- min(
      60,
      5 * 2^(attempt - 1)
    )

    message(
      "[", year, "][", data_type, "] Download incomplete; ",
      "retrying in ",
      wait_seconds,
      " seconds."
    )

    Sys.sleep(wait_seconds)
  }
}


download_year <- function(year) {

  selected_months <- months_to_download(year)

  if (length(selected_months) == 0L) {
    message(
      "[", year, "] No complete calendar months are available yet."
    )

    return(invisible(year))
  }

  expected_timesteps_year <- expected_timesteps_for_months(
    year,
    selected_months
  )

  message(
    "\n[", year, "] Starting hourly ERA5 downloads for month",
    if (length(selected_months) == 1L) " " else "s ",
    paste(
      sprintf("%02d", selected_months),
      collapse = ", "
    ),
    " (expected timesteps: ",
    expected_timesteps_year,
    ")"
  )

  temperature_file <- file.path(
    output_dir,
    sprintf(
      "ERA5_tropics_hourly_temperature_%04d.nc",
      year
    )
  )

  precipitation_file <- file.path(
    output_dir,
    sprintf(
      "ERA5_tropics_hourly_precipitation_%04d.nc",
      year
    )
  )

  common_request <- list(
    product_type = "reanalysis",
    data_format = "netcdf",
    download_format = "unarchived",
    year = as.character(year),
    month = as.list(sprintf("%02d", selected_months)),
    day = as.list(sprintf("%02d", 1:31)),
    time = as.list(times),
    area = area,
    grid = grid
  )

  # Hourly 2-m air temperature.
  temperature_request <- c(
    common_request,
    list(
      variable = as.list(temperature_variables)
    )
  )

  retrieve_with_retry(
    client = cds_client,
    request = temperature_request,
    target = temperature_file,
    expected_variables = temperature_short_names,
    expected_timesteps = expected_timesteps_year,
    year = year,
    data_type = "temperature"
  )

  # Hourly total precipitation.
  precipitation_request <- c(
    common_request,
    list(
      variable = as.list(precipitation_variables)
    )
  )

  retrieve_with_retry(
    client = cds_client,
    request = precipitation_request,
    target = precipitation_file,
    expected_variables = precipitation_short_names,
    expected_timesteps = expected_timesteps_year,
    year = year,
    data_type = "precipitation"
  )

  message("[", year, "] All downloads completed")

  invisible(year)
}

# ---------------------------------------------------------------------------
# Create the parallel cluster
# ---------------------------------------------------------------------------

message(
  "Using ",
  n_workers,
  " parallel worker",
  if (n_workers == 1L) "" else "s",
  "."
)

cluster <- makePSOCKcluster(
  n_workers,
  outfile = ""
)

# Send the selected Python executable to every worker before initializing
# reticulate there.
clusterExport(
  cluster,
  varlist = "python_bin",
  envir = .GlobalEnv
)

# A reticulate Python object cannot safely be shared between R processes.
# Consequently, each worker creates and retains its own CDS API client.
clusterEvalQ(
  cluster,
  {
    reticulate::use_python(
      python_bin,
      required = TRUE
    )

    cdsapi_worker <- reticulate::import(
      "cdsapi",
      delay_load = FALSE
    )

    cds_client <- cdsapi_worker$Client()

    NULL
  }
)

clusterExport(
  cluster,
  varlist = c(
    "download_year",
    "retrieve_with_retry",
    "netcdf_is_valid",
    "cdo_values",
    "months_to_download",
    "expected_timesteps_for_months",
    "output_dir",
    "times",
    "area",
    "grid",
    "temperature_variables",
    "temperature_short_names",
    "precipitation_variables",
    "precipitation_short_names",
    "max_attempts",
    "cdo"
  ),
  envir = .GlobalEnv
)

# ---------------------------------------------------------------------------
# Download years in parallel
# ---------------------------------------------------------------------------

results <- tryCatch(
  {
    parLapplyLB(
      cluster,
      years,
      function(year) {

        tryCatch(
          {
            download_year(year)

            list(
              year = year,
              success = TRUE,
              error = NA_character_
            )
          },
          error = function(e) {
            list(
              year = year,
              success = FALSE,
              error = conditionMessage(e)
            )
          }
        )
      }
    )
  },
  finally = {
    stopCluster(cluster)
  }
)

# ---------------------------------------------------------------------------
# Final report
# ---------------------------------------------------------------------------

success <- vapply(
  results,
  function(result) isTRUE(result$success),
  logical(1)
)

if (any(!success)) {

  failed_messages <- vapply(
    results[!success],
    function(result) {
      sprintf(
        "%04d: %s",
        result$year,
        result$error
      )
    },
    character(1)
  )

  stop(
    "Some annual downloads failed:\n",
    paste(
      failed_messages,
      collapse = "\n"
    )
  )
}

message(
  "\nAll ERA5 downloads completed successfully for ",
  min(years),
  "-",
  max(years),
  ". The current-year files contain complete calendar months only."
)

# scp /Users/felicien/Documents/projects/Drying.CB/scripts/download.ERA5.hourly.robust.R hpc:/kyukon/data/gent/vo/000/gvo00074/felicien/R/

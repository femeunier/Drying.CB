rm(list = ls())
gc()

library(dplyr)
library(lubridate)
library(terra)
library(Drying.CB)

# -------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------

input_dir <- "/data/gent/vo/000/gvo00074/felicien/R/outputs/all.climate"

output_dir <- "/data/gent/vo/000/gvo00074/felicien/R/outputs/Drying.CB"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

summary_dir <- file.path(output_dir, "timeseries_rainfor")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

terra_tmp <- "/data/gent/vo/000/gvo00074/felicien/R/terra_tmp"
dir.create(terra_tmp, recursive = TRUE, showWarnings = FALSE)

stopifnot(
  dir.exists(input_dir),
  dir.exists(output_dir),
  dir.exists(summary_dir),
  dir.exists(terra_tmp)
)

# Set this below the memory allocated to each Slurm task.
# For example, use memmax = 24 with --mem=32G.
terra::terraOptions(
  tempdir  = terra_tmp,
  memfrac  = 0.25,
  memmax   = 24,
  todisk   = TRUE,
  progress = 1
)

baseline_start <- as.Date("1985-01-01")
baseline_end   <- as.Date("2014-12-31")

suffix <- "rainfor"

datasets_to_keep <- c(
  "3IMERG", "Berk", "CAMS", "chirps", "chirpsv3", "CRU",
  "ERA5", "GLDAS", "GPCC", "MSWEP", "NCEP", "GSMaP",
  "CPC", "JRA3Q"
)

datasets_to_exclude <- c(
  "CHELSA", "CRUJRA", "MERRA2", "W5E5", "GSWP"
)

vars_to_keep <- c(
  "tas", "tasmin", "tasmax", "pre"
)

# -------------------------------------------------------------------------
# Identify input files
# -------------------------------------------------------------------------

files <- list.files(
  input_dir,
  pattern = "\\.tif$",
  full.names = TRUE
)

file_info <- lapply(basename(files), function(x) {
  parts <- strsplit(x, "_", fixed = TRUE)[[1]]

  data.frame(
    product = if (length(parts) >= 1L) parts[1] else NA_character_,
    variable = if (length(parts) >= 2L) parts[2] else NA_character_,
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows()

keep <- (
  file_info$product %in% datasets_to_keep &
    !(file_info$product %in% datasets_to_exclude) &
    file_info$variable %in% vars_to_keep
)

files <- files[keep]
file_info <- file_info[keep, , drop = FALSE]

ord <- order(file_info$product, file_info$variable)
files <- files[ord]
file_info <- file_info[ord, , drop = FALSE]

if (length(files) == 0L) {
  stop("No matching input files were found.")
}

message("Number of input files: ", length(files))

# Allow the number of jobs to be queried before submitting the array
if ("--count" %in% commandArgs(trailingOnly = TRUE)) {
  cat(length(files), "\n")
  quit(save = "no", status = 0)
}

# Under Slurm, each array task processes one file.
# Outside Slurm, files are processed sequentially.
task_id <- suppressWarnings(
  as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA_character_))
)

if (is.na(task_id)) {
  file_indices <- seq_along(files)
} else {
  if (task_id < 1L || task_id > length(files)) {
    stop(
      "SLURM_ARRAY_TASK_ID = ", task_id,
      " but only ", length(files), " files were found."
    )
  }

  file_indices <- task_id
}

# -------------------------------------------------------------------------
# Load mask
# -------------------------------------------------------------------------

# Mask <- read_sf("./data/Rainforests.shp")
# Mask <- vect(st_as_sfc(
#   st_bbox(c(
#     xmin = -20,
#     ymin = -15,
#     xmax = 55,
#     ymax = 15),
#   crs = 4326)
# ))
Mask <- rast("./data/shapefiles/LC_evergreen_01deg_mask.tif")

if (nlyr(Mask) != 1L) {
  Mask <- Mask[[1]]
}

# -------------------------------------------------------------------------
# Output-writing helper
# -------------------------------------------------------------------------

write_anomaly_raster <- function(x, product, variable, label, dates = NULL) {

  if (!is.null(dates)) {
    time(x) <- dates
  }

  output_file <- file.path(
    output_dir,
    paste0(
      product, "_", variable, "_",
      label, "_", suffix, ".tif"
    )
  )

  writeRaster(
    x,
    output_file,
    overwrite = TRUE,
    datatype = "FLT4S",
    gdal = c(
      "COMPRESS=DEFLATE",
      "PREDICTOR=3",
      "BIGTIFF=YES",
      "TILED=YES"
    )
  )

  invisible(output_file)
}

# -------------------------------------------------------------------------
# Process one file
# -------------------------------------------------------------------------

process_file <- function(ifile) {

  cfile <- files[ifile]
  cproduct <- file_info$product[ifile]
  cvar <- file_info$variable[ifile]

  message(
    "\nProcessing ", cproduct, "-", cvar,
    " [", ifile, "/", length(files), "]"
  )
  message("Input: ", cfile)

  temporary_files <- character()

  on.exit(
    {
      rm(list = intersect(
        c(
          "cdata", "cdata_crop", "mask_aligned",
          "cdata_msk", "anomalies"
        ),
        ls()
      ))

      gc()

      existing_temp <- temporary_files[file.exists(temporary_files)]
      if (length(existing_temp) > 0L) {
        unlink(existing_temp, force = TRUE)
      }

      try(
        terra::tmpFiles(current = TRUE, remove = TRUE),
        silent = TRUE
      )

      gc()
    },
    add = TRUE
  )

  cdata <- rast(cfile)

  input_dates <- time(cdata)

  if (length(input_dates) != nlyr(cdata)) {
    stop("Invalid or missing time dimension in ", cfile)
  }

  # Crop first so that the subsequent mask and anomaly calculations
  # do not retain cells outside the mask extent.
  cdata_crop <- crop(
    cdata,
    ext(Mask),
    snap = "near"
  )

  # Use a single-layer target when aligning the mask.
  mask_aligned <- resample(
    Mask,
    cdata_crop[[1]],
    method = "near"
  )

  # Materialize the masked raster on disk. This is important: subsequent
  # calculations then refer to a file-backed raster rather than to a long
  # chain of unevaluated crop/resample/mask operations.
  masked_file <- tempfile(
    pattern = paste0(cproduct, "_", cvar, "_masked_"),
    tmpdir = terra_tmp,
    fileext = ".tif"
  )

  temporary_files <- c(temporary_files, masked_file)

  cdata_msk <- mask(
    cdata_crop,
    mask_aligned,
    maskvalues = 0,
    updatevalue = NA,
    filename = masked_file,
    overwrite = TRUE,
    wopt = list(
      datatype = "FLT4S",
      gdal = c(
        "COMPRESS=DEFLATE",
        "PREDICTOR=3",
        "BIGTIFF=YES",
        "TILED=YES"
      )
    )
  )

  time(cdata_msk) <- input_dates

  # Release the original global raster and crop before calculating
  # the anomalies.
  rm(cdata, cdata_crop, mask_aligned)
  gc()

  # -----------------------------------------------------------------------
  # Spatial means
  # -----------------------------------------------------------------------

  ts <- global(
    cdata_msk,
    fun = "mean",
    na.rm = TRUE
  )

  values_ts <- ts[[1]]

  temperatures_in_kelvin <- (
    cvar %in% c("tas", "tasmin", "tasmax") &&
      mean(values_ts, na.rm = TRUE) > 200
  )

  if (temperatures_in_kelvin) {

    message("Converting temperature from K to degrees C")

    values_ts <- values_ts - 273.15

    converted_file <- tempfile(
      pattern = paste0(cproduct, "_", cvar, "_celsius_"),
      tmpdir = terra_tmp,
      fileext = ".tif"
    )

    temporary_files <- c(temporary_files, converted_file)

    cdata_msk <- writeRaster(
      cdata_msk - 273.15,
      converted_file,
      overwrite = TRUE,
      datatype = "FLT4S",
      gdal = c(
        "COMPRESS=DEFLATE",
        "PREDICTOR=3",
        "BIGTIFF=YES",
        "TILED=YES"
      )
    )

    time(cdata_msk) <- input_dates

    gc()

    # The converted raster no longer depends on masked_file.
    if (file.exists(masked_file)) {
      unlink(masked_file, force = TRUE)
    }
  }

  cdf2save <- data.frame(
    time = input_dates,
    value = values_ts
  ) %>%
    mutate(
      year = year(time),
      month = month(time),
      var = cvar,
      product = cproduct
    )

  summary_file <- file.path(
    summary_dir,
    paste0(
      tools::file_path_sans_ext(basename(cfile)),
      "_", suffix, ".RDS"
    )
  )

  saveRDS(cdf2save, summary_file)

  rm(ts, values_ts, cdf2save)
  gc()

  # -----------------------------------------------------------------------
  # Anomalies
  # -----------------------------------------------------------------------

  message("Calculating anomalies")

  anomalies <- anomalies_spatraster_roll(
    input = cdata_msk,
    baseline_start = baseline_start,
    baseline_end = baseline_end,
    detrend = FALSE,
    SD_threshold = 0.001,
    Z_anom_threshold = 7
  )

  # The input raster is no longer needed after the function has completed.
  rm(cdata_msk)
  gc()

  # Write one result and immediately remove it from the returned list.
  write_anomaly_raster(
    anomalies$trend,
    cproduct,
    cvar,
    "trends"
  )
  anomalies$trend <- NULL
  gc()

  write_anomaly_raster(
    anomalies$anom,
    cproduct,
    cvar,
    "anomalies",
    dates = input_dates
  )
  anomalies$anom <- NULL
  gc()

  write_anomaly_raster(
    anomalies$z_anom,
    cproduct,
    cvar,
    "Zanomalies",
    dates = input_dates
  )
  anomalies$z_anom <- NULL
  gc()

  roll_dates <- as.Date(anomalies$roll_times)

  write_anomaly_raster(
    anomalies$roll_mean_input,
    cproduct,
    cvar,
    "Rollmeaninput",
    dates = roll_dates
  )
  anomalies$roll_mean_input <- NULL
  anomalies$roll_times <- NULL
  rm(roll_dates)
  gc()

  write_anomaly_raster(
    anomalies$trend_anom,
    cproduct,
    cvar,
    "trendsanomalies"
  )
  anomalies$trend_anom <- NULL
  gc()

  write_anomaly_raster(
    anomalies$trend_z_anom,
    cproduct,
    cvar,
    "trendsZanomalies"
  )
  anomalies$trend_z_anom <- NULL
  gc()

  rm(anomalies)
  gc()

  message("Finished ", cproduct, "-", cvar)

  invisible(summary_file)
}

# -------------------------------------------------------------------------
# Run
# -------------------------------------------------------------------------

for (ifile in file_indices) {
  process_file(ifile)

  try(
    terra::tmpFiles(current = TRUE, remove = TRUE),
    silent = TRUE
  )

  gc()
}

# scp /Users/felicien/Documents/projects/Drying.CB/scripts/Extract.allclimate.var.R hpc:/data/gent/vo/000/gvo00074/felicien/R/

rm(list = ls())
gc()

library(dplyr)
library(lubridate)
library(terra)
library(Drying.CB)

# Load the blockwise function when it is not already available from the package
# or the current R session.
if (!exists("anomalies_spatraster_roll_blockwise", mode = "function")) {
  blockwise_function_file <- "./anomalies_spatraster_roll_blockwise.R"

  if (!file.exists(blockwise_function_file)) {
    stop(
      "Could not find ", blockwise_function_file,
      ". Copy it next to this extraction script or install it in Drying.CB."
    )
  }

  source(blockwise_function_file)
}

stopifnot(
  exists("anomalies_spatraster_roll_blockwise", mode = "function")
)

# -------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------

input_dir <- "/data/gent/vo/000/gvo00074/felicien/R/outputs/all.climate"

output_dir <- "/data/gent/vo/000/gvo00074/felicien/R/outputs/Drying.CB"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

summary_dir <- file.path(output_dir, "timeseries_rainfor")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(
  dir.exists(input_dir),
  dir.exists(output_dir),
  dir.exists(summary_dir)
)

# Use node-local scratch for all intermediates. Writing Terra's implicit
# temporary rasters to the shared /data filesystem caused files to disappear
# between crop/resample and mask operations.
node_tmp <- Sys.getenv("TMPDIR", unset = "")
if (!nzchar(node_tmp)) node_tmp <- tempdir()

terra_tmp <- file.path(
  node_tmp,
  paste0(
    "DryingCB_",
    Sys.getenv("SLURM_JOB_ID", unset = "local"),
    "_task_",
    Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "0"),
    "_pid_",
    Sys.getpid()
  )
)

dir.create(terra_tmp, recursive = TRUE, showWarnings = FALSE)
stopifnot(dir.exists(terra_tmp))

message("Terra temporary directory: ", terra_tmp)

terra::terraOptions(
  tempdir = terra_tmp,
  memfrac = 0.25,
  memmax = 24,
  todisk = FALSE,
  progress = 1
)

baseline_start <- as.Date("1985-01-01")
baseline_end <- as.Date("2014-12-31")

suffix <- "rainfor"

# Set TRUE to resume safely after an interrupted run. A product is skipped only
# when its summary and all six expected raster outputs exist and can be opened
# with the expected number of layers.
skip_completed <- TRUE

datasets_to_keep <- c(
  "3IMERG", "Berk", "CAMS", "chirps", "chirpsv3", "CRU",
  "ERA5", "GLDAS", "GPCC", "MSWEP", "NCEP", "GSMaP",
  "CPC", "JRA3Q"
)

datasets_to_exclude <- c(
  "CHELSA", "CRUJRA", "MERRA2", "W5E5", "GSWP"
)

vars_to_keep <- c("tas", "tasmin", "tasmax", "pre")

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
}) |>
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

if (!length(files)) stop("No matching input files were found.")

message("Number of input files: ", length(files))

# Allow the number of jobs to be queried before submitting a Slurm array.
if ("--count" %in% commandArgs(trailingOnly = TRUE)) {
  cat(length(files), "\n")
  quit(save = "no", status = 0)
}

# Under Slurm, each array task processes one file. Outside Slurm, the selected
# files are processed sequentially in the current R process.
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

Mask <- rast("./data/shapefiles/LC_evergreen_01deg_mask.tif")
if (nlyr(Mask) != 1L) Mask <- Mask[[1]]

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

make_anomaly_files <- function(product, variable) {
  c(
    trend = file.path(
      output_dir,
      paste0(product, "_", variable, "_trends_", suffix, ".tif")
    ),
    anom = file.path(
      output_dir,
      paste0(product, "_", variable, "_anomalies_", suffix, ".tif")
    ),
    z_anom = file.path(
      output_dir,
      paste0(product, "_", variable, "_Zanomalies_", suffix, ".tif")
    ),
    roll_mean_input = file.path(
      output_dir,
      paste0(product, "_", variable, "_Rollmeaninput_", suffix, ".tif")
    ),
    trend_anom = file.path(
      output_dir,
      paste0(product, "_", variable, "_trendsanomalies_", suffix, ".tif")
    ),
    trend_z_anom = file.path(
      output_dir,
      paste0(product, "_", variable, "_trendsZanomalies_", suffix, ".tif")
    )
  )
}

raster_has_nlayers <- function(filename, expected_nlayers) {
  if (!file.exists(filename)) return(FALSE)

  tryCatch(
    {
      x <- rast(filename)
      identical(as.integer(nlyr(x)), as.integer(expected_nlayers))
    },
    error = function(e) FALSE
  )
}

product_is_complete <- function(
    completion_file,
    summary_file,
    anomaly_files,
    input_nlayers) {

  # Raster headers and layer counts can look valid even when a job stopped
  # halfway through its spatial blocks. The marker is written only after the
  # final validation at the end of process_file().
  if (!file.exists(completion_file)) return(FALSE)
  if (!file.exists(summary_file)) return(FALSE)

  expected_nlayers <- c(
    trend = 2L,
    anom = input_nlayers,
    z_anom = input_nlayers,
    roll_mean_input = input_nlayers,
    trend_anom = 2L,
    trend_z_anom = 2L
  )

  checks <- vapply(
    names(anomaly_files),
    function(key) {
      raster_has_nlayers(
        anomaly_files[[key]],
        expected_nlayers[[key]]
      )
    },
    logical(1)
  )

  all(checks)
}

raster_write_options <- list(
  datatype = "FLT4S",
  gdal = c(
    "COMPRESS=DEFLATE",
    "PREDICTOR=3",
    "BIGTIFF=YES",
    "TILED=YES"
  )
)

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

  file_stem <- tools::file_path_sans_ext(basename(cfile))

  summary_file <- file.path(
    summary_dir,
    paste0(file_stem, "_", suffix, ".RDS")
  )

  anomaly_files <- make_anomaly_files(cproduct, cvar)
  completion_file <- file.path(
    output_dir,
    paste0(cproduct, "_", cvar, "_", suffix, ".complete")
  )
  temporary_files <- character()

  # Only files explicitly registered in temporary_files are removed. Never use
  # terra::tmpFiles(remove=TRUE), as a live SpatRaster may still reference one.
  on.exit(
    {
      local_objects <- intersect(
        c(
          "cdata", "cdata_crop", "mask_aligned", "cdata_msk",
          "cdata_celsius", "anomalies", "ts", "cdf2save"
        ),
        ls(envir = environment())
      )

      if (length(local_objects)) {
        rm(list = local_objects, envir = environment())
      }

      gc(verbose = FALSE)

      existing_temp <- temporary_files[file.exists(temporary_files)]
      if (length(existing_temp)) unlink(existing_temp, force = TRUE)

      gc(verbose = FALSE)
    },
    add = TRUE
  )

  cdata <- rast(cfile)
  input_dates <- time(cdata)
  input_nlayers <- nlyr(cdata)

  if (is.null(input_dates) || length(input_dates) != input_nlayers) {
    stop("Invalid or missing time dimension in ", cfile)
  }

  if (
    isTRUE(skip_completed) &&
    product_is_complete(
      completion_file,
      summary_file,
      anomaly_files,
      input_nlayers
    )
  ) {
    message("Already complete and readable; skipping ", cproduct, "-", cvar)
    return(invisible(summary_file))
  }

  # A new attempt invalidates any stale marker before output files are opened.
  if (file.exists(completion_file)) unlink(completion_file, force = TRUE)

  crop_file <- file.path(terra_tmp, paste0(file_stem, "_crop.tif"))
  aligned_mask_file <- file.path(
    terra_tmp,
    paste0(file_stem, "_mask_aligned.tif")
  )
  masked_file <- file.path(terra_tmp, paste0(file_stem, "_masked.tif"))

  temporary_files <- c(
    temporary_files,
    crop_file,
    aligned_mask_file,
    masked_file
  )

  # Every potentially large spatial intermediate is explicitly materialized on
  # node-local disk. No implicit shared-filesystem spat_*.tif is needed.
  cdata_crop <- crop(
    cdata,
    ext(Mask),
    snap = "near",
    filename = crop_file,
    overwrite = TRUE,
    wopt = raster_write_options
  )

  mask_aligned <- resample(
    Mask,
    cdata_crop[[1]],
    method = "near",
    filename = aligned_mask_file,
    overwrite = TRUE,
    wopt = list(
      datatype = "INT1U",
      gdal = c(
        "COMPRESS=DEFLATE",
        "BIGTIFF=YES",
        "TILED=YES"
      )
    )
  )

  cdata_msk <- mask(
    cdata_crop,
    mask_aligned,
    maskvalues = 0,
    updatevalue = NA,
    filename = masked_file,
    overwrite = TRUE,
    wopt = raster_write_options
  )

  time(cdata_msk) <- input_dates

  rm(cdata, cdata_crop, mask_aligned)
  gc(verbose = FALSE)

  # The completed masked raster no longer depends on the crop or aligned mask.
  unlink(c(crop_file, aligned_mask_file), force = TRUE)

  # -----------------------------------------------------------------------
  # Spatial means and temperature-unit conversion
  # -----------------------------------------------------------------------

  ts <- global(cdata_msk, fun = "mean", na.rm = TRUE)
  values_ts <- ts[[1]]

  temperatures_in_kelvin <- (
    cvar %in% c("tas", "tasmin", "tasmax") &&
      mean(values_ts, na.rm = TRUE) > 200
  )

  if (temperatures_in_kelvin) {
    message("Converting temperature from K to degrees C")

    values_ts <- values_ts - 273.15

    converted_file <- file.path(
      terra_tmp,
      paste0(file_stem, "_celsius.tif")
    )
    temporary_files <- c(temporary_files, converted_file)

    cdata_celsius <- writeRaster(
      cdata_msk - 273.15,
      converted_file,
      overwrite = TRUE,
      datatype = "FLT4S",
      gdal = raster_write_options$gdal
    )
    time(cdata_celsius) <- input_dates

    rm(cdata_msk)
    gc(verbose = FALSE)
    unlink(masked_file, force = TRUE)

    cdata_msk <- cdata_celsius
    rm(cdata_celsius)
  }

  cdf2save <- data.frame(
    time = input_dates,
    value = values_ts
  ) |>
    mutate(
      year = year(time),
      month = month(time),
      var = cvar,
      product = cproduct
    )

  saveRDS(cdf2save, summary_file)

  rm(ts, values_ts, cdf2save)
  gc(verbose = FALSE)

  # -----------------------------------------------------------------------
  # Blockwise anomalies: the function writes directly to anomaly_files.
  # -----------------------------------------------------------------------

  message("Calculating anomalies")

  anomalies <- anomalies_spatraster_roll_blockwise(
    input = cdata_msk,
    baseline_start = baseline_start,
    baseline_end = baseline_end,
    detrend = FALSE,
    roll_window = 12L,
    roll_align = "right",
    SD_threshold = 0.001,
    Z_anom_threshold = 7,
    outputs = names(anomaly_files),
    filenames = anomaly_files,
    overwrite = TRUE,
    # Larger blocks substantially reduce R and GeoTIFF write overhead while
    # remaining well below the 24 GB process limit.
    block_mem_mb = 2048,
    finalize_retries = 60L,
    finalize_wait_seconds = 2
  )

  # Do not call writeRaster() here: every returned SpatRaster is already backed
  # by its final target file.
  rm(anomalies, cdata_msk)
  gc(verbose = FALSE)

  # Validate the final files before declaring the product complete.
  expected_nlayers <- c(
    trend = 2L,
    anom = input_nlayers,
    z_anom = input_nlayers,
    roll_mean_input = input_nlayers,
    trend_anom = 2L,
    trend_z_anom = 2L
  )

  final_checks <- vapply(
    names(anomaly_files),
    function(key) {
      raster_has_nlayers(anomaly_files[[key]], expected_nlayers[[key]])
    },
    logical(1)
  )

  if (!file.exists(summary_file) || !all(final_checks)) {
    stop("One or more final outputs failed validation for ", cproduct, "-", cvar)
  }

  writeLines(
    c(
      paste0("product=", cproduct),
      paste0("variable=", cvar),
      paste0("completed=", format(Sys.time(), tz = "UTC", usetz = TRUE)),
      paste0("input=", cfile)
    ),
    completion_file
  )

  message("Finished ", cproduct, "-", cvar)
  invisible(summary_file)
}

# -------------------------------------------------------------------------
# Run
# -------------------------------------------------------------------------

for (ifile in file_indices) {
  process_file(ifile)
  gc(verbose = FALSE)
}


# scp /Users/felicien/Documents/projects/Drying.CB/scripts/Extract.allclimate.var.block.R hpc:/data/gent/vo/000/gvo00074/felicien/R/

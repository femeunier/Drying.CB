rm(list = ls())

library(terra)
library(data.table)

# -------------------------------------------------------------------------
# Settings
# -------------------------------------------------------------------------

input_dir  <- "./outputs"
output_dir <- "./outputs/all.climate"

overwrite <- TRUE

datasets_to_keep <- paste(
  c(
    "3IMERG", "Berk", "CAMS", "chirps", "chirpsv3", "CRU",
    "ERA5", "GLDAS", "GPCC", "MSWEP", "NCEP", "GSMaP",
    "CPC", "JRA3Q"
  ),
  collapse = "|"
)

datasets_to_exclude <- paste(
  c("CHELSA", "CRUJRA", "MERRA2", "W5E5", "GSWP"),
  collapse = "|"
)

# Maximum number of target cells passed to RANN at once
nn_chunk_size <- 500000L

# Maximum accepted nearest-neighbour distance relative to cell diagonal.
# This prevents missing CHIRPS cells from being filled from distant cells.
nn_distance_factor <- 0.55

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

terraOptions(
  memfrac = 0.25,
  todisk = TRUE,
  progress = 1
)

# -------------------------------------------------------------------------
# Helper: standardise variable names
# -------------------------------------------------------------------------

rename_to_standard <- function(
    x,
    target,
    aliases,
    filename
) {

  found <- intersect(
    c(target, aliases),
    names(x)
  )

  if (length(found) == 0L) {
    return(invisible(NULL))
  }

  if (length(found) > 1L) {
    stop(
      "Multiple columns represent '",
      target,
      "' in ",
      basename(filename),
      ": ",
      paste(found, collapse = ", ")
    )
  }

  if (found != target) {
    setnames(
      x,
      old = found,
      new = target
    )
  }

  invisible(NULL)
}

# -------------------------------------------------------------------------
# Spatial helper functions
# -------------------------------------------------------------------------

nearest_coordinate <- function(query, reference) {

  reference <- sort(unique(reference))

  if (length(reference) == 1L) {
    return(rep(1L, length(query)))
  }

  midpoints <- (
    reference[-1] +
      reference[-length(reference)]
  ) / 2

  findInterval(query, midpoints) + 1L
}


infer_resolution <- function(coordinates) {

  coordinates <- sort(unique(coordinates))

  differences <- diff(coordinates)

  differences <- differences[
    is.finite(differences) &
      differences > sqrt(.Machine$double.eps)
  ]

  if (length(differences) == 0L) {
    stop("Could not infer spatial resolution.")
  }

  # Identify the most frequent spacing after removing floating-point noise
  rounded_differences <- signif(
    differences,
    digits = 8
  )

  possible_values <- unique(rounded_differences)

  counts <- tabulate(
    match(
      rounded_differences,
      possible_values
    )
  )

  resolution <- possible_values[
    which.max(counts)
  ]

  as.numeric(resolution)
}


make_regular_grid_mapping <- function(
    first_month,
    chunk_size = 500000L,
    distance_factor = 0.55
) {

  source_grid <- unique(
    first_month[, .(lon, lat)]
  )

  setorder(
    source_grid,
    lat,
    lon
  )

  source_grid[, source_id := .I]

  source_lons <- sort(unique(source_grid$lon))
  source_lats <- sort(unique(source_grid$lat))

  nx_source <- length(source_lons)
  ny_source <- length(source_lats)

  if (nx_source < 2L || ny_source < 2L) {
    stop(
      "Cannot construct a grid from fewer than two coordinates."
    )
  }

  complete_rectilinear <- (
    nrow(source_grid) == nx_source * ny_source
  )

  if (complete_rectilinear) {

    # For complete Gaussian grids, retain the original number of rows
    # and columns but distribute them regularly.
    xres <- (
      max(source_lons) - min(source_lons)
    ) / (nx_source - 1L)

    yres <- (
      max(source_lats) - min(source_lats)
    ) / (ny_source - 1L)

    nx_target <- nx_source
    ny_target <- ny_source

  } else {

    # For incomplete regular grids such as masked CHIRPS grids,
    # infer the nominal spatial resolution.
    xres <- infer_resolution(source_lons)
    yres <- infer_resolution(source_lats)

    nx_target <- as.integer(
      round(
        (max(source_lons) - min(source_lons)) / xres
      ) + 1L
    )

    ny_target <- as.integer(
      round(
        (max(source_lats) - min(source_lats)) / yres
      ) + 1L
    )
  }

  target_xmax <- min(source_lons) +
    (nx_target - 1L) * xres

  target_ymax <- min(source_lats) +
    (ny_target - 1L) * yres

  target <- rast(
    ncols = nx_target,
    nrows = ny_target,
    xmin = min(source_lons) - xres / 2,
    xmax = target_xmax + xres / 2,
    ymin = min(source_lats) - yres / 2,
    ymax = target_ymax + yres / 2,
    crs = "EPSG:4326"
  )

  message(
    "  Source grid contains ",
    format(nrow(source_grid), big.mark = ","),
    " cells; complete rectilinear grid: ",
    complete_rectilinear
  )

  if (complete_rectilinear) {

    # Fast lookup for complete rectilinear grids
    source_lookup <- matrix(
      NA_integer_,
      nrow = ny_source,
      ncol = nx_source
    )

    source_lon_id <- match(
      source_grid$lon,
      source_lons
    )

    source_lat_id <- match(
      source_grid$lat,
      source_lats
    )

    source_lookup[
      cbind(source_lat_id, source_lon_id)
    ] <- source_grid$source_id

    target_xy <- crds(
      target,
      df = TRUE
    )

    target_lon_id <- nearest_coordinate(
      target_xy$x,
      source_lons
    )

    target_lat_id <- nearest_coordinate(
      target_xy$y,
      source_lats
    )

    nearest_source_id <- source_lookup[
      cbind(target_lat_id, target_lon_id)
    ]

    rm(
      source_lookup,
      target_xy
    )

  } else {

    # General fallback for incomplete grids
    if (!requireNamespace("RANN", quietly = TRUE)) {
      stop(
        paste0(
          "Package 'RANN' is required for incomplete grids. ",
          "Install it using install.packages('RANN')."
        )
      )
    }

    n_target <- ncell(target)

    nearest_source_id <- rep(
      NA_integer_,
      n_target
    )

    nearest_distance <- rep(
      NA_real_,
      n_target
    )

    source_xy <- as.matrix(
      source_grid[, .(lon, lat)]
    )

    source_original_id <- source_grid$source_id

    # Add only points near the dateline as wrapped copies
    longitude_range <- diff(
      range(source_grid$lon)
    )

    if (longitude_range > 300) {

      edge_width <- max(
        2 * xres,
        0.1
      )

      left_edge <- which(
        source_grid$lon <=
          min(source_grid$lon) + edge_width
      )

      right_edge <- which(
        source_grid$lon >=
          max(source_grid$lon) - edge_width
      )

      source_xy <- rbind(
        source_xy,
        cbind(
          source_grid$lon[left_edge] + 360,
          source_grid$lat[left_edge]
        ),
        cbind(
          source_grid$lon[right_edge] - 360,
          source_grid$lat[right_edge]
        )
      )

      source_original_id <- c(
        source_original_id,
        source_grid$source_id[left_edge],
        source_grid$source_id[right_edge]
      )
    }

    chunk_starts <- seq.int(
      from = 1L,
      to = n_target,
      by = chunk_size
    )

    message(
      "  Calculating nearest-neighbour mapping in ",
      length(chunk_starts),
      " chunks."
    )

    for (ichunk in seq_along(chunk_starts)) {

      cell_start <- chunk_starts[ichunk]

      cell_end <- min(
        cell_start + chunk_size - 1L,
        n_target
      )

      cells <- cell_start:cell_end

      query_xy <- xyFromCell(
        target,
        cells
      )

      nn <- RANN::nn2(
        data = source_xy,
        query = query_xy,
        k = 1
      )

      nearest_source_id[cells] <- source_original_id[
        nn$nn.idx[, 1]
      ]

      nearest_distance[cells] <- nn$nn.dists[, 1]

      message(
        "    Mapping chunk ",
        ichunk,
        "/",
        length(chunk_starts)
      )

      rm(
        query_xy,
        nn
      )

      gc()
    }

    # Do not extrapolate across substantial gaps
    maximum_distance <- distance_factor *
      sqrt(xres^2 + yres^2)

    too_far <- nearest_distance > maximum_distance

    nearest_source_id[too_far] <- NA_integer_

    message(
      "  Target cells retained: ",
      round(
        100 * mean(!is.na(nearest_source_id)),
        2
      ),
      "%"
    )

    rm(
      source_xy,
      source_original_id,
      nearest_distance
    )

    gc()
  }

  list(
    target = target,
    source_grid = source_grid,
    nearest_source_id = nearest_source_id,
    xres = xres,
    yres = yres,
    complete_rectilinear = complete_rectilinear
  )
}


remap_one_variable <- function(
    month_data,
    variable,
    mapping
) {

  month_coordinates <- month_data[
    ,
    .(lon, lat)
  ]

  source_id <- mapping$source_grid[
    month_coordinates,
    on = .(lon, lat),
    source_id
  ]

  if (anyNA(source_id)) {
    stop(
      "Monthly coordinates differ from the reference spatial grid."
    )
  }

  if (anyDuplicated(source_id)) {
    stop(
      "Duplicate spatial cells found within a month."
    )
  }

  source_values <- rep(
    NA_real_,
    nrow(mapping$source_grid)
  )

  source_values[source_id] <-
    month_data[[variable]]

  target_values <- source_values[
    mapping$nearest_source_id
  ]

  output <- rast(mapping$target)

  values(output) <- target_values
  names(output) <- variable

  output
}

# -------------------------------------------------------------------------
# Find input files
# -------------------------------------------------------------------------

files <- list.files(
  input_dir,
  pattern = "\\.climate\\.RDS$",
  full.names = TRUE
)

files <- files[
  grepl(
    datasets_to_keep,
    basename(files)
  ) &
    !grepl(
      datasets_to_exclude,
      basename(files)
    )
]

# files <- rev(sort(files))
files <- sort(files)

if (length(files) == 0L) {
  stop(
    "No matching climate RDS files found in ",
    input_dir
  )
}

# -------------------------------------------------------------------------
# Process datasets
# -------------------------------------------------------------------------

for (ifile in seq_along(files)) {

  cfile <- files[ifile]

  filename_parts <- strsplit(
    basename(cfile),
    "\\."
  )[[1]]

  if (length(filename_parts) < 2L) {
    warning(
      "Cannot determine dataset name from ",
      basename(cfile)
    )
    next
  }

  cdataset <- filename_parts[2]

  message(
    "\nDataset ",
    ifile,
    "/",
    length(files),
    ": ",
    cdataset
  )

  # -----------------------------------------------------------------------
  # Read and standardise data
  # -----------------------------------------------------------------------

  cdf <- readRDS(cfile)
  setDT(cdf)

  # Precipitation
  rename_to_standard(
    cdf,
    target = "pre",
    aliases = c("MAP", "Pmm", "prate"),
    filename = cfile
  )

  # Mean/near-surface air temperature
  rename_to_standard(
    cdf,
    target = "tas",
    aliases = c("MAT", "tmp"),
    filename = cfile
  )

  # Minimum temperature
  rename_to_standard(
    cdf,
    target = "tasmin",
    aliases = c("tmin", "tmn", "Tmin"),
    filename = cfile
  )

  # Maximum temperature
  rename_to_standard(
    cdf,
    target = "tasmax",
    aliases = c("tmax", "tmx", "Tmax"),
    filename = cfile
  )

  required_columns <- c(
    "year",
    "month",
    "lon",
    "lat"
  )

  missing_columns <- setdiff(
    required_columns,
    names(cdf)
  )

  if (length(missing_columns) > 0L) {
    warning(
      "Skipping ",
      basename(cfile),
      ": missing ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
    next
  }

  variables <- setdiff(
    names(cdf),
    required_columns
  )

  variables <- variables[
    vapply(
      cdf[, ..variables],
      is.numeric,
      logical(1)
    )
  ]

  if (length(variables) == 0L) {
    warning(
      "No numeric climate variables found in ",
      basename(cfile)
    )
    next
  }

  variables <- variables[variables %in% c("tas","pre","tasmin","tasmax")]

  output_files <- file.path(
    output_dir,
    paste0(
      cdataset,
      "_",
      variables,
      "_all.years.tif"
    )
  )

  variables_to_process <- variables[
    overwrite | !file.exists(output_files)
  ]

  if (length(variables_to_process) == 0L) {

    message(
      "  All outputs already exist; skipping."
    )

    rm(cdf)
    gc()

    next
  }

  message(
    "  Variables: ",
    paste(
      variables_to_process,
      collapse = ", "
    )
  )

  columns_to_keep <- c(
    required_columns,
    variables_to_process
  )

  columns_to_drop <- setdiff(
    names(cdf),
    columns_to_keep
  )

  if (length(columns_to_drop) > 0L) {
    cdf[, (columns_to_drop) := NULL]
  }

  cdf <- cdf[
    is.finite(lon) &
      is.finite(lat) &
      is.finite(year) &
      is.finite(month) &
      month >= 1 &
      month <= 12
  ]

  cdf[, year := as.integer(year)]
  cdf[, month := as.integer(month)]
  cdf[, lon := as.numeric(lon)]
  cdf[, lat := as.numeric(lat)]

  setorder(
    cdf,
    year,
    month
  )

  # -----------------------------------------------------------------------
  # Identify monthly row ranges
  # -----------------------------------------------------------------------

  time_id <- cdf$year * 100L + cdf$month
  time_runs <- rle(time_id)

  time_end <- cumsum(
    time_runs$lengths
  )

  time_start <- c(
    1L,
    head(time_end, -1L) + 1L
  )

  times <- data.table(
    time_id = time_runs$values,
    start = time_start,
    end = time_end
  )

  times[, year := time_id %/% 100L]
  times[, month := time_id %% 100L]

  times[, date := as.Date(
    sprintf(
      "%04d-%02d-01",
      year,
      month
    )
  )]

  rm(
    time_id,
    time_runs,
    time_start,
    time_end
  )

  gc()

  # -----------------------------------------------------------------------
  # Test direct raster construction
  # -----------------------------------------------------------------------

  first_month <- cdf[
    times$start[1]:times$end[1]
  ]

  first_xyz <- as.data.frame(
    first_month[
      ,
      c(
        "lon",
        "lat",
        variables_to_process
      ),
      with = FALSE
    ]
  )

  direct_test <- tryCatch(
    {
      rast(
        first_xyz,
        type = "xyz",
        crs = "EPSG:4326"
      )
    },
    error = function(e) e
  )

  if (inherits(direct_test, "error")) {

    use_remapping <- TRUE

    message(
      "  Direct conversion failed: ",
      conditionMessage(direct_test)
    )

    message(
      "  Using regular-grid nearest-neighbour remapping."
    )

    mapping <- make_regular_grid_mapping(
      first_month = first_month,
      chunk_size = nn_chunk_size,
      distance_factor = nn_distance_factor
    )

    message(
      "  Target grid: ",
      ncol(mapping$target),
      " × ",
      nrow(mapping$target),
      "; resolution ",
      round(mapping$xres, 6),
      "° × ",
      round(mapping$yres, 6),
      "°"
    )

    rm(direct_test)

  } else {

    use_remapping <- FALSE
    reference_geometry <- direct_test

    message(
      "  Direct raster conversion is possible."
    )
  }

  rm(
    first_month,
    first_xyz
  )

  gc()

  # -----------------------------------------------------------------------
  # Create temporary directory
  # -----------------------------------------------------------------------

  temporary_dir <- tempfile(
    pattern = paste0(
      ".tmp_convert_",
      cdataset,
      "_"
    ),
    tmpdir = output_dir
  )

  dir.create(
    temporary_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  terraOptions(
    memfrac = 0.25,
    todisk = TRUE,
    tempdir = temporary_dir,
    progress = 1
  )

  temporary_files <- matrix(
    NA_character_,
    nrow = nrow(times),
    ncol = length(variables_to_process),
    dimnames = list(
      NULL,
      variables_to_process
    )
  )

  # -----------------------------------------------------------------------
  # Process one month at a time
  # -----------------------------------------------------------------------

  for (itime in seq_len(nrow(times))) {

    message(
      "  Processing month ",
      sprintf(
        "%03d/%03d",
        itime,
        nrow(times)
      ),
      ": ",
      format(
        times$date[itime],
        "%Y-%m"
      )
    )

    rows <- times$start[itime]:times$end[itime]

    month_data <- cdf[rows]

    if (!use_remapping) {

      xyz <- as.data.frame(
        month_data[
          ,
          c(
            "lon",
            "lat",
            variables_to_process
          ),
          with = FALSE
        ]
      )

      monthly_raster <- rast(
        xyz,
        type = "xyz",
        crs = "EPSG:4326"
      )

      if (!compareGeom(
        reference_geometry,
        monthly_raster,
        stopOnError = FALSE,
        messages = FALSE
      )) {
        stop(
          "Grid geometry changes in ",
          cdataset,
          " during ",
          format(
            times$date[itime],
            "%Y-%m"
          )
        )
      }
    }

    for (ivar in seq_along(variables_to_process)) {

      cvar <- variables_to_process[ivar]

      temporary_file <- file.path(
        temporary_dir,
        sprintf(
          "%04d_%s.tif",
          itime,
          cvar
        )
      )

      if (use_remapping) {

        current_raster <- remap_one_variable(
          month_data = month_data,
          variable = cvar,
          mapping = mapping
        )

      } else {

        current_raster <- monthly_raster[[cvar]]
      }

      writeRaster(
        current_raster,
        temporary_file,
        overwrite = TRUE,
        datatype = "FLT4S",
        gdal = c(
          "COMPRESS=LZW",
          "TILED=YES",
          "BIGTIFF=YES"
        )
      )

      temporary_files[itime, ivar] <- temporary_file

      rm(current_raster)
      gc()
    }

    if (exists("monthly_raster")) {
      rm(monthly_raster)
    }

    if (exists("xyz")) {
      rm(xyz)
    }

    rm(month_data)
    gc()
  }

  rm(cdf)

  if (exists("mapping")) {
    rm(mapping)
  }

  if (exists("reference_geometry")) {
    rm(reference_geometry)
  }

  gc()

  # -----------------------------------------------------------------------
  # Assemble final multiband GeoTIFFs
  # -----------------------------------------------------------------------

  for (ivar in seq_along(variables_to_process)) {

    cvar <- variables_to_process[ivar]

    output_file <- file.path(
      output_dir,
      paste0(
        cdataset,
        "_",
        cvar,
        "_all.years.tif"
      )
    )

    message(
      "  Assembling ",
      basename(output_file)
    )

    variable_raster <- rast(
      temporary_files[, ivar]
    )

    time(variable_raster) <- times$date

    names(variable_raster) <- paste0(
      cvar,
      "_",
      format(
        times$date,
        "%Y_%m"
      )
    )

    writeRaster(
      variable_raster,
      output_file,
      overwrite = overwrite,
      datatype = "FLT4S",
      gdal = c(
        "COMPRESS=DEFLATE",
        "PREDICTOR=3",
        "TILED=YES",
        "BIGTIFF=YES"
      )
    )

    rm(variable_raster)
    gc()
  }

  # Remove dataset-specific temporary files
  unlink(
    temporary_dir,
    recursive = TRUE,
    force = TRUE
  )

  rm(
    temporary_files,
    times
  )

  gc()

  message(
    "  Finished ",
    cdataset
  )
}

# scp /Users/felicien/Documents/projects/Drying.CB/scripts/convert.climate.to.rast_fast.R hpc:/data/gent/vo/000/gvo00074/felicien/R

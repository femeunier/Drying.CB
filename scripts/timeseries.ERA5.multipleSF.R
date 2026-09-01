rm(list = ls())

library(ncdf4)
library(terra)
library(data.table)
library(matrixStats)
library(dplyr)


# =============================================================================
# Settings
# =============================================================================

climate_dir <- paste0(
  "/data/gent/vo/000/gvo00074/",
  "ED_common_data/met/Tropics"
)

output_monthly <- "./outputs/ERA5_regions_monthly_climate.csv"

# FALSE excludes a month if it does not contain every expected hourly timestep.
# This still allows complete months from the incomplete year 2026 to be included.
include_incomplete_months <- FALSE

shapefiles <- list(
  Congo_Basin = "./data/shapefiles/Rainforests.shp",
  Amazon      = "./data/shapefiles/Amazon.shp"
)

dir.create(
  dirname(output_monthly),
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# General functions
# =============================================================================

extract_year <- function(file) {

  as.integer(
    sub(
      ".*_([0-9]{4})\\.nc$",
      "\\1",
      basename(file)
    )
  )
}


weighted_mean_safe <- function(
    x,
    weights
) {

  valid <- (
    is.finite(x) &
      is.finite(weights) &
      weights > 0
  )

  if (!any(valid)) {
    return(NA_real_)
  }

  weighted.mean(
    x[valid],
    weights[valid]
  )
}


days_in_month <- function(
    year,
    month
) {

  first_day <- as.Date(
    sprintf(
      "%04d-%02d-01",
      year,
      month
    )
  )

  if (month == 12) {

    next_month <- as.Date(
      sprintf(
        "%04d-01-01",
        year + 1
      )
    )

  } else {

    next_month <- as.Date(
      sprintf(
        "%04d-%02d-01",
        year,
        month + 1
      )
    )
  }

  as.integer(
    next_month -
      first_day
  )
}


# =============================================================================
# NetCDF time functions
# =============================================================================

read_netcdf_time <- function(file) {

  nc <- nc_open(file)
  on.exit(nc_close(nc))

  dimension_names <- names(nc$dim)

  time_name <- intersect(
    c("valid_time", "time"),
    dimension_names
  )

  if (length(time_name) == 0) {
    stop(
      "No valid_time or time dimension found in:\n",
      file
    )
  }

  time_name <- time_name[1]

  time_values <- nc$dim[[time_name]]$vals
  time_units <- nc$dim[[time_name]]$units

  if (
    is.null(time_units) ||
    length(time_units) == 0 ||
    is.na(time_units)
  ) {

    time_attribute <- ncatt_get(
      nc,
      time_name,
      "units"
    )

    if (isTRUE(time_attribute$hasatt)) {
      time_units <- time_attribute$value
    }
  }

  if (
    is.null(time_units) ||
    !grepl(
      " since ",
      time_units,
      fixed = TRUE
    )
  ) {
    stop(
      "Unsupported NetCDF time units in:\n",
      file,
      "\nUnits: ",
      time_units
    )
  }

  unit_name <- sub(
    " since .*",
    "",
    time_units
  )

  time_origin <- sub(
    "^.* since ",
    "",
    time_units
  )

  time_origin <- gsub(
    "T",
    " ",
    time_origin,
    fixed = TRUE
  )

  time_origin <- sub(
    " UTC$",
    "",
    time_origin
  )

  time_origin <- sub(
    "Z$",
    "",
    time_origin
  )

  time_origin <- sub(
    "([0-9]{2}:[0-9]{2}:[0-9]{2})\\.0+$",
    "\\1",
    time_origin
  )

  origin <- as.POSIXct(
    time_origin,
    tz = "UTC"
  )

  if (is.na(origin)) {
    stop(
      "Could not interpret time origin: ",
      time_origin
    )
  }

  multiplier <- switch(
    tolower(unit_name),

    "second"  = 1,
    "seconds" = 1,

    "minute"  = 60,
    "minutes" = 60,

    "hour"  = 3600,
    "hours" = 3600,

    "day"  = 86400,
    "days" = 86400,

    stop(
      "Unsupported time unit: ",
      unit_name
    )
  )

  origin + time_values * multiplier
}


# =============================================================================
# Check monthly temporal completeness
# =============================================================================

get_month_information <- function(dates) {

  dates <- as.POSIXct(
    dates,
    origin = "1970-01-01",
    tz = "UTC"
  )

  date_table <- data.table(
    datetime = dates,
    year = as.integer(
      format(
        dates,
        "%Y",
        tz = "UTC"
      )
    ),
    month = as.integer(
      format(
        dates,
        "%m",
        tz = "UTC"
      )
    )
  )

  month_information <- date_table[
    ,
    .(
      n_hours = .N,

      first_time = min(
        datetime
      ),

      last_time = max(
        datetime
      ),

      unique_hours = (
        uniqueN(datetime) == .N
      ),

      regular_hourly = if (.N > 1) {
        all(
          as.numeric(
            diff(datetime),
            units = "hours"
          ) == 1
        )
      } else {
        FALSE
      }
    ),
    by = .(
      year,
      month
    )
  ]

  month_information[
    ,
    expected_hours := vapply(
      seq_len(.N),
      function(i) {
        24L *
          days_in_month(
            year[i],
            month[i]
          )
      },
      integer(1)
    )
  ]

  month_information[
    ,
    expected_first_time := as.POSIXct(
      sprintf(
        "%04d-%02d-01 00:00:00",
        year,
        month
      ),
      tz = "UTC"
    )
  ]

  month_information[
    ,
    expected_last_time := (
      expected_first_time +
        expected_hours * 3600 -
        3600
    )
  ]

  month_information[
    ,
    complete_month := (
      n_hours == expected_hours &
        unique_hours &
        regular_hourly &
        first_time == expected_first_time &
        last_time == expected_last_time
    )
  ]

  setorder(
    month_information,
    year,
    month
  )

  month_information
}


# =============================================================================
# NetCDF coordinate functions
# =============================================================================

find_dimension_name <- function(
    nc,
    candidates
) {

  result <- intersect(
    candidates,
    names(nc$dim)
  )

  if (length(result) == 0) {
    stop(
      "Could not find NetCDF dimension: ",
      paste(
        candidates,
        collapse = " or "
      )
    )
  }

  result[1]
}


get_coordinates <- function(nc) {

  longitude_name <- find_dimension_name(
    nc,
    c("longitude", "lon")
  )

  latitude_name <- find_dimension_name(
    nc,
    c("latitude", "lat")
  )

  list(
    longitude = nc$dim[[longitude_name]]$vals,
    latitude = nc$dim[[latitude_name]]$vals
  )
}


check_variable_dimensions <- function(
    nc,
    variable
) {

  if (!variable %in% names(nc$var)) {
    stop(
      "Variable ",
      variable,
      " is absent from the NetCDF file."
    )
  }

  dimension_names <- vapply(
    nc$var[[variable]]$dim,
    function(x) {
      x$name
    },
    character(1)
  )

  valid <- (
    length(dimension_names) == 3 &&
      dimension_names[1] %in% c(
        "longitude",
        "lon"
      ) &&
      dimension_names[2] %in% c(
        "latitude",
        "lat"
      ) &&
      dimension_names[3] %in% c(
        "valid_time",
        "time"
      )
  )

  if (!valid) {
    stop(
      "Unexpected dimensions for ",
      variable,
      ": ",
      paste(
        dimension_names,
        collapse = ", "
      )
    )
  }
}


# =============================================================================
# Prepare polygon weights on the NetCDF grid
# =============================================================================

prepare_region_grid <- function(
    longitude,
    latitude,
    polygon
) {

  polygon <- project(
    polygon,
    "EPSG:4326"
  )

  polygon_extent <- as.vector(
    ext(polygon)
  )

  longitude_resolution <- median(
    abs(diff(longitude))
  )

  latitude_resolution <- median(
    abs(diff(latitude))
  )

  longitude_indices <- which(
    longitude >=
      polygon_extent[1] -
      longitude_resolution / 2 &
      longitude <=
      polygon_extent[2] +
      longitude_resolution / 2
  )

  latitude_indices <- which(
    latitude >=
      polygon_extent[3] -
      latitude_resolution / 2 &
      latitude <=
      polygon_extent[4] +
      latitude_resolution / 2
  )

  if (
    length(longitude_indices) == 0 ||
    length(latitude_indices) == 0
  ) {
    stop(
      "The polygon does not overlap the NetCDF grid."
    )
  }

  longitude_indices <- seq(
    min(longitude_indices),
    max(longitude_indices)
  )

  latitude_indices <- seq(
    min(latitude_indices),
    max(latitude_indices)
  )

  longitude_subset_raw <- longitude[
    longitude_indices
  ]

  latitude_subset_raw <- latitude[
    latitude_indices
  ]

  reverse_longitude <- (
    longitude_subset_raw[1] >
      longitude_subset_raw[
        length(longitude_subset_raw)
      ]
  )

  reverse_latitude <- (
    latitude_subset_raw[1] <
      latitude_subset_raw[
        length(latitude_subset_raw)
      ]
  )

  longitude_subset <- if (
    reverse_longitude
  ) {
    rev(longitude_subset_raw)
  } else {
    longitude_subset_raw
  }

  # Terra raster rows run from north to south
  latitude_subset <- if (
    reverse_latitude
  ) {
    rev(latitude_subset_raw)
  } else {
    latitude_subset_raw
  }

  number_longitudes <- length(
    longitude_subset
  )

  number_latitudes <- length(
    latitude_subset
  )

  template <- rast(
    ncols = number_longitudes,
    nrows = number_latitudes,

    xmin = min(longitude_subset) -
      longitude_resolution / 2,

    xmax = max(longitude_subset) +
      longitude_resolution / 2,

    ymin = min(latitude_subset) -
      latitude_resolution / 2,

    ymax = max(latitude_subset) +
      latitude_resolution / 2,

    crs = "EPSG:4326"
  )

  cell_id <- init(
    template,
    "cell"
  )

  names(cell_id) <- "cell_id"

  cell_area <- cellSize(
    template,
    unit = "km"
  )

  names(cell_area) <- "cell_area_km2"

  extracted <- extract(
    c(
      cell_id,
      cell_area
    ),
    polygon,
    exact = TRUE,
    ID = FALSE
  )

  fraction_column <- intersect(
    c("fraction", "weight"),
    names(extracted)
  )

  if (length(fraction_column) != 1) {
    stop(
      "Could not identify polygon coverage fractions."
    )
  }

  weight_table <- data.table(
    cell_id = as.integer(
      extracted[["cell_id"]]
    ),

    weight = (
      extracted[["cell_area_km2"]] *
        extracted[[fraction_column]]
    )
  )

  weight_table <- weight_table[
    is.finite(cell_id) &
      is.finite(weight) &
      weight > 0,
    .(
      weight = sum(weight)
    ),
    by = cell_id
  ]

  setorder(
    weight_table,
    cell_id
  )

  if (nrow(weight_table) == 0) {
    stop(
      "No NetCDF cells overlap the polygon."
    )
  }

  cell_latitude <- rep(
    latitude_subset,
    each = number_longitudes
  )

  list(
    longitude_start = min(
      longitude_indices
    ),

    latitude_start = min(
      latitude_indices
    ),

    number_longitudes = number_longitudes,
    number_latitudes = number_latitudes,

    reverse_longitude = reverse_longitude,
    reverse_latitude = reverse_latitude,

    selected_cells = weight_table$cell_id,
    weights = weight_table$weight,

    selected_latitudes = cell_latitude[
      weight_table$cell_id
    ]
  )
}


# =============================================================================
# Read a regional NetCDF subset directly into memory
# =============================================================================

read_netcdf_subset <- function(
    nc,
    variable,
    grid,
    time_indices
) {

  if (
    length(time_indices) > 1 &&
    any(diff(time_indices) != 1)
  ) {
    stop(
      "Requested NetCDF time indices are not contiguous."
    )
  }

  values <- ncvar_get(
    nc,
    variable,

    start = c(
      grid$longitude_start,
      grid$latitude_start,
      min(time_indices)
    ),

    count = c(
      grid$number_longitudes,
      grid$number_latitudes,
      length(time_indices)
    ),

    collapse_degen = FALSE
  )

  if (grid$reverse_longitude) {
    values <- values[
      grid$number_longitudes:1,
      ,
      ,
      drop = FALSE
    ]
  }

  if (grid$reverse_latitude) {
    values <- values[
      ,
      grid$number_latitudes:1,
      ,
      drop = FALSE
    ]
  }

  values <- matrix(
    values,
    nrow = (
      grid$number_longitudes *
        grid$number_latitudes
    ),
    ncol = length(time_indices)
  )

  values[
    grid$selected_cells,
    ,
    drop = FALSE
  ]
}


# =============================================================================
# Extraterrestrial radiation for Hargreaves PET
# =============================================================================

extraterrestrial_radiation <- function(
    latitude,
    day_of_year
) {

  latitude_radians <- (
    latitude *
      pi /
      180
  )

  inverse_distance <- (
    1 +
      0.033 *
      cos(
        2 *
          pi *
          day_of_year /
          365
      )
  )

  solar_declination <- (
    0.409 *
      sin(
        2 *
          pi *
          day_of_year /
          365 -
          1.39
      )
  )

  sunset_argument <- (
    -tan(latitude_radians) *
      tan(solar_declination)
  )

  sunset_argument <- pmax(
    -1,
    pmin(
      1,
      sunset_argument
    )
  )

  sunset_angle <- acos(
    sunset_argument
  )

  radiation_MJ <- (
    (24 * 60 / pi) *
      0.0820 *
      inverse_distance *
      (
        sunset_angle *
          sin(latitude_radians) *
          sin(solar_declination) +
          cos(latitude_radians) *
          cos(solar_declination) *
          sin(sunset_angle)
      )
  )

  # Convert MJ m-2 day-1 to equivalent evaporation in mm day-1
  0.408 * radiation_MJ
}


# =============================================================================
# Read shapefiles
# =============================================================================

regions <- lapply(
  shapefiles,
  function(file) {

    if (!file.exists(file)) {
      stop(
        "Shapefile does not exist:\n",
        file
      )
    }

    polygon <- vect(file)
    polygon <- makeValid(polygon)

    # Dissolve all features into one regional polygon
    aggregate(polygon)
  }
)

if (
  is.null(names(regions)) ||
  any(names(regions) == "")
) {
  stop(
    "Every shapefile must have a region name."
  )
}


# =============================================================================
# Locate and match ERA5 files
# =============================================================================

temperature_files <- list.files(
  climate_dir,
  pattern = paste0(
    "^ERA5_tropics_hourly_temperature_",
    "[0-9]{4}\\.nc$"
  ),
  full.names = TRUE
)

precipitation_files <- list.files(
  climate_dir,
  pattern = paste0(
    "^ERA5_tropics_hourly_precipitation_",
    "[0-9]{4}\\.nc$"
  ),
  full.names = TRUE
)

temperature_table <- data.table(
  year = extract_year(
    temperature_files
  ),
  temperature_file = temperature_files
)

precipitation_table <- data.table(
  year = extract_year(
    precipitation_files
  ),
  precipitation_file = precipitation_files
)

files_by_year <- merge(
  temperature_table,
  precipitation_table,
  by = "year"
)

setDT(files_by_year)

setorder(
  files_by_year,
  year
)

if (nrow(files_by_year) == 0) {
  stop(
    "No matching ERA5 files were found for 2020–2026."
  )
}

message(
  "Years with matching ERA5 files: ",
  min(files_by_year$year),
  "–",
  max(files_by_year$year)
)


# =============================================================================
# Process all years and regions
# =============================================================================

monthly_results <- list()
result_index <- 0L

for (iyear in seq_len(nrow(files_by_year))) {

  current_year <- files_by_year$year[iyear]

  temperature_file <- (
    files_by_year$temperature_file[iyear]
  )

  precipitation_file <- (
    files_by_year$precipitation_file[iyear]
  )

  message(
    "\nProcessing ",
    current_year,
    " (",
    iyear,
    "/",
    nrow(files_by_year),
    ")"
  )

  # ---------------------------------------------------------------------------
  # Read and compare timestamps
  # ---------------------------------------------------------------------------

  temperature_dates <- read_netcdf_time(
    temperature_file
  )

  precipitation_dates <- read_netcdf_time(
    precipitation_file
  )

  if (
    length(temperature_dates) !=
    length(precipitation_dates) ||
    !all(
      temperature_dates ==
      precipitation_dates
    )
  ) {
    stop(
      "Temperature and precipitation timestamps differ in ",
      current_year
    )
  }

  if (anyDuplicated(temperature_dates)) {
    stop(
      "Duplicated timestamps found in ",
      current_year
    )
  }

  years_in_file <- unique(
    as.integer(
      format(
        temperature_dates,
        "%Y",
        tz = "UTC"
      )
    )
  )

  if (
    length(years_in_file) != 1 ||
    years_in_file != current_year
  ) {
    stop(
      "Unexpected years in file for ",
      current_year,
      ": ",
      paste(
        years_in_file,
        collapse = ", "
      )
    )
  }

  month_information <- get_month_information(
    temperature_dates
  )

  if (any(!month_information$complete_month)) {

    message(
      "  Incomplete months: ",
      paste(
        sprintf(
          "%04d-%02d",
          month_information[
            complete_month == FALSE,
            year
          ],
          month_information[
            complete_month == FALSE,
            month
          ]
        ),
        collapse = ", "
      )
    )
  }

  if (!include_incomplete_months) {

    month_information <- month_information[
      complete_month == TRUE
    ]
  }

  if (nrow(month_information) == 0) {

    message(
      "  No complete months available; skipping ",
      current_year
    )

    next
  }

  # ---------------------------------------------------------------------------
  # Open NetCDF files
  # ---------------------------------------------------------------------------

  temperature_nc <- nc_open(
    temperature_file
  )

  precipitation_nc <- nc_open(
    precipitation_file
  )

  check_variable_dimensions(
    temperature_nc,
    "t2m"
  )

  check_variable_dimensions(
    precipitation_nc,
    "tp"
  )

  temperature_coordinates <- get_coordinates(
    temperature_nc
  )

  precipitation_coordinates <- get_coordinates(
    precipitation_nc
  )

  same_grid <- (
    isTRUE(
      all.equal(
        temperature_coordinates$longitude,
        precipitation_coordinates$longitude
      )
    ) &&
      isTRUE(
        all.equal(
          temperature_coordinates$latitude,
          precipitation_coordinates$latitude
        )
      )
  )

  if (!same_grid) {

    nc_close(temperature_nc)
    nc_close(precipitation_nc)

    stop(
      "Temperature and precipitation use different grids in ",
      current_year
    )
  }

  # ===========================================================================
  # Process every shapefile
  # ===========================================================================

  for (region_name in names(regions)) {

    message(
      "  Region: ",
      region_name
    )

    grid <- prepare_region_grid(
      longitude = temperature_coordinates$longitude,
      latitude = temperature_coordinates$latitude,
      polygon = regions[[region_name]]
    )

    message(
      "    Using ",
      length(grid$selected_cells),
      " polygon-intersecting cells"
    )

    number_cells <- length(
      grid$selected_cells
    )

    number_months <- nrow(
      month_information
    )

    monthly_temperature <- rep(
      NA_real_,
      number_months
    )

    monthly_precipitation <- rep(
      NA_real_,
      number_months
    )

    monthly_pet <- rep(
      NA_real_,
      number_months
    )

    # =========================================================================
    # Process every available month
    # =========================================================================

    for (imonth in seq_len(number_months)) {

      selected_year <- (
        month_information$year[imonth]
      )

      selected_month <- (
        month_information$month[imonth]
      )

      message(
        "    Month ",
        sprintf(
          "%04d-%02d",
          selected_year,
          selected_month
        )
      )

      time_indices <- which(
        as.integer(
          format(
            temperature_dates,
            "%Y",
            tz = "UTC"
          )
        ) == selected_year &
          as.integer(
            format(
              temperature_dates,
              "%m",
              tz = "UTC"
            )
          ) == selected_month
      )

      # -----------------------------------------------------------------------
      # Read hourly temperature for this polygon and month
      # -----------------------------------------------------------------------

      temperature_matrix <- read_netcdf_subset(
        temperature_nc,
        "t2m",
        grid,
        time_indices
      )

      # Kelvin to degrees Celsius
      temperature_matrix <- (
        temperature_matrix -
          273.15
      )

      monthly_temperature_cells <- rowMeans(
        temperature_matrix,
        na.rm = TRUE
      )

      monthly_temperature[imonth] <- (
        weighted_mean_safe(
          monthly_temperature_cells,
          grid$weights
        )
      )

      # -----------------------------------------------------------------------
      # Calculate daily Hargreaves PET
      # -----------------------------------------------------------------------

      dates_month <- temperature_dates[
        time_indices
      ]

      daily_dates <- as.Date(
        dates_month,
        tz = "UTC"
      )

      unique_daily_dates <- sort(
        unique(daily_dates)
      )

      pet_month_cells <- numeric(
        number_cells
      )

      for (
        iday in seq_along(unique_daily_dates)
      ) {

        current_date <- unique_daily_dates[iday]

        daily_indices <- which(
          daily_dates == current_date
        )

        daily_temperature <- (
          temperature_matrix[
            ,
            daily_indices,
            drop = FALSE
          ]
        )

        daily_mean <- rowMeans(
          daily_temperature,
          na.rm = TRUE
        )

        daily_min <- rowMins(
          daily_temperature,
          na.rm = TRUE
        )

        daily_max <- rowMaxs(
          daily_temperature,
          na.rm = TRUE
        )

        temperature_range <- pmax(
          0,
          daily_max -
            daily_min
        )

        day_of_year <- as.integer(
          format(
            current_date,
            "%j"
          )
        )

        radiation <- extraterrestrial_radiation(
          latitude = grid$selected_latitudes,
          day_of_year = day_of_year
        )

        daily_pet <- (
          0.0023 *
            (daily_mean + 17.8) *
            sqrt(temperature_range) *
            radiation
        )

        daily_pet <- pmax(
          0,
          daily_pet
        )

        pet_month_cells <- (
          pet_month_cells +
            daily_pet
        )
      }

      monthly_pet[imonth] <- weighted_mean_safe(
        pet_month_cells,
        grid$weights
      )

      rm(
        temperature_matrix,
        monthly_temperature_cells,
        daily_temperature,
        daily_mean,
        daily_min,
        daily_max,
        temperature_range,
        radiation,
        daily_pet,
        pet_month_cells
      )

      gc()

      # -----------------------------------------------------------------------
      # Read and sum hourly precipitation
      # -----------------------------------------------------------------------

      precipitation_matrix <- read_netcdf_subset(
        precipitation_nc,
        "tp",
        grid,
        time_indices
      )

      # Remove small negative artefacts and convert m to mm
      precipitation_matrix <- (
        pmax(
          precipitation_matrix,
          0
        ) *
          1000
      )

      precipitation_month_cells <- rowSums(
        precipitation_matrix,
        na.rm = TRUE
      )

      monthly_precipitation[imonth] <- (
        weighted_mean_safe(
          precipitation_month_cells,
          grid$weights
        )
      )

      rm(
        precipitation_matrix,
        precipitation_month_cells
      )

      gc()
    }

    # -------------------------------------------------------------------------
    # Store region-year monthly results
    # -------------------------------------------------------------------------

    result_index <- result_index + 1L

    monthly_results[[result_index]] <- data.table(
      region = region_name,
      year = month_information$year,
      month = month_information$month,

      mean_monthly_temperature_C = (
        monthly_temperature
      ),

      monthly_precipitation_mm = (
        monthly_precipitation
      ),

      monthly_PET_Hargreaves_mm = (
        monthly_pet
      ),

      complete_month = (
        month_information$complete_month
      )
    )

    rm(
      grid,
      monthly_temperature,
      monthly_precipitation,
      monthly_pet
    )

    gc()
  }

  nc_close(
    temperature_nc
  )

  nc_close(
    precipitation_nc
  )

  rm(
    temperature_nc,
    precipitation_nc,
    temperature_coordinates,
    precipitation_coordinates
  )

  gc()
}


# =============================================================================
# Combine monthly results
# =============================================================================

monthly_climate <- rbindlist(
  monthly_results,
  use.names = TRUE,
  fill = TRUE
)

if (nrow(monthly_climate) == 0) {
  stop(
    "No monthly region results were produced."
  )
}

setorder(
  monthly_climate,
  region,
  year,
  month
)


# =============================================================================
# Continuous monthly climatic water deficit
# =============================================================================

monthly_climate[
  ,
  CWD_mm := {

    result <- rep(
      NA_real_,
      .N
    )

    previous_cwd <- 0
    previous_month_id <- NA_integer_

    month_id <- (
      year * 12L +
        month
    )

    for (i in seq_len(.N)) {

      # Restart from zero following a gap
      if (
        i == 1L ||
        is.na(previous_month_id) ||
        month_id[i] !=
        previous_month_id + 1L
      ) {
        previous_cwd <- 0
      }

      if (
        is.finite(
          monthly_precipitation_mm[i]
        ) &&
        is.finite(
          monthly_PET_Hargreaves_mm[i]
        )
      ) {

        previous_cwd <- min(
          0,
          previous_cwd +
            monthly_precipitation_mm[i] -
            monthly_PET_Hargreaves_mm[i]
        )

        result[i] <- previous_cwd
      }

      previous_month_id <- month_id[i]
    }

    result
  },
  by = region
]


# =============================================================================
# Add a monthly date and reorder columns
# =============================================================================

monthly_climate[
  ,
  date := as.Date(
    sprintf(
      "%04d-%02d-01",
      year,
      month
    )
  )
]

setcolorder(
  monthly_climate,
  c(
    "region",
    "date",
    "year",
    "month",
    "mean_monthly_temperature_C",
    "monthly_precipitation_mm",
    "monthly_PET_Hargreaves_mm",
    "CWD_mm",
    "complete_month"
  )
)


# =============================================================================
# Save
# =============================================================================

fwrite(
  monthly_climate,
  output_monthly
)

print(monthly_climate)

message(
  "\nMonthly output written to: ",
  output_monthly
)

# scp /Users/felicien/Documents/projects/Drying.CB/scripts/timeseries.ERA5.multipleSF.R hpc:/kyukon/data/gent/vo/000/gvo00074/felicien/R/

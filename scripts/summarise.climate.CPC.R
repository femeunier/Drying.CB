rm(list = ls())

library(dplyr)
library(ncdf4)
library(lubridate)

# ============================================================
# Settings
# ============================================================

init.years <- 1979:as.integer(format(Sys.Date(), "%Y"))

latmin <- -30
latmax <-  30

WD <- paste0(
  "/data/gent/vo/000/gvo00074/",
  "ED_common_data/met/CPC"
)

output.file <- "./outputs/df.CPC.Tropics.climate.RDS"

dir.create(
  dirname(output.file),
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# Functions
# ============================================================

# Convert the NetCDF time coordinate to Date
read_nc_dates <- function(nc) {

  time.values <- ncvar_get(nc, "time")
  time.units  <- ncatt_get(nc, "time", "units")$value

  if (is.null(time.units) || is.na(time.units)) {
    stop("The NetCDF time variable has no units attribute")
  }

  time.unit <- trimws(
    sub("\\s+since.*$", "", time.units)
  )

  origin.text <- trimws(
    sub("^.*?since\\s+", "", time.units)
  )

  # Remove some common suffixes that interfere with parsing
  origin.text <- sub("\\s+UTC$", "", origin.text)
  origin.text <- sub("\\.0+$", "", origin.text)

  origin <- suppressWarnings(
    parse_date_time(
      origin.text,
      orders = c(
        "Ymd HMS",
        "Ymd HM",
        "Ymd"
      ),
      tz = "UTC"
    )
  )

  if (is.na(origin)) {
    stop("Could not interpret time units: ", time.units)
  }

  multiplier <- switch(
    tolower(time.unit),
    "second"  = 1,
    "seconds" = 1,
    "hour"    = 3600,
    "hours"   = 3600,
    "day"     = 86400,
    "days"    = 86400,
    stop("Unsupported time unit: ", time.unit)
  )

  as.Date(
    origin + time.values * multiplier,
    tz = "UTC"
  )
}


# Check that dimensions are ordered as longitude, latitude, time
check_dimensions <- function(nc, variable) {

  dimensions <- vapply(
    nc$var[[variable]]$dim,
    function(x) x$name,
    character(1)
  )

  expected <- c("lon", "lat", "time")

  if (!identical(dimensions, expected)) {
    stop(
      "Unexpected dimensions for ",
      variable,
      ": ",
      paste(dimensions, collapse = ", "),
      ". Expected lon, lat, time."
    )
  }
}


# Read one month for a selected latitude range
read_month_array <- function(nc, variable, poslat, time.positions) {

  check_dimensions(nc, variable)

  # The daily positions belonging to a month should be consecutive
  if (length(time.positions) > 1 &&
      any(diff(time.positions) != 1)) {
    stop("Non-consecutive daily indices for ", variable)
  }

  ncvar_get(
    nc,
    variable,
    start = c(
      1,
      min(poslat),
      min(time.positions)
    ),
    count = c(
      -1,
      length(poslat),
      length(time.positions)
    ),
    collapse_degen = FALSE
  )
}


# Sum daily precipitation while preserving completely missing cells
monthly_sum <- function(x) {

  number.valid <- apply(
    !is.na(x),
    MARGIN = c(1, 2),
    FUN = sum
  )

  result <- apply(
    x,
    MARGIN = c(1, 2),
    FUN = sum,
    na.rm = TRUE
  )

  result[number.valid == 0] <- NA_real_

  result
}


# Average daily temperatures while preserving missing cells
monthly_mean <- function(x) {

  number.valid <- apply(
    !is.na(x),
    MARGIN = c(1, 2),
    FUN = sum
  )

  result <- apply(
    x,
    MARGIN = c(1, 2),
    FUN = mean,
    na.rm = TRUE
  )

  result[number.valid == 0] <- NA_real_

  result
}


# Convert a longitude-latitude matrix to a data frame
matrix_to_dataframe <- function(x, lons, lats, value.name) {

  result <- expand.grid(
    lon = lons,
    lat = lats,
    KEEP.OUT.ATTRS = FALSE
  )

  result[[value.name]] <- as.vector(x)

  result
}


# Convert temperature to degrees Celsius if stored in Kelvin
temperature_to_celsius <- function(x, nc, variable) {

  units <- ncatt_get(nc, variable, "units")$value

  if (!is.null(units) && !is.na(units)) {

    units.clean <- tolower(trimws(units))

    if (units.clean %in% c(
      "k",
      "kelvin",
      "degree_k",
      "degrees_k"
    )) {
      x <- x - 273.15
    }
  }

  x
}


# Convert precipitation rate to daily mm if necessary
precipitation_to_daily_mm <- function(x, nc, variable) {

  units <- ncatt_get(nc, variable, "units")$value

  if (is.null(units) || is.na(units)) {
    warning("No precipitation units were found; assuming mm/day")
    return(x)
  }

  units.clean <- tolower(
    gsub(" ", "", units)
  )

  message("Precipitation units: ", units)

  # kg m-2 s-1 is equivalent to mm s-1
  if (
    grepl("kg", units.clean) &&
    grepl("s-1|/s", units.clean)
  ) {
    x <- x * 86400
  }

  x
}

# ============================================================
# Process one year
# ============================================================

process_cpc_year <- function(year) {

  message("")
  message("Processing CPC year ", year)

  precip.file <- file.path(
    WD,
    sprintf("precip.%d.nc", year)
  )

  tmin.file <- file.path(
    WD,
    sprintf("tmin.%d.nc", year)
  )

  tmax.file <- file.path(
    WD,
    sprintf("tmax.%d.nc", year)
  )

  required.files <- c(
    precip.file,
    tmin.file,
    tmax.file
  )

  missing.files <- required.files[
    !file.exists(required.files)
  ]

  if (length(missing.files) > 0) {

    warning(
      "Skipping ",
      year,
      "; missing files: ",
      paste(basename(missing.files), collapse = ", ")
    )

    return(NULL)
  }

  # Open files
  nc.pr   <- nc_open(precip.file)
  nc.tmin <- nc_open(tmin.file)
  nc.tmax <- nc_open(tmax.file)

  # Ensure that files are closed if processing fails
  on.exit({
    nc_close(nc.pr)
    nc_close(nc.tmin)
    nc_close(nc.tmax)
  }, add = TRUE)

  # Coordinates
  lons <- ncvar_get(nc.pr, "lon")
  lats <- ncvar_get(nc.pr, "lat")

  # Confirm that all products use the same grid
  if (
    !isTRUE(all.equal(lons, ncvar_get(nc.tmin, "lon"))) ||
    !isTRUE(all.equal(lons, ncvar_get(nc.tmax, "lon"))) ||
    !isTRUE(all.equal(lats, ncvar_get(nc.tmin, "lat"))) ||
    !isTRUE(all.equal(lats, ncvar_get(nc.tmax, "lat")))
  ) {
    stop("Precipitation, Tmin and Tmax grids differ in ", year)
  }

  # Tropical latitude positions
  poslat <- which(
    lats >= latmin &
      lats <= latmax
  )

  select.lat <- lats[poslat]

  # Convert 0–360 longitudes to -180–180
  select.lon <- ifelse(
    lons > 180,
    lons - 360,
    lons
  )

  # Time coordinates can differ because some CPC temperature
  # dates are missing
  dates.pr   <- read_nc_dates(nc.pr)
  dates.tmin <- read_nc_dates(nc.tmin)
  dates.tmax <- read_nc_dates(nc.tmax)

  available.months <- Reduce(
    intersect,
    list(
      unique(month(dates.pr)),
      unique(month(dates.tmin)),
      unique(month(dates.tmax))
    )
  )

  available.months <- sort(available.months)

  year.results <- vector(
    mode = "list",
    length = length(available.months)
  )

  for (imonth in seq_along(available.months)) {

    current.month <- available.months[imonth]

    message(
      sprintf(
        "Year %d, month %02d",
        year,
        current.month
      )
    )

    pos.time.pr <- which(
      year(dates.pr) == year &
        month(dates.pr) == current.month
    )

    pos.time.tmin <- which(
      year(dates.tmin) == year &
        month(dates.tmin) == current.month
    )

    pos.time.tmax <- which(
      year(dates.tmax) == year &
        month(dates.tmax) == current.month
    )

    if (
      length(pos.time.pr) == 0 ||
      length(pos.time.tmin) == 0 ||
      length(pos.time.tmax) == 0
    ) {
      next
    }

    # Avoid processing an incomplete current month
    month.start <- as.Date(
      sprintf("%d-%02d-01", year, current.month)
    )

    month.end <- ceiling_date(
      month.start,
      unit = "month"
    ) - days(1)

    if (year == as.integer(format(Sys.Date(), "%Y"))) {

      latest.available.date <- min(
        max(dates.pr[pos.time.pr]),
        max(dates.tmin[pos.time.tmin]),
        max(dates.tmax[pos.time.tmax])
      )

      if (latest.available.date < month.end) {

        message(
          "Skipping incomplete month ",
          format(month.start, "%Y-%m"),
          "; data currently end on ",
          latest.available.date
        )

        next
      }
    }

    # --------------------------------------------------------
    # Precipitation
    # --------------------------------------------------------

    precip.daily <- read_month_array(
      nc             = nc.pr,
      variable       = "precip",
      poslat          = poslat,
      time.positions  = pos.time.pr
    )

    precip.daily <- precipitation_to_daily_mm(
      precip.daily,
      nc.pr,
      "precip"
    )

    precip.monthly <- monthly_sum(
      precip.daily
    )

    rm(precip.daily)

    # --------------------------------------------------------
    # Minimum temperature
    # --------------------------------------------------------

    tmin.daily <- read_month_array(
      nc             = nc.tmin,
      variable       = "tmin",
      poslat          = poslat,
      time.positions  = pos.time.tmin
    )

    tmin.daily <- temperature_to_celsius(
      tmin.daily,
      nc.tmin,
      "tmin"
    )

    tmin.monthly <- monthly_mean(
      tmin.daily
    )

    rm(tmin.daily)

    # --------------------------------------------------------
    # Maximum temperature
    # --------------------------------------------------------

    tmax.daily <- read_month_array(
      nc             = nc.tmax,
      variable       = "tmax",
      poslat          = poslat,
      time.positions  = pos.time.tmax
    )

    tmax.daily <- temperature_to_celsius(
      tmax.daily,
      nc.tmax,
      "tmax"
    )

    tmax.monthly <- monthly_mean(
      tmax.daily
    )

    rm(tmax.daily)

    # Mean temperature
    temperature.monthly <- (
      tmin.monthly + tmax.monthly
    ) / 2

    # --------------------------------------------------------
    # Convert monthly matrices to one data frame
    # --------------------------------------------------------

    monthly.df <- matrix_to_dataframe(
      precip.monthly,
      select.lon,
      select.lat,
      "MAP"
    )

    monthly.df$Tmin <- as.vector(tmin.monthly)
    monthly.df$Tmax <- as.vector(tmax.monthly)
    monthly.df$MAT  <- as.vector(temperature.monthly)

    monthly.df$year  <- year
    monthly.df$month <- current.month

    monthly.df <- monthly.df %>%
      filter(
        !is.na(MAP) |
          !is.na(MAT)
      ) %>%
      select(
        year,
        month,
        lat,
        lon,
        MAP,
        MAT,
        Tmin,
        Tmax
      ) %>%
      arrange(
        lat,
        lon
      )

    year.results[[imonth]] <- monthly.df

    rm(
      precip.monthly,
      tmin.monthly,
      tmax.monthly,
      temperature.monthly,
      monthly.df
    )

    gc(verbose = FALSE)
  }

  bind_rows(year.results)
}

# ============================================================
# Process all years
# ============================================================

all.results <- vector(
  mode = "list",
  length = length(init.years)
)

for (iyear in seq_along(init.years)) {

  current.year <- init.years[iyear]

  message("")
  message(
    sprintf(
      "Overall progress: %d/%d (%.1f%%)",
      iyear,
      length(init.years),
      100 * iyear / length(init.years)
    )
  )

  all.results[[iyear]] <- process_cpc_year(
    current.year
  )
}

# Combine only once to avoid repeatedly copying a very large object
df.CPC <- bind_rows(all.results)

rm(all.results)
gc()

saveRDS(
  df.CPC,
  output.file
)

message("")
message("Saved: ", output.file)
message("Rows: ", format(nrow(df.CPC), big.mark = ","))
message("Size: ", round(file.size(output.file) / 1024^2, 1), " MB")

# scp /Users/felicien/Documents/projects/Drying.CB/scripts/summarise.climate.CPC.R hpc:/data/gent/vo/000/gvo00074/felicien/R

rm(list = ls())

library(dplyr)
library(ncdf4)
library(lubridate)

# ============================================================
# Settings
# ============================================================

WD <- paste0(
  "/data/gent/vo/000/gvo00074/ED_common_data/met/JRA3Q"
)

temperature.dir <- file.path(
  WD,
  "temperature"
)

precipitation.dir <- file.path(
  WD,
  "precipitation"
)

output.file <- "./outputs/df.JRA3Q.Tropics.climate.RDS"

dir.create(
  dirname(output.file),
  recursive = TRUE,
  showWarnings = FALSE
)

latmin <- -30
latmax <-  30

init.years <- 1947:as.integer(format(Sys.Date(), "%Y"))

temperature.variable <- "tmp2m-hgt-an-ll125-mn"
precipitation.variable <- "tprate1have-sfc-fc-ll125-mn"

# ============================================================
# Functions
# ============================================================

# Extract YYYYMM from the JRA-3Q filename
extract_year_month <- function(filename) {

  filename <- basename(filename)

  position <- regexpr(
    "[12][0-9]{5}0100_",
    filename
  )

  if (position < 0) {
    return(NA_character_)
  }

  substr(
    filename,
    position,
    position + 5
  )
}


# Read a JRA-3Q variable and put dimensions in:
# longitude × latitude × time order
read_jra_variable <- function(nc, variable) {

  if (!variable %in% names(nc$var)) {
    stop(
      "Variable '",
      variable,
      "' not found. Available variables: ",
      paste(names(nc$var), collapse = ", ")
    )
  }

  dimension.names <- vapply(
    nc$var[[variable]]$dim,
    function(x) x$name,
    character(1)
  )

  expected.dimensions <- c(
    "lon",
    "lat",
    "time"
  )

  if (!all(expected.dimensions %in% dimension.names)) {
    stop(
      "Unexpected dimensions for ",
      variable,
      ": ",
      paste(dimension.names, collapse = ", ")
    )
  }

  values <- ncvar_get(
    nc,
    variable,
    collapse_degen = FALSE
  )

  # Reorder dimensions if necessary
  dimension.order <- match(
    expected.dimensions,
    dimension.names
  )

  values <- aperm(
    values,
    dimension.order
  )

  values
}


# Convert temperature to degrees Celsius if needed
temperature_to_celsius <- function(values, nc, variable) {

  units <- ncatt_get(
    nc,
    variable,
    "units"
  )$value

  if (
    !is.null(units) &&
    !is.na(units) &&
    tolower(trimws(units)) %in%
    c("k", "kelvin", "degree_k", "degrees_k")
  ) {
    values <- values - 273.15
  }

  values
}


# Convert the mean precipitation rate to monthly mm
precipitation_to_monthly_mm <- function(
    values,
    date,
    nc,
    variable
) {

  units <- ncatt_get(
    nc,
    variable,
    "units"
  )$value

  if (!is.null(units) && !is.na(units)) {
    message("Precipitation units: ", units)
  }

  seconds.in.month <- days_in_month(date) * 86400

  # kg m-2 is equivalent to mm of water
  values * seconds.in.month
}

# ============================================================
# Locate downloaded files
# ============================================================

temperature.files <- list.files(
  temperature.dir,
  pattern = paste0(
    "^jra3q-ms-mn\\.anl_surf125\\.",
    "0_0_0\\.tmp2m-hgt-an-ll125-mn\\.",
    "[0-9]{10}_[0-9]{10}\\.nc$"
  ),
  full.names = TRUE,
  recursive = TRUE
)

precipitation.files <- list.files(
  precipitation.dir,
  pattern = paste0(
    "^jra3q-ms-mn\\.fcst_phy2m125\\.",
    "0_1_52\\.tprate1have-sfc-fc-ll125-mn\\.",
    "[0-9]{10}_[0-9]{10}\\.nc$"
  ),
  full.names = TRUE,
  recursive = TRUE
)

message(
  "Temperature files found: ",
  length(temperature.files)
)

message(
  "Precipitation files found: ",
  length(precipitation.files)
)

if (length(temperature.files) == 0) {
  stop("No JRA-3Q temperature files found in ", temperature.dir)
}

if (length(precipitation.files) == 0) {
  stop("No JRA-3Q precipitation files found in ", precipitation.dir)
}

# ============================================================
# Match temperature and precipitation by YYYYMM
# ============================================================

temperature.table <- data.frame(
  ym = vapply(
    temperature.files,
    extract_year_month,
    character(1)
  ),
  temperature.file = temperature.files,
  stringsAsFactors = FALSE
)

precipitation.table <- data.frame(
  ym = vapply(
    precipitation.files,
    extract_year_month,
    character(1)
  ),
  precipitation.file = precipitation.files,
  stringsAsFactors = FALSE
)

temperature.table <- temperature.table %>%
  filter(!is.na(ym)) %>%
  distinct(ym, .keep_all = TRUE)

precipitation.table <- precipitation.table %>%
  filter(!is.na(ym)) %>%
  distinct(ym, .keep_all = TRUE)

# Report missing file pairs
missing.precipitation <- setdiff(
  temperature.table$ym,
  precipitation.table$ym
)

missing.temperature <- setdiff(
  precipitation.table$ym,
  temperature.table$ym
)

if (length(missing.precipitation) > 0) {
  warning(
    "Temperature exists but precipitation is missing for: ",
    paste(missing.precipitation, collapse = ", ")
  )
}

# scp /Users/felicien/Documents/projects/Drying.CB/scripts/summarise.climate.JRA3Q.R hpc:/data/gent/vo/000/gvo00074/felicien/R


if (length(missing.temperature) > 0) {
  warning(
    "Precipitation exists but temperature is missing for: ",
    paste(missing.temperature, collapse = ", ")
  )
}

file.table <- inner_join(
  temperature.table,
  precipitation.table,
  by = "ym"
) %>%
  mutate(
    year = as.integer(substr(ym, 1, 4)),
    month = as.integer(substr(ym, 5, 6)),
    date = as.Date(
      paste0(ym, "01"),
      format = "%Y%m%d"
    )
  ) %>%
  filter(year %in% init.years) %>%
  arrange(date)

message(
  "Matched months to process: ",
  nrow(file.table)
)

# ============================================================
# Process one month
# ============================================================

process_jra_month <- function(
    temperature.file,
    precipitation.file,
    date
) {

  nc.temperature <- nc_open(
    temperature.file
  )

  nc.precipitation <- nc_open(
    precipitation.file
  )

  on.exit(
    {
      nc_close(nc.temperature)
      nc_close(nc.precipitation)
    },
    add = TRUE
  )

  # ----------------------------------------------------------
  # Coordinates
  # ----------------------------------------------------------

  lons.temperature <- ncvar_get(
    nc.temperature,
    "lon"
  )

  lats.temperature <- ncvar_get(
    nc.temperature,
    "lat"
  )

  lons.precipitation <- ncvar_get(
    nc.precipitation,
    "lon"
  )

  lats.precipitation <- ncvar_get(
    nc.precipitation,
    "lat"
  )

  if (
    !isTRUE(all.equal(
      lons.temperature,
      lons.precipitation
    )) ||
    !isTRUE(all.equal(
      lats.temperature,
      lats.precipitation
    ))
  ) {
    stop(
      "Temperature and precipitation grids differ for ",
      format(date, "%Y-%m")
    )
  }

  lons <- lons.temperature
  lats <- lats.temperature

  poslat <- which(
    lats >= latmin &
      lats <= latmax
  )

  selected.lats <- lats[poslat]

  # Convert longitude from 0–360 to -180–180
  selected.lons <- ifelse(
    lons > 180,
    lons - 360,
    lons
  )

  # ----------------------------------------------------------
  # Temperature
  # ----------------------------------------------------------

  temperature <- read_jra_variable(
    nc.temperature,
    temperature.variable
  )

  temperature <- temperature[
    ,
    poslat,
    1,
    drop = TRUE
  ]

  temperature <- temperature_to_celsius(
    temperature,
    nc.temperature,
    temperature.variable
  )

  # ----------------------------------------------------------
  # Precipitation
  # ----------------------------------------------------------

  precipitation.rate <- read_jra_variable(
    nc.precipitation,
    precipitation.variable
  )

  precipitation.rate <- precipitation.rate[
    ,
    poslat,
    1,
    drop = TRUE
  ]

  precipitation.monthly <- precipitation_to_monthly_mm(
    precipitation.rate,
    date,
    nc.precipitation,
    precipitation.variable
  )

  # ----------------------------------------------------------
  # Convert matrices to data frame
  # ----------------------------------------------------------

  result <- expand.grid(
    lon = selected.lons,
    lat = selected.lats,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  result$MAP <- as.vector(
    precipitation.monthly
  )

  result$MAT <- as.vector(
    temperature
  )

  result$year <- year(date)
  result$month <- month(date)

  result %>%
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
      MAT
    ) %>%
    arrange(
      lat,
      lon
    )
}

# ============================================================
# Process all monthly files
# ============================================================

monthly.results <- vector(
  mode = "list",
  length = nrow(file.table)
)

for (ifile in seq_len(nrow(file.table))) {

  message("")
  message(
    sprintf(
      "Processing %s — %d/%d (%.1f%%)",
      format(file.table$date[ifile], "%Y-%m"),
      ifile,
      nrow(file.table),
      100 * ifile / nrow(file.table)
    )
  )

  monthly.results[[ifile]] <- tryCatch(
    {
      process_jra_month(
        temperature.file =
          file.table$temperature.file[ifile],
        precipitation.file =
          file.table$precipitation.file[ifile],
        date =
          file.table$date[ifile]
      )
    },
    error = function(e) {

      warning(
        "Failed to process ",
        file.table$ym[ifile],
        ": ",
        conditionMessage(e)
      )

      NULL
    }
  )
}

# Combine only once
df.JRA3Q <- bind_rows(
  monthly.results
)

rm(monthly.results)
gc()

# ============================================================
# Save
# ============================================================

saveRDS(
  df.JRA3Q,
  output.file
)

message("")
message("Saved: ", output.file)

message(
  "Rows: ",
  format(
    nrow(df.JRA3Q),
    big.mark = ","
  )
)

message(
  "Period: ",
  min(df.JRA3Q$year),
  "-",
  sprintf("%02d", min(df.JRA3Q$month[df.JRA3Q$year ==
                                       min(df.JRA3Q$year)])),
  " to ",
  max(df.JRA3Q$year),
  "-",
  sprintf("%02d", max(df.JRA3Q$month[df.JRA3Q$year ==
                                       max(df.JRA3Q$year)]))
)

message(
  "File size: ",
  round(
    file.size(output.file) / 1024^2,
    1
  ),
  " MB"
)

# scp /Users/felicien/Documents/projects/Drying.CB/scripts/summarise.climate.JRA3Q.R hpc:/data/gent/vo/000/gvo00074/felicien/R

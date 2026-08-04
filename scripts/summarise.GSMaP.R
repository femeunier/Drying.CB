rm(list = ls())

library(terra)
library(dplyr)
library(lubridate)

#===========================================================
# Configuration
#===========================================================

input_dir <- paste0(
  "/data/gent/vo/000/gvo00074/ED_common_data/met/GSMaP"
)

output_dir <- "./outputs/GSMaP_monthly_by_year"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

final_output <- "./outputs/df.GSMaP_Gauge_v8.Tropics.climate.RDS"

latmin <- -30
latmax <-  30

# GSMaP native resolution = 0.1 degrees
# factor 2 produces a 0.2-degree grid, matching your CHIRPS output
spatial_factor <- 2

# Set to 1 to retain native 0.1-degree resolution
# spatial_factor <- 1

#===========================================================
# GSMaP grid definition
#===========================================================

nlon <- 3600
nlat <- 1200
ncell_gsmap <- nlon * nlat

# Global domain:
# longitude: 0–360 degrees
# latitude: 60°S–60°N
gsmap_template <- rast(
  nrows = nlat,
  ncols = nlon,
  xmin = 0,
  xmax = 360,
  ymin = -60,
  ymax = 60,
  crs = "EPSG:4326"
)

#===========================================================
# Locate files
#===========================================================

files <- list.files(
  input_dir,
  pattern = paste0(
    "^gsmap_gauge\\.[0-9]{6}",
    "\\.0\\.1d\\.monthly.*\\.dat\\.gz$"
  ),
  full.names = TRUE,
  recursive = TRUE
)

if (length(files) == 0) {
  stop("No GSMaP monthly files found in: ", input_dir)
}

yyyymm <- sub(
  ".*gsmap_gauge\\.([0-9]{6}).*",
  "\\1",
  basename(files)
)

file_info <- data.frame(
  file = files,
  date = as.Date(paste0(yyyymm, "01"), format = "%Y%m%d")
) %>%
  mutate(
    year = year(date),
    month = month(date)
  ) %>%
  arrange(date)

message(
  "Found ",
  nrow(file_info),
  " monthly GSMaP files covering ",
  format(min(file_info$date), "%Y-%m"),
  " to ",
  format(max(file_info$date), "%Y-%m")
)

# Check for duplicate months
duplicate_dates <- file_info$date[duplicated(file_info$date)]

if (length(duplicate_dates) > 0) {
  stop(
    "Multiple GSMaP files found for: ",
    paste(unique(duplicate_dates), collapse = ", ")
  )
}

#===========================================================
# Process files year by year
#===========================================================

for (year_i in unique(file_info$year)) {

  yearly_output <- file.path(
    output_dir,
    sprintf("GSMaP_Gauge_v8_monthly_%d.rds", year_i)
  )

  # Skip completed historical years.
  # Recalculate the current year because new months may have appeared.
  if (
    file.exists(yearly_output) &&
    year_i < year(Sys.Date())
  ) {
    message("Already processed: ", year_i)
    next
  }

  files_year <- file_info %>%
    filter(year == year_i)

  message(
    "Processing ",
    year_i,
    " — ",
    nrow(files_year),
    " monthly files"
  )

  monthly_list <- vector("list", nrow(files_year))

  for (i in seq_len(nrow(files_year))) {

    file_i <- files_year$file[i]
    date_i <- files_year$date[i]

    message(
      "  ",
      format(date_i, "%Y-%m"),
      " (",
      i,
      "/",
      nrow(files_year),
      ")"
    )

    #-------------------------------------------------------
    # Read compressed little-endian binary file
    #-------------------------------------------------------

    binary_data <- local({

      connection <- gzfile(file_i, open = "rb")
      on.exit(close(connection))

      # Field 1: monthly mean precipitation rate [mm/hour]
      rain_rate <- readBin(
        connection,
        what = "numeric",
        n = ncell_gsmap,
        size = 4,
        endian = "little"
      )

      # Field 2: number of valid hourly samples
      n_samples <- readBin(
        connection,
        what = "numeric",
        n = ncell_gsmap,
        size = 4,
        endian = "little"
      )

      list(
        rain_rate = rain_rate,
        n_samples = n_samples
      )
    })

    # Verify that both complete fields were read
    if (
      length(binary_data$rain_rate) != ncell_gsmap ||
      length(binary_data$n_samples) != ncell_gsmap
    ) {
      stop(
        "Unexpected binary file size for: ",
        file_i,
        "\nExpected two fields of ",
        format(ncell_gsmap, big.mark = ","),
        " values."
      )
    }

    #-------------------------------------------------------
    # Convert mean rate to monthly precipitation
    #-------------------------------------------------------

    rain_rate <- binary_data$rain_rate
    n_samples <- binary_data$n_samples

    # Negative precipitation values indicate missing data.
    invalid <- rain_rate < 0 | n_samples <= 0

    # mm/hour × number of valid hours = mm/month
    monthly_precip <- rain_rate * n_samples
    monthly_precip[invalid] <- NA_real_

    rm(binary_data, rain_rate, n_samples, invalid)
    gc(FALSE)

    #-------------------------------------------------------
    # Convert to SpatRaster
    #-------------------------------------------------------

    precipitation_raster <- gsmap_template
    values(precipitation_raster) <- monthly_precip

    rm(monthly_precip)
    gc(FALSE)

    # Extract tropical region
    precipitation_raster <- crop(
      precipitation_raster,
      ext(0, 360, latmin, latmax)
    )

    #-------------------------------------------------------
    # Spatial aggregation
    #-------------------------------------------------------

    if (spatial_factor > 1) {

      precipitation_raster <- aggregate(
        precipitation_raster,
        fact = spatial_factor,
        fun = "mean",
        na.rm = TRUE
      )
    }

    #-------------------------------------------------------
    # Convert raster to data frame
    #-------------------------------------------------------

    df_month <- as.data.frame(
      precipitation_raster,
      xy = TRUE,
      na.rm = TRUE
    )

    names(df_month) <- c("lon", "lat", "Pmm")

    monthly_list[[i]] <- df_month %>%
      mutate(
        # Convert 0–360° to −180–180°
        lon = if_else(lon > 180, lon - 360, lon),
        lon = round(lon, 3),
        lat = round(lat, 3),
        month = month(date_i),
        year = year(date_i)
      ) %>%
      select(lon, lat, Pmm, month, year)

    rm(precipitation_raster, df_month)
    gc(FALSE)
  }

  df_year <- bind_rows(monthly_list)

  saveRDS(
    df_year,
    yearly_output,
    compress = FALSE
  )

  message(
    "Saved: ",
    yearly_output,
    " — ",
    format(nrow(df_year), big.mark = ","),
    " rows"
  )

  rm(monthly_list, df_year)
  gc()
}

#===========================================================
# Combine yearly files
#===========================================================

yearly_files <- list.files(
  output_dir,
  pattern = "^GSMaP_Gauge_v8_monthly_[0-9]{4}\\.rds$",
  full.names = TRUE
)

yearly_years <- as.integer(
  sub(
    "^GSMaP_Gauge_v8_monthly_([0-9]{4})\\.rds$",
    "\\1",
    basename(yearly_files)
  )
)

yearly_files <- yearly_files[order(yearly_years)]

message(
  "Combining ",
  length(yearly_files),
  " yearly files"
)

df.all <- bind_rows(
  lapply(yearly_files, readRDS)
)

saveRDS(
  df.all,
  final_output,
  compress = FALSE
)

message("Final dataset saved to: ", final_output)

message(
  "Number of rows: ",
  format(nrow(df.all), big.mark = ",")
)

message(
  "Period: ",
  sprintf(
    "%d-%02d",
    min(df.all$year),
    min(df.all$month[df.all$year == min(df.all$year)])
  ),
  " to ",
  sprintf(
    "%d-%02d",
    max(df.all$year),
    max(df.all$month[df.all$year == max(df.all$year)])
  )
)

# scp /Users/felicien/Documents/projects/Drying.CB/scripts/summarise.GSMaP.R hpc:/data/gent/vo/000/gvo00074/felicien/R

rm(list = ls())

# ============================================================
# Download CPC precipitation, Tmin and Tmax
# ============================================================

output_dir <- paste0(
  "/data/gent/vo/000/gvo00074/",
  "ED_common_data/met/CPC"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

options(timeout = 3600)

# Available from 1979 onward
years <- 1979:as.integer(format(Sys.Date(), "%Y"))
current_year <- as.integer(format(Sys.Date(), "%Y"))

# Refresh the current-year files because they are incomplete
# and updated regularly by NOAA
refresh_current_year <- TRUE

max_attempts <- 3

# Product definitions
products <- data.frame(
  variable = c(
    "precip",
    "tmin",
    "tmax"
  ),
  base_url = c(
    "https://downloads.psl.noaa.gov/Datasets/cpc_global_precip",
    "https://downloads.psl.noaa.gov/Datasets/cpc_global_temp",
    "https://downloads.psl.noaa.gov/Datasets/cpc_global_temp"
  ),
  stringsAsFactors = FALSE
)

# Total number of requested files
total_downloads <- length(years) * nrow(products)
download_number <- 0L

# ============================================================
# Download loop
# ============================================================

for (year in years) {

  for (iproduct in seq_len(nrow(products))) {

    download_number <- download_number + 1L

    variable <- products$variable[iproduct]
    base_url <- products$base_url[iproduct]

    filename <- sprintf(
      "%s.%d.nc",
      variable,
      year
    )

    url <- paste(
      base_url,
      filename,
      sep = "/"
    )

    dest <- file.path(
      output_dir,
      filename
    )

    temp <- paste0(
      dest,
      ".part"
    )

    is_current_year <- year == current_year

    # Skip complete historical files
    if (
      file.exists(dest) &&
      file.size(dest) > 1000 &&
      !(is_current_year && refresh_current_year)
    ) {

      message("Already downloaded: ", filename)
      next
    }

    message("")
    message(
      sprintf(
        "Downloading %s (%d/%d)",
        filename,
        download_number,
        total_downloads
      )
    )

    success <- FALSE

    for (attempt in seq_len(max_attempts)) {

      if (file.exists(temp)) {
        unlink(temp)
      }

      message(
        sprintf(
          "Attempt %d/%d",
          attempt,
          max_attempts
        )
      )

      success <- tryCatch(
        {
          status <- download.file(
            url      = url,
            destfile = temp,
            method   = "libcurl",
            mode     = "wb",
            quiet    = FALSE
          )

          if (!identical(status, 0L)) {
            stop("download.file returned status ", status)
          }

          if (
            !file.exists(temp) ||
            file.size(temp) <= 1000
          ) {
            stop("Downloaded file is missing or empty")
          }

          TRUE
        },
        error = function(e) {

          warning(
            "Failed to download ",
            filename,
            ": ",
            conditionMessage(e)
          )

          FALSE
        }
      )

      if (isTRUE(success)) {
        break
      }

      if (attempt < max_attempts) {
        Sys.sleep(5 * attempt)
      }
    }

    if (isTRUE(success)) {

      # Only remove the previous current-year file after the
      # replacement has been downloaded successfully
      if (file.exists(dest)) {
        unlink(dest)
      }

      renamed <- file.rename(
        from = temp,
        to   = dest
      )

      if (!renamed) {
        stop(
          "Downloaded but could not rename ",
          temp,
          " to ",
          dest
        )
      }

      message(
        "Completed: ",
        filename,
        " — ",
        round(file.size(dest) / 1024^2, 1),
        " MB"
      )

    } else {

      if (file.exists(temp)) {
        unlink(temp)
      }

      warning(
        "Giving up after ",
        max_attempts,
        " attempts: ",
        filename
      )
    }
  }
}

message("")
message("Downloads finished.")

# scp /Users/felicien/Documents/projects/Drying.CB/scripts/download.CPC.R hpc:/data/gent/vo/000/gvo00074/felicien/R

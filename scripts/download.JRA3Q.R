rm(list = ls())

# ============================================================
# Download monthly JRA-3Q temperature and precipitation
#
# Temperature:
#   Product:  anl_surf125
#   Variable: 2-metre temperature
#
# Precipitation:
#   Product:  fcst_phy2m125
#   Variable: Total precipitation rate
# ============================================================

options(timeout = 3600)

# ---------- Settings -----------------------------------------

base.url <- paste0(
  "https://tds.gdex.ucar.edu/",
  "thredds/fileServer/files/d640002"
)

output.root <- paste0(
  "/data/gent/vo/000/gvo00074/ED_common_data/met/JRA3Q"
)

start.date <- as.Date("1947-09-01")

# Try through the previous completed month.
# Months not yet available are skipped.
this.month <- as.Date(format(Sys.Date(), "%Y-%m-01"))
end.date   <- seq(this.month, by = "-1 month", length.out = 2)[2]

max.tries <- 3

# ---------- Output directories -------------------------------

temp.dir   <- file.path(output.root, "temperature")
precip.dir <- file.path(output.root, "precipitation")

dir.create(temp.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(precip.dir, recursive = TRUE, showWarnings = FALSE)

# ---------- Functions ----------------------------------------

last_day_of_month <- function(x) {

  next.month <- seq(x, by = "1 month", length.out = 2)[2]
  next.month - 1
}


make_filename <- function(date, variable) {

  ym <- format(date, "%Y%m")
  ymd.last <- format(last_day_of_month(date), "%Y%m%d")

  if (variable == "temperature") {

    paste0(
      "jra3q-ms-mn.anl_surf125.",
      "0_0_0.tmp2m-hgt-an-ll125-mn.",
      ym, "0100_", ymd.last, "18.nc"
    )

  } else if (variable == "precipitation") {

    paste0(
      "jra3q-ms-mn.fcst_phy2m125.",
      "0_1_52.tprate1have-sfc-fc-ll125-mn.",
      ym, "0100_", ymd.last, "23.nc"
    )

  } else {

    stop("Unknown variable: ", variable)
  }
}


download_jra3q <- function(date, variable, max.tries = 3) {

  ym <- format(date, "%Y%m")

  if (variable == "temperature") {

    product <- "anl_surf125"
    output.dir <- temp.dir

  } else {

    product <- "fcst_phy2m125"
    output.dir <- precip.dir
  }

  filename <- make_filename(date, variable)

  url <- paste(
    base.url,
    product,
    ym,
    filename,
    sep = "/"
  )

  destination <- file.path(output.dir, filename)
  temporary <- paste0(destination, ".part")

  # Do not download an existing, non-empty file again
  if (file.exists(destination) &&
      file.info(destination)$size > 1000) {

    message("Already present: ", filename)
    return("existing")
  }

  if (file.exists(temporary)) {
    unlink(temporary)
  }

  for (attempt in seq_len(max.tries)) {

    message(
      "[", variable, "] ",
      ym,
      " — attempt ",
      attempt,
      "/",
      max.tries
    )

    success <- tryCatch({

      status <- suppressWarnings(
        download.file(
          url       = url,
          destfile  = temporary,
          mode      = "wb",
          method    = "libcurl",
          quiet     = FALSE
        )
      )

      identical(status, 0L) &&
        file.exists(temporary) &&
        file.info(temporary)$size > 1000

    }, error = function(e) {

      message("Download error: ", conditionMessage(e))
      FALSE
    })

    if (success) {

      if (!file.rename(temporary, destination)) {
        stop("Could not rename temporary file: ", temporary)
      }

      message("Saved: ", destination)
      return("downloaded")
    }

    if (file.exists(temporary)) {
      unlink(temporary)
    }

    if (attempt < max.tries) {
      Sys.sleep(2 * attempt)
    }
  }

  warning("Unavailable or failed: ", url)
  "failed"
}

# ---------- Download all months -------------------------------

months <- rev(seq(
  from = start.date,
  to   = end.date,
  by   = "month",
))

download.log <- vector("list", length(months) * 2)
ilog <- 0L

for (date in months) {

  # Required because looping over Date objects can remove the Date class
  date <- as.Date(date, origin = "1970-01-01")

  for (variable in c("temperature", "precipitation")) {

    ilog <- ilog + 1L

    status <- download_jra3q(
      date       = date,
      variable   = variable,
      max.tries  = max.tries
    )

    download.log[[ilog]] <- data.frame(
      date     = date,
      year     = as.integer(format(date, "%Y")),
      month    = as.integer(format(date, "%m")),
      variable = variable,
      status   = status
    )
  }
}

download.log <- do.call(rbind, download.log)

write.csv(
  download.log,
  file.path(output.root, "JRA3Q_download_log.csv"),
  row.names = FALSE
)

print(table(download.log$variable, download.log$status))

# scp /Users/felicien/Documents/projects/Drying.CB/scripts/download.JRA3Q.R hpc:/data/gent/vo/000/gvo00074/felicien/R

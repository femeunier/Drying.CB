rm(list = ls())

# install.packages(c("terra","rvest","xml2","lubridate"))

library(terra)
library(dplyr)
library(rvest)
library(xml2)
library(lubridate)
library(dplyr)
library(sf)

A <- read_sf("/home/femeunier/Downloads/africa_rfe.20010101.shp/africa_rfe.20010101.shp")

plot(A)

# rfe2_monthly_sum <- function(
#     base_url   = "https://ftp.cpc.ncep.noaa.gov/fews/fewsdata/africa/rfe2/geotiff/",
#     from       = as.Date("2001-01-01"),
#     to         = as.Date("2001-12-31"),
#     out_dir    = "rfe2_monthly",
#     overwrite  = FALSE,
#     na_to_zero = TRUE,
#     mask       = NULL,
#     quiet      = FALSE,
#     coverage.name = "_coverage.csv",
#     max_retries = 3,
#     timeout_sec = 120,          # longer than default 60
#     method      = "libcurl"     # more robust than default on many systems
# ) {
#   dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
#
#   # 1) Scrape directory for .tif.zip links
#   pg <- rvest::read_html(base_url)
#   hrefs <- rvest::html_attr(rvest::html_elements(pg, "a"), "href")
#   hrefs <- hrefs[grepl("^africa_rfe\\.[0-9]{8}\\.tif\\.zip$", hrefs)]
#   if (!length(hrefs)) stop("No .tif.zip files found at: ", base_url)
#
#   # 2) Extract dates and filter by range
#   dates_all <- as.Date(sub("^africa_rfe\\.(\\d{8})\\.tif\\.zip$", "\\1", hrefs), "%Y%m%d")
#   keep <- which(dates_all >= from & dates_all <= to)
#   if (!length(keep)) stop("No files in the requested date range.")
#   hrefs <- hrefs[keep]
#   dates <- dates_all[keep]
#   urls  <- paste0(base_url, hrefs)
#
#   # Helper: robust downloader returning path to a valid extracted .tif or NULL
#   read_daily_rast <- function(u) {
#     zfile  <- tempfile(fileext = ".zip")
#     tifdir <- tempfile()
#     on.exit({
#       try(unlink(zfile), silent = TRUE)
#       try(unlink(tifdir, recursive = TRUE, force = TRUE), silent = TRUE)
#       terra::tmpFiles(remove = TRUE)
#     }, add = TRUE)
#
#     # retry loop
#     ok <- FALSE
#     for (k in seq_len(max_retries)) {
#       if (!quiet) message("  [try ", k, "/", max_retries, "] ", basename(u))
#       # extend timeout for this call
#       old_to <- getOption("timeout")
#       options(timeout = max(timeout_sec, old_to %||% 60))
#       res <- try(utils::download.file(u, zfile, mode = "wb", quiet = quiet, method = method), silent = TRUE)
#       options(timeout = old_to)
#       if (inherits(res, "try-error")) next
#       # check unzip validity
#       tf <- try(utils::unzip(zfile, exdir = tifdir), silent = TRUE)
#       if (inherits(tf, "try-error") || length(tf) == 0) next
#       ok <- TRUE
#       break
#     }
#     if (!ok) return(NULL)
#
#     r <- try(terra::rast(tf[1]), silent = TRUE)
#     if (inherits(r, "try-error")) return(NULL)
#
#     if (!is.null(mask)) {
#       if (inherits(mask, "SpatVector"))    r <- terra::crop(r, mask) |> terra::mask(mask)
#       else if (inherits(mask, "SpatRaster")) r <- terra::crop(r, mask) |> terra::mask(mask)
#     }
#     if (na_to_zero) r[is.na(r)] <- 0
#     r
#   }
#
#   # 3) Iterate by month and keep coverage
#   ym <- format(dates, "%Y%m")
#   months <- sort(unique(ym))
#   out_paths <- character(length(months)); names(out_paths) <- months
#
#   coverage <- data.frame(
#     month = months,
#     days_expected = as.integer(NA),
#     days_ok       = as.integer(NA),
#     days_missing  = as.integer(NA),
#     missing_dates = I(vector("list", length(months))),
#     stringsAsFactors = FALSE
#   )
#
#   for (ii in seq_along(months)) {
#     m <- months[ii]
#     idx <- which(ym == m)
#     coverage$days_expected[ii] <- length(idx)
#     if (!length(idx)) next
#
#     if (!quiet) message("\n-- Summing month ", m, " (", length(idx), " days)")
#
#     # read one successful day to anchor grid; otherwise skip month
#     r_first <- NULL; first_j <- NA_integer_
#     for (j in idx) {
#       rj <- read_daily_rast(urls[j])
#       if (!is.null(rj)) { r_first <- rj; first_j <- j; break }
#     }
#     if (is.null(r_first)) {
#       if (!quiet) message("  No valid days for ", m, " → skipping write.")
#       coverage$days_ok[ii] <- 0
#       coverage$days_missing[ii] <- length(idx)
#       coverage$missing_dates[[ii]] <- as.list(as.character(dates[idx]))
#       next
#     }
#
#     sumr <- r_first
#     ok_days <- 1L
#     missed <- character(0)
#
#     # loop remaining days of the month
#     for (j in setdiff(idx, first_j)) {
#       rj <- read_daily_rast(urls[j])
#       if (is.null(rj)) {
#         missed <- c(missed, as.character(dates[j]))
#         next
#       }
#       if (!terra::compareGeom(sumr, rj, stopOnError = FALSE)) {
#         # resample just in case (should rarely happen)
#         rj <- terra::resample(rj, sumr, method = "bilinear")
#       }
#       sumr <- sumr + rj
#       ok_days <- ok_days + 1L
#       rm(rj); gc()
#     }
#
#     out_file <- file.path(out_dir, paste0("rfe2_monthsum_", m, ".tif"))
#     if (!overwrite && file.exists(out_file)) {
#       if (!quiet) message("  Exists, skipping write: ", out_file)
#     } else {
#       if (!quiet) message("  Writing: ", out_file, " [", ok_days, "/", length(idx), " days]")
#       terra::writeRaster(sumr, out_file, overwrite = TRUE)
#     }
#     out_paths[m] <- out_file
#
#     coverage$days_ok[ii]      <- ok_days
#     coverage$days_missing[ii] <- length(idx) - ok_days
#     coverage$missing_dates[[ii]] <- if (length(missed)) as.list(missed) else list(character(0))
#
#     rm(sumr, r_first); gc()
#   }
#
#   # 4) Save a simple CSV coverage (with semicolon-separated missing dates)
#   cov_csv <- coverage
#   cov_csv$missing_dates <- vapply(coverage$missing_dates, function(x) paste(x, collapse = ";"), "")
#   utils::write.csv(cov_csv, file.path(out_dir, coverage.name), row.names = FALSE)
#
#   if (!quiet) {
#     message("\n== Done ==")
#     message("Coverage summary written to: ", file.path(out_dir, "_coverage.csv"))
#   }
#
#   # return both paths and coverage for programmatic use
#   invisible(list(out_files = out_paths, coverage = coverage))
# }


# res <- rfe2_monthly_sum(
#   from = as.Date("2001-01-01"),
#   to   = as.Date("2024-12-31"),
#   out_dir = "rfe2_monthly_2001",
#   quiet = FALSE
# )
#
# res$out_files

################################################################################

# cov <- read.csv("./rfe2_monthly_2001/_coverage.csv")
# missing <- cov %>%
#   dplyr::filter(days_missing>0)
#
# all.res <- list(out_files = data.frame(),
#                 coverage = data.frame())
# for (imonth in seq(5,length(missing$month))){
#   cmonth <- missing$month[imonth]
#   cdays_expected <- missing$days_expected[imonth]
#
#   cdate.init <- as.Date(paste0(cmonth,"01"),format = "%Y%m%d")
#   cdate.end <- as.Date(paste0(cmonth,cdays_expected),format = "%Y%m%d")
#
#   res <- rfe2_monthly_sum(
#     from = cdate.init,
#     to   = cdate.end,
#     out_dir = "rfe2_monthly_2001",
#     coverage.name = "_coverage2.csv",
#     overwrite = TRUE,
#     quiet = FALSE
#   )
#
#   all.res[["out_files"]] <- bind_rows(all.res[["out_files"]],
#                                       res$out_files)
#   all.res[["coverage"]] <- bind_rows(all.res[["coverage"]],
#                                      res$coverage)
#
# }

# utils::write.csv(all.res$coverage,
#                  "rfe2_monthly_2001/_coverage2.csv",
#                  row.names = FALSE)

files <- list.files("./data/RFE2/",pattern = "*.tif",full.names = TRUE)

ym  <- sub("^.*_(\\d{6})\\.tif$", "\\1", basename(files), perl = TRUE)

A <- rast(files)
A_aligned <- project(A, "EPSG:4326")
time(A_aligned) <- as.Date(paste0(ym, "01"), "%Y%m%d")

coord <- expand.grid(lon = seq(-179.75,179.75,0.5),
                     lat = seq(-30.25,30.25,0.5)) %>%
  mutate(value = 1)
craster <- rast(raster::rasterFromXYZ(coord))
A.rspld <- resample(A_aligned,craster)
A.rspld_aligned <- project(A.rspld,
                           "EPSG:4326")

writeRaster(A.rspld,
            paste0("./outputs/RFE2_pre_all.years.tif"),
            overwrite=TRUE, gdal=c("COMPRESS=NONE", "TFW=YES"))

system2("scp",
        c(paste0("./outputs/RFE2_pre_all.years.tif"),
          "hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/all.climate/"))



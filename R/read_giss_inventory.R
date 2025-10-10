
read_giss_inventory <- function(path) {
  stopifnot(file.exists(path))
  # needs: stringr, readr, dplyr
  if (!requireNamespace("stringr", quietly = TRUE)) stop("Install 'stringr'")
  if (!requireNamespace("readr",   quietly = TRUE)) stop("Install 'readr'")
  if (!requireNamespace("dplyr",   quietly = TRUE)) stop("Install 'dplyr'")

  lines <- readr::read_lines(path)
  lines <- lines[nzchar(lines)]                      # drop empty lines
  if (length(lines) && grepl("^ID\\s+Lat\\s+Lon", lines[1])) lines <- lines[-1]  # drop header

  # ID, LAT, LON, ELEV, NAME (greedy), BI
  re <- "^(\\S+)\\s+([+-]?[0-9.]+)\\s+([+-]?[0-9.]+)\\s+([+-]?[0-9.]+)\\s+(.+?)\\s+(-?\\d+)\\s*$"
  m <- stringr::str_match(lines, re)

  # Guard against non-matching lines
  bad <- which(is.na(m[,1]))
  if (length(bad)) warning("Skipped ", length(bad), " malformed line(s).")

  dplyr::tibble(
    station_id      = m[,2],
    lat  = as.numeric(m[,3]),
    lon = as.numeric(m[,4]),
    elevation_m = as.numeric(m[,5]),
    name    = stringr::str_squish(m[,6]),
    BI      = as.integer(m[,7])
  ) |>
    dplyr::mutate(
      # common sentinel for unknown elevation
      elevation_m = dplyr::if_else(elevation_m >= 9998.9, NA_real_, elevation_m)
    )
}

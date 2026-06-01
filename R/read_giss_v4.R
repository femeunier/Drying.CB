
read_giss_v4 <- function(path,
                         na_value = -9999L,
                         scale = 1) {

  stopifnot(requireNamespace("stringi"),
            requireNamespace("data.table"))

  # read all lines (fast in base)
  L <- readLines(path, warn = FALSE)

  # Regexes (compiled by stringi)
  # Year line: "YYYY v1 v2 ... v12" (13 integers, allowing leading/trailing spaces)
  year_re <- paste0("^\\s*(\\d{4})",
                    paste(rep("\\s+(-?\\d+)", 12), collapse = ""),
                    "\\s*$")

  # Header: ID (must contain a letter to avoid matching pure years),
  #         lat lon elev, then a name (greedy, can include spaces),
  #         and a trailing integer (count)
  header_re <- "^(\\S*[A-Za-z]\\S*)\\s+([+-]?\\d+\\.?\\d*)\\s+([+-]?\\d+\\.?\\d*)\\s+([+-]?\\d+\\.?\\d*)\\s+(.+?)\\s+(\\d+)\\s*$"

  is_year   <- stringi::stri_detect_regex(L, year_re)
  is_header <- stringi::stri_detect_regex(L, header_re)

  # Extract all headers at once
  H <- stringi::stri_match_first_regex(L[is_header], header_re)
  # H columns: [1]=full, [2]=ID, [3]=lat, [4]=lon, [5]=elev, [6]=name, [7]=count
  headers <- data.table::data.table(
    line_idx   = which(is_header),
    station_id = H[, 2],
    latitude   = as.numeric(H[, 3]),
    longitude  = as.numeric(H[, 4]),
    elevation_m= as.numeric(H[, 5]),
    name       = trimws(H[, 6]),
    count      = as.integer(H[, 7])
  )

  # Extract all year lines at once
  Y <- stringi::stri_match_first_regex(L[is_year], year_re)
  # Y columns: [1]=full, [2]=year, [3:14]=12 monthly integers
  year_rows <- data.table::data.table(
    line_idx = which(is_year),
    year     = as.integer(Y[, 2]),
    m01 = as.integer(Y[, 3]),  m02 = as.integer(Y[, 4]),
    m03 = as.integer(Y[, 5]),  m04 = as.integer(Y[, 6]),
    m05 = as.integer(Y[, 7]),  m06 = as.integer(Y[, 8]),
    m07 = as.integer(Y[, 9]),  m08 = as.integer(Y[,10]),
    m09 = as.integer(Y[,11]),  m10 = as.integer(Y[,12]),
    m11 = as.integer(Y[,13]),  m12 = as.integer(Y[,14])
  )

  if (nrow(year_rows) == 0L || nrow(headers) == 0L) {
    return(data.frame(station_id=character(), name=character(), latitude=numeric(),
                      longitude=numeric(), elevation_m=numeric(),
                      year=integer(), month=integer(), value=numeric()))
  }

  # Map each year line to the most recent preceding header (by line index)
  # For efficiency: sort headers by line, and use cumulative join
  data.table::setkey(headers, line_idx)
  data.table::setkey(year_rows, line_idx)

  # Find, for each year line, the last header with line_idx <= year line
  # (rolling join backward)
  yr_with_meta <- headers[year_rows, roll = TRUE]

  # Now yr_with_meta contains header fields + the year-row columns
  # Replace NA codes and apply scaling
  month_cols <- paste0("m", sprintf("%02d", 1:12))
  for (mc in month_cols) {
    v <- yr_with_meta[[mc]]
    v[v == na_value] <- NA_integer_
    yr_with_meta[[mc]] <- as.numeric(v) * scale
  }

  # Long/tidy format: one row per month
  DT <- data.table::melt(
    yr_with_meta,
    id.vars = c("station_id", "name", "latitude", "longitude", "elevation_m", "year"),
    measure.vars = month_cols,
    variable.name = "month",
    value.name = "value",
    variable.factor = FALSE
  )

  # month "m01" -> 1..12
  DT[, month := as.integer(sub("^m", "", month))]

  # Order and return as data.frame (if you prefer data.table, drop as.data.frame)
  data.table::setorder(DT, station_id, year, month)
  as.data.frame(DT[, .(station_id, name, latitude, longitude, elevation_m, year, month, value)])
}


read_hydroweb <- function(file) {
  # Read all lines
  lines <- readLines(file, warn = FALSE)

  ## ---- METADATA ----
  meta_lines <- grep("^#", lines, value = TRUE)
  # Drop the big separator line(s)
  meta_lines <- meta_lines[!grepl("^#+\\s*$", meta_lines)]
  # Remove leading "# " or "#"
  meta_clean <- sub("^#\\s*", "", meta_lines)

  # Keep only lines with '::'
  kv_lines <- grep("::", meta_clean, value = TRUE)

  # Split into key / value on first '::'
  kv_split <- strsplit(kv_lines, "::", fixed = TRUE)
  keys <- vapply(kv_split, function(x) trimws(x[1]), character(1))
  vals <- vapply(kv_split, function(x) {
    # Recombine in case value also contains '::'
    trimws(paste(x[-1], collapse = "::"))
  }, character(1))

  # Build metadata list with clean names
  meta <- as.list(vals)
  names(meta) <- make.names(keys, unique = TRUE)

  # Try to convert numeric-like values
  meta <- lapply(meta, function(x) {
    if (identical(x, "") || is.na(x)) return(NA)
    suppressWarnings(nx <- as.numeric(x))
    if (!all(is.na(nx))) nx else x
  })

  ## ---- DATA ----
  # First non-# line = start of data
  data_start <- which(!grepl("^#", lines))[1]
  if (is.na(data_start)) {
    stop("No data block found in file.")
  }

  data_lines <- lines[data_start:length(lines)]
  data_lines <- data_lines[nzchar(trimws(data_lines))]  # drop empty lines

  # Read with whitespace sep (':' will appear as its own column)
  dat <- read.table(text = data_lines,
                    header = FALSE,
                    sep = "",
                    stringsAsFactors = FALSE,
                    na.strings = c("NA", "NaN"))

  # Remove pure ':' column if present
  colon_cols <- which(sapply(dat, function(z) all(z %in% c(":", NA, ""))))
  if (length(colon_cols) == 1L) {
    dat <- dat[, -colon_cols, drop = FALSE]
  }

  # Expect 15 columns after removing ':' (DATE + TIME + 13 others)
  expected_ncol <- 15
  if (ncol(dat) != expected_ncol) {
    stop(
      sprintf("Unexpected number of columns in data block: got %d, expected %d.",
              ncol(dat), expected_ncol)
    )
  }

  colnames(dat) <- c(
    "date",
    "time",
    "orthometric_height_m",
    "uncertainty_m",
    "lon",
    "lat",
    "ellipsoid_height_m",
    "geoid_ondulation_m",
    "distance_km",
    "satellite",
    "orbit_mission",
    "ground_track_number",
    "cycle_number",
    "retracking_algorithm",
    "gdr_version"
  )

  # Types
  num_cols <- c(
    "orthometric_height_m",
    "uncertainty_m",
    "lon",
    "lat",
    "ellipsoid_height_m",
    "geoid_ondulation_m",
    "distance_km",
    "ground_track_number",
    "cycle_number"
  )
  dat[num_cols] <- lapply(dat[num_cols], function(x) suppressWarnings(as.numeric(x)))

  # Datetime
  dat$datetime <- as.POSIXct(
    paste(dat$date, dat$time),
    format = "%Y-%m-%d %H:%M",
    tz = "UTC"
  )

  # Reorder columns: datetime first
  dat <- dat[, c("datetime", setdiff(names(dat), "datetime"))]
  dat[dat == 9999.999] <- NA
  dat[,"gdr_version"] <- as.character(dat[,"gdr_version"])

  # Return both
  list(
    metadata = meta,
    data = dat
  )
}

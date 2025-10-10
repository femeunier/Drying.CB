read_cru_ts <- function(file,
                        na_codes = c(-9999L, 9999L),
                        zero_is_na = FALSE,
                        warn_on_skip = FALSE,
                        verbose = TRUE) {
  if (!file.exists(file)) stop("File not found: ", file)
  lines <- readLines(file, warn = FALSE)
  n <- length(lines); i <- 1L

  # --- helpers ---------------------------------------------------------------
  is_year_line <- function(s) grepl("^\\s*(?:18|19|20)\\d{2}\\b", s)

  # Pre-tokenization repair: separate glued -9999 from neighbours
  # e.g., "1080-999911080" -> "1080 -9999 11080"
  unglue_neg9999 <- function(s) {
    s <- gsub("(?<=\\d)-9999(?=\\d)", " -9999 ", s, perl = TRUE)
    s <- gsub("(?<=\\d)-9999\\b",       " -9999",   s, perl = TRUE)
    s <- gsub("\\b-9999(?=\\d)",        "-9999 ",   s, perl = TRUE)
    s
  }

  # Extract signed integer tokens (after -9999 unglue)
  ints_chr <- function(s) regmatches(s, gregexpr("-?\\d+", s, perl = TRUE))[[1]]

  # Split an unsigned long run into i5-ish chunks, right-biased (prefer 5-digit tails).
  # Returns a character vector of chunks.
  split_to_i5 <- function(tok) {
    if (!grepl("^[0-9]{6,}$", tok)) return(tok)
    L <- nchar(tok); out <- character(0)
    pos <- 1L
    # We want tail chunks of 5 digits; compute how many leading digits to keep before the 5s
    k5 <- L %/% 5; rem <- L %% 5
    lead <- if (rem == 0) 5 else rem  # keep a first chunk of 'rem' (1..4), else 5 when divisible
    # But avoid silly 1–2 digit leading pieces if we don't need them:
    if (lead <= 2 && k5 >= 1) { lead <- lead + 5; k5 <- k5 - 1 }
    # first chunk
    if (lead > 0 && lead < L) {
      out <- c(out, substr(tok, pos, pos + lead - 1)); pos <- pos + lead
    }
    # then 5-digit chunks
    while (pos <= L) {
      take <- min(5L, L - pos + 1L)
      out <- c(out, substr(tok, pos, pos + take - 1)); pos <- pos + take
    }
    out
  }

  # Ensure exactly 12 month tokens (character → integer), by splitting long runs if needed.
  # Returns integer(12) or NULL if impossible.
  to_12_months <- function(month_chr) {
    # First pass: split any long unsigned runs
    expanded <- unlist(lapply(month_chr, split_to_i5), use.names = FALSE)
    # If still under 12, keep splitting the longest unsigned piece until we reach 12 or can’t.
    while (length(expanded) < 12) {
      lens <- nchar(expanded)
      # pick longest unsigned numeric piece with len >= 6
      cand_idx <- which(grepl("^[0-9]{6,}$", expanded))
      if (!length(cand_idx)) break
      j <- cand_idx[ which.max(lens[cand_idx]) ]
      repl <- split_to_i5(expanded[j])
      if (length(repl) == 1L) break
      expanded <- c(expanded[seq_len(j-1)], repl, expanded[seq.int(j+1, length(expanded))])
    }
    # Trim or pad to 12
    if (length(expanded) > 12) expanded <- expanded[1:12]
    if (length(expanded) < 12) expanded <- c(expanded, rep("NA", 12 - length(expanded)))

    mm <- suppressWarnings(as.integer(expanded))
    mm[mm %in% na_codes] <- NA_integer_
    if (length(mm) != 12) return(NULL)
    mm
  }

  # Header per spec: WMO LAT LON ALT STATION COUNTRY START END (single spaces; middle has letters)
  hdr_rx <- paste0(
    "^\\s*(\\d{1,7})\\s+([+-]?\\d{1,5})\\s+([+-]?\\d{1,6})\\s+([+-]?\\d{1,4})\\s+",
    "([^\\d]*[A-Za-z][A-Za-z0-9 .,'/\\-]*?)\\s+",
    "((?:18|19|20)\\d{2})\\s+((?:18|19|20)\\d{2})\\s*$"
  )
  parse_header <- function(s) {
    hdr_rx <- paste0(
      "^\\s*(\\d{1,7})\\s+([+-]?\\d{1,5})\\s+([+-]?\\d{1,6})\\s+([+-]?\\d{1,4})\\s+",
      "([^\\d]*[A-Za-z][A-Za-z0-9 .,'/\\-]*?)\\s+",
      "((?:18|19|20)\\d{2})\\s+((?:18|19|20)\\d{2})\\s*$"
    )
    m <- regexec(hdr_rx, s, perl = TRUE); g <- regmatches(s, m)[[1]]
    if (length(g) != 8) return(NULL)

    wmo     <- as.integer(g[2])
    lat_raw <- as.integer(g[3])   # ×100
    lon_raw <- as.integer(g[4])   # ×100
    alt_raw <- as.integer(g[5])
    mid     <- trimws(g[6])
    start_y <- as.integer(g[7]); end_y <- as.integer(g[8])

    # Station / country split (country may be absent)
    parts <- strsplit(mid, " {2,}", perl = TRUE)[[1]]
    if (length(parts) >= 2) {
      cand_country <- trimws(parts[length(parts)])
      cand_name    <- trimws(paste(parts[1:(length(parts)-1)], collapse = "  "))
      if (nchar(cand_country) <= 13 && grepl("^[A-Z][A-Z .'-]*$", cand_country)) {
        station_name <- cand_name; country <- cand_country
      } else {
        station_name <- mid; country <- ""
      }
    } else {
      station_name <- mid; country <- ""
    }

    # ---- MISSING-CODE AWARE SCALING ----
    # Common CRU sentinels
    na_lat <- c(-9999L,  9999L, 32767L)
    na_lon <- c(-99999L, 99999L, 32767L)
    na_alt <- c(-999L,    9999L, 32767L)

    lat <- if (is.na(lat_raw) || lat_raw %in% na_lat || abs(lat_raw) >  9000L) NA_real_ else lat_raw / 100
    lon <- if (is.na(lon_raw) || lon_raw %in% na_lon || abs(lon_raw) > 18000L) NA_real_ else {
      x <- lon_raw / 100
      if (x > 180) x <- x - 360
      if (x < -180) x <- x + 360
      x
    }
    elev_m <- if (is.na(alt_raw) || alt_raw %in% na_alt) NA_integer_ else alt_raw

    list(
      wmo = wmo, station_name = station_name, country = country,
      lat = lat, lon = lon, elev_m = elev_m,
      start_year = start_y, end_year = end_y
    )
  }


  out <- list(); n_st <- 0L
  while (i <= n) {
    # find next header
    hdr <- NULL
    while (i <= n && is.null(hdr)) { hdr <- parse_header(lines[i]); if (is.null(hdr)) i <- i + 1L }
    if (is.null(hdr)) break
    i <- i + 1L

    # skip exactly one normals line (if present)
    if (i <= n && !is_year_line(lines[i])) i <- i + 1L

    years <- integer(0); rows <- list()

    while (i <= n) {
      # stop on next header
      if (!is.null(parse_header(lines[i]))) break
      if (!is_year_line(lines[i])) { i <- i + 1L; next }

      # repair "-9999" glue then tokenize
      L <- unglue_neg9999(lines[i])
      toks <- ints_chr(L); i <- i + 1L
      if (!length(toks)) next

      yr <- suppressWarnings(as.integer(toks[1]))
      mm <- to_12_months(toks[-1])
      if (is.null(mm)) {
        if (warn_on_skip) warning("Skipping year ", yr, " (could not form 12 months)")
        next
      }
      if (zero_is_na) mm[mm == 0L] <- NA_integer_

      years <- c(years, yr)
      rows[[length(rows)+1L]] <- mm
    }

    if (length(rows)) {
      M <- do.call(rbind, rows)
      out[[length(out)+1L]] <- data.frame(
        wmo          = rep.int(hdr$wmo,          nrow(M) * 12),
        station_name = rep.int(hdr$station_name, nrow(M) * 12),
        country      = rep.int(hdr$country,      nrow(M) * 12),
        lat          = rep.int(hdr$lat,          nrow(M) * 12),
        lon          = rep.int(hdr$lon,          nrow(M) * 12),
        elev_m       = rep.int(hdr$elev_m,       nrow(M) * 12),
        start_year   = rep.int(hdr$start_year,   nrow(M) * 12),
        end_year     = rep.int(hdr$end_year,     nrow(M) * 12),
        year         = rep(years, each = 12),
        month        = rep(1:12,  times = length(years)),
        value        = as.numeric(t(M)),
        stringsAsFactors = FALSE
      )
      n_st <- n_st + 1L
    }
  }

  if (verbose) message("Stations parsed: ", n_st)
  if (!length(out)) {
    return(data.frame(
      wmo=integer(), station_name=character(), country=character(),
      lat=numeric(), lon=numeric(), elev_m=integer(),
      start_year=integer(), end_year=integer(),
      year=integer(), month=integer(), value=numeric(), stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, out)
}

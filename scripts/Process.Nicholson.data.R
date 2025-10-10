library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(stringi)
library(ggplot2)

# Robust Nicholson multi-station parser (variable headers, sentinels, 9999/99999 as NA)
# Output:
#   $metadata    : data.table(station_id, station_name, lat_dd, lon_dd, elev_m, source, var_code)
#   $climatology : data.table(station_id, month, mean_mm, sd_mm, annual_mean, annual_sd)
#   $monthly     : data.table(station_id, year, month, value_mm, annual_total)

library(data.table)
library(stringi)

read_nicholson_all <- function(file,
                               na_month_code = 9999,
                               na_annual_code = 99999,
                               progress = FALSE) {

  # --- dependencies (explicit namespaces used below) ---
  if (!requireNamespace("stringi", quietly = TRUE))
    stop("Package 'stringi' is required.")
  if (!requireNamespace("data.table", quietly = TRUE))
    stop("Package 'data.table' is required.")

  lines <- readLines(file, warn = FALSE)
  n <- length(lines)
  if (!n) stop("Empty file.")

  # index helpers
  mean_idx <- which(stringi::stri_detect_regex(lines, "\\sMEAN(\\s|$)"))
  stdv_idx <- which(stringi::stri_detect_regex(lines, "\\sSTDV(\\s|$)"))
  if (length(mean_idx) == 0L || length(mean_idx) != length(stdv_idx))
    stop("Could not match MEAN and STDV lines consistently.")
  header_idx <- mean_idx - 1L
  k <- length(header_idx)

  # Extract exact station_id from the MEAN line (safe even when it contains spaces)
  station_id_vec <- stringi::stri_trim_right(
    stringi::stri_replace_first_regex(lines[mean_idx], "\\sMEAN.*$", "")
  )

  # ------------------- parse_header_africa (self-contained) -------------------
  parse_header_africa <- local({
    .AF_LAT_MIN <- -36; .AF_LAT_MAX <- 38
    .AF_LON_MIN <- -26; .AF_LON_MAX <- 60
    .NA_CODES_RAW   <- c(-999, 999, 9999, 32767)
    .NA_SCALED_LAT  <- c(-9.99, 9.99, -99.9, 99.9)
    .NA_SCALED_LON  <- c(-9.99, 9.99, -99.9, 99.9)

    is_int <- function(x) grepl("^-?\\d+$", x)

    africa_pick_lat <- function(lat_raw) {
      if (is.na(lat_raw) || lat_raw %in% .NA_CODES_RAW) return(NA_real_)
      cand <- c(lat_raw/100, lat_raw/10, lat_raw/1)
      cand <- cand[is.finite(cand) & !(round(cand, 2) %in% .NA_SCALED_LAT)]
      inside <- cand[cand >= .AF_LAT_MIN & cand <= .AF_LAT_MAX]
      if (length(inside)) inside[1] else if (length(cand)) cand[which.min(abs(cand))] else NA_real_
    }

    africa_pick_lon <- function(lon_raw) {
      if (is.na(lon_raw) || lon_raw %in% .NA_CODES_RAW) return(NA_real_)
      cand <- c(lon_raw/100, lon_raw/10, lon_raw/1)
      # wrap 0–360 → [-180,180] when needed
      cand <- c(cand, ifelse(cand > 180, cand - 360, cand))
      cand <- cand[is.finite(cand) & !(round(cand, 2) %in% .NA_SCALED_LON)]
      inside <- cand[cand >= .AF_LON_MIN & cand <= .AF_LON_MAX]
      if (length(inside)) inside[1] else if (length(cand)) cand[which.min(abs(cand))] else NA_real_
    }

    function(hline, stid, elev_na_codes = c(-999, 999, 9999, 32767)) {
      line <- trimws(hline)

      # Detect glued lat±lon anywhere (allow spaces around +/-)
      m <- regexpr("(?<![A-Za-z0-9])(\\d{2,5})\\s*([\\+\\-])\\s*(\\d{2,5})(?![A-Za-z0-9])",
                   line, perl = TRUE)
      lat_raw <- lon_raw <- NA_integer_
      line_wo_glue <- line
      if (m[1] != -1) {
        gl <- regmatches(line, m)[1]
        parts <- regmatches(gl, gregexpr("\\d{2,5}", gl, perl = TRUE))[[1]]
        signc <- regmatches(gl, regexpr("[\\+\\-]", gl, perl = TRUE))
        lat_raw <- suppressWarnings(as.integer(parts[1]))
        lon_raw <- suppressWarnings(as.integer(parts[2])) * ifelse(signc == "-", -1L, 1L)
        line_wo_glue <- sub("(\\d{2,5})\\s*([\\+\\-])\\s*(\\d{2,5})", " ", line, perl = TRUE)
      }

      # Tokenize (after removing glue)
      tok <- strsplit(trimws(line_wo_glue), "\\s+", perl = TRUE)[[1]]
      tok <- tok[tok != ""]
      ints_idx <- which(is_int(tok))
      ints     <- suppressWarnings(as.integer(tok[ints_idx]))
      texts    <- tok[setdiff(seq_along(tok), ints_idx)]

      # var_code = last "small" integer (≤ 200) on the line
      var_code <- NA_integer_; var_pos <- NA_integer_
      if (length(ints)) {
        small <- which(ints <= 200)
        if (length(small)) {
          var_code <- ints[tail(small, 1)]
          var_pos  <- ints_idx[tail(small, 1)]
        }
      }

      # elevation = nearest plausible integer *left of var_pos*
      elev_m <- NA_integer_; elev_pos <- NA_integer_
      pick_plaus_elev <- function(idx_vec) {
        for (p in rev(idx_vec)) {
          v <- suppressWarnings(as.integer(tok[p]))
          if (!is.na(v) && v %in% elev_na_codes) { elev_pos <<- p; elev_m <<- NA_integer_; return(invisible()) }
          if (!is.na(v) && v >= -100 && v <= 9000) { elev_pos <<- p; elev_m <<- v; return(invisible()) }
        }
        invisible()
      }
      if (!is.na(var_pos)) {
        left <- ints_idx[ints_idx < var_pos]
        if (length(left)) pick_plaus_elev(left)
      }
      if (is.na(elev_pos) && length(ints_idx)) pick_plaus_elev(ints_idx)

      # If lat/lon not provided via glue, derive from integers
      if (is.na(lat_raw) || is.na(lon_raw)) {
        boundary <- if (!is.na(elev_pos)) elev_pos else if (!is.na(var_pos)) var_pos else Inf
        cand_idx <- ints_idx[ints_idx < boundary]
        if (!length(cand_idx)) cand_idx <- ints_idx

        cand_vals <- suppressWarnings(as.integer(tok[cand_idx]))
        cand_vals <- cand_vals[!(cand_vals %in% .NA_CODES_RAW) & is.finite(cand_vals)]

        if (length(cand_vals) < 2 && length(ints)) {
          base_vals <- ints[!(ints %in% .NA_CODES_RAW)]
          if (length(base_vals) >= 2) cand_vals <- base_vals[1:2]
        }

        if (length(cand_vals) >= 2) {
          a <- cand_vals[1]; b <- cand_vals[2]
          lat1 <- africa_pick_lat(a); lon1 <- africa_pick_lon(b)
          lat2 <- africa_pick_lat(b); lon2 <- africa_pick_lon(a)

          pick1_ok <- is.finite(lat1) && is.finite(lon1) &&
            lat1 >= .AF_LAT_MIN && lat1 <= .AF_LAT_MAX &&
            lon1 >= .AF_LON_MIN && lon1 <= .AF_LON_MAX
          pick2_ok <- is.finite(lat2) && is.finite(lon2) &&
            lat2 >= .AF_LAT_MIN && lat2 <= .AF_LAT_MAX &&
            lon2 >= .AF_LON_MIN && lon2 <= .AF_LON_MAX

          if (pick1_ok) { lat_raw <- a; lon_raw <- b }
          else if (pick2_ok) { lat_raw <- b; lon_raw <- a }
          else {
            bad1 <- sum(is.na(c(lat1, lon1))) + sum(round(c(lat1, lon1), 2) %in% c(.NA_SCALED_LAT, .NA_SCALED_LON))
            bad2 <- sum(is.na(c(lat2, lon2))) + sum(round(c(lat2, lon2), 2) %in% c(.NA_SCALED_LAT, .NA_SCALED_LON))
            if (bad1 <= bad2) { lat_raw <- a; lon_raw <- b } else { lat_raw <- b; lon_raw <- a }
          }
        } else {
          stop(sprintf("Could not locate lat/lon in header: '%s'", hline))
        }
      }

      # Final scale to degrees
      lat_dd <- africa_pick_lat(lat_raw)
      lon_dd <- africa_pick_lon(lon_raw)

      # Station name = tokens before the first integer (minus station_id)
      toks_h <- strsplit(trimws(hline), "\\s+", perl = TRUE)[[1]]
      first_int_pos <- which(is_int(toks_h))[1]
      lead_text <- if (!is.na(first_int_pos)) paste(toks_h[seq_len(first_int_pos - 1)], collapse = " ") else paste(toks_h, collapse = " ")
      station_name <- trimws(sub(paste0("^", stid), "", lead_text, fixed = TRUE))

      # Source: pick the last plausible text token
      source <- ""
      if (length(texts)) {
        patt <- "(?i)^(GHCN\\s*V?\\d*|GHCNV\\d+|AGRHYMET|CANAR|SUDAN|SOSUD|SWAZILAND|VF)$"
        hits <- grep(patt, texts, perl = TRUE)
        if (length(hits)) {
          source <- texts[hits[length(hits)]]
        } else {
          caps <- grep("^[A-Z0-9][A-Z0-9]+$", texts)
          if (length(caps)) source <- texts[caps[length(caps)]]
        }
      }

      list(
        station_name = station_name,
        lat = lat_dd,
        lon = lon_dd,
        elev = if (length(elev_m)) elev_m else NA_integer_,
        source = source,
        var_code = if (length(var_code)) var_code else NA_integer_
      )
    }
  })
  # ----------------- end parse_header_africa -----------------

  # Build metadata table
  meta_list <- vector("list", k)
  for (i in seq_len(k)) {
    stid <- station_id_vec[i]
    hdr  <- parse_header_africa(lines[header_idx[i]], stid)
    meta_list[[i]] <- data.table::data.table(
      station_id   = stid,
      station_name = hdr$station_name,
      lat_dd       = hdr$lat,
      lon_dd       = hdr$lon,
      elev_m       = hdr$elev,
      source       = hdr$source,
      var_code     = hdr$var_code
    )
  }
  metadata <- data.table::rbindlist(meta_list, use.names = TRUE)

  # Climatology (MEAN/STDV): take first 13 integers after the keyword
  num_extract <- function(s) {
    out <- stringi::stri_extract_all_regex(s, "-?\\d+", simplify = TRUE)
    as.numeric(out[!is.na(out)])
  }

  clim_list <- vector("list", k)
  for (i in seq_len(k)) {
    mv <- num_extract(lines[mean_idx[i]])
    sv <- num_extract(lines[stdv_idx[i]])
    if (length(mv) < 13L || length(sv) < 13L)
      stop(sprintf("Climatology line lacks 13 values for station_id '%s'", station_id_vec[i]))
    mv <- mv[1:13]; sv <- sv[1:13]
    mv[13] <- ifelse(mv[13] == na_annual_code, NA_real_, mv[13])
    sv[13] <- ifelse(sv[13] == na_annual_code, NA_real_, sv[13])
    clim_list[[i]] <- data.table::data.table(
      station_id   = station_id_vec[i],
      month        = 1:12,
      mean_mm      = mv[1:12],
      sd_mm        = sv[1:12],
      annual_mean  = mv[13],
      annual_sd    = sv[13]
    )
  }
  climatology <- data.table::rbindlist(clim_list, use.names = TRUE)

  # Yearly blocks (fixed for glued station_id+year)
  block_start <- stdv_idx + 1L
  block_end   <- c(header_idx[-1] - 1L, n)
  monthly_list <- vector("list", k)

  if (progress) {
    pb <- utils::txtProgressBar(min = 0, max = k, style = 3)
    on.exit(close(pb), add = TRUE)
  }

  for (i in seq_len(k)) {
    if (progress) utils::setTxtProgressBar(pb, i)
    stid <- station_id_vec[i]
    s <- block_start[i]; e <- block_end[i]
    if (s > e) { monthly_list[[i]] <- data.table::data.table(); next }

    blk <- lines[s:e]

    # Cut at sentinel: line starts with exact station_id then only 9s
    trim_blk <- stringi::stri_trim_both(blk)
    stid_esc <- gsub("([\\^$.|?*+(){}\\[\\]\\\\])", "\\\\\\1", stid)
    sent_pos <- which(stringi::stri_startswith_fixed(trim_blk, stid) &
                        grepl(paste0("^", stid_esc, "9+$"), trim_blk))
    if (length(sent_pos)) {
      cut_at <- min(sent_pos) - 1L
      blk <- if (cut_at >= 1L) blk[seq_len(cut_at)] else character(0)
    }
    if (!length(blk)) { monthly_list[[i]] <- data.table::data.table(); next }

    years <- integer(0)
    vals  <- list()

    for (L in blk) {
      Ltrim <- trimws(L)

      # Prefer: year glued right after station id (with or without a space)
      m_glued <- regexec(paste0("^", stid_esc, "\\s*(\\d{4})\\b"), Ltrim, perl = TRUE)
      reg <- regmatches(Ltrim, m_glued)[[1]]

      yr <- NA_integer_
      rest <- Ltrim
      if (length(reg)) {
        yr <- as.integer(reg[2])
        rest <- sub(paste0("^", stid_esc, "\\s*\\d{4}\\b"), "", Ltrim, perl = TRUE)
      } else {
        # Fallback: generic scan of integers and pick first 1800–2099
        all_nums <- stringi::stri_extract_all_regex(Ltrim, "-?\\d+", simplify = FALSE)[[1]]
        if (length(all_nums)) {
          nums <- suppressWarnings(as.integer(all_nums))
          yi <- which(nums >= 1800 & nums <= 2099)[1]
          if (!is.na(yi)) {
            yr <- nums[yi]
            rest <- sub("(.*?\\b(?:18|19|20)\\d{2}\\b)", "", Ltrim, perl = TRUE)
          }
        }
      }

      if (is.na(yr)) next  # cannot find a year on this line

      # Extract up to 13 integers after the year: 12 months + annual
      after_nums <- stringi::stri_extract_all_regex(rest, "-?\\d+", simplify = FALSE)[[1]]
      if (!length(after_nums)) next
      an <- suppressWarnings(as.integer(after_nums))
      an <- an[is.finite(an)]
      if (length(an) < 13L) an <- c(an, rep(NA_integer_, 13L - length(an)))
      if (length(an) > 13L) an <- an[1:13]

      # Apply Nicholson missing codes
      m12 <- an[1:12]; m12[m12 == na_month_code] <- NA_integer_
      ann <- an[13];   if (!is.na(ann) && ann == na_annual_code) ann <- NA_integer_

      years <- c(years, yr)
      vals[[length(vals) + 1L]] <- c(m12, ann)
    }

    if (!length(vals)) { monthly_list[[i]] <- data.table::data.table(); next }

    V <- matrix(unlist(vals, use.names = FALSE), ncol = 13, byrow = TRUE)
    monthly_list[[i]] <- data.table::data.table(
      station_id    = stid,
      year          = rep(years, each = 12),
      month         = rep(1:12, times = length(years)),
      value_mm      = as.numeric(t(V[, 1:12, drop = FALSE])),
      annual_total  = rep(as.numeric(V[, 13]), each = 12)
    )
  }

  monthly <- data.table::rbindlist(monthly_list, use.names = TRUE, fill = TRUE)
  data.table::setkey(monthly, station_id, year, month)
  data.table::setkey(climatology, station_id, month)

  list(
    metadata    = metadata[],
    climatology = climatology[],
    monthly     = monthly[]
  )
}


res <- read_nicholson_all("./data/station_data/Nicholson/NIC131_03_2021.txt",
                          progress = TRUE)

res$monthly %>%
  dplyr::filter(station_id == "GUINBBISS2")

metadata <- (res$metadata) %>%
  rename(lon = lon_dd,
         lat = lat_dd)

sites.selected <- metadata %>%
  filter(lat <= 25, lat>=-25,
         lon >= -25, lon <= 65) %>%
  pull(station_id)


world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

raw.monthly <- as.data.frame(res$monthly) %>%
  left_join(metadata %>%
              dplyr::select(station_id,lon,lat,elev_m) %>%
              rename(elevation = elev_m) %>%
              distinct(),
            by = c("station_id"))

raw.sum <- raw.monthly %>%
  group_by(station_id) %>%
  summarise(N = n(),
            year.min = min(year),
            year.max = max(year),
            .groups = "keep")


ggplot() +
  stat_bin2d(data = raw.monthly, alpha = 0.8,
             aes(lon, lat, fill = stat(log10(count))),
             binwidth = c(1, 1)) +
  labs(fill = "Obs count", x = "Longitude", y = "Latitude") +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  scale_fill_gradient2(mid = "grey") +
  scale_x_continuous(limits = range(raw.monthly$lon,
                                    na.rm = TRUE) + c(-1, 1)) +
  scale_y_continuous(limits = range(raw.monthly$lat,
                                    na.rm = TRUE) + c(-1, 1)) +
  theme_minimal()


N.year <- raw.monthly %>%
  dplyr::select(station_id,year) %>%
  distinct() %>%
  group_by(year) %>%
  summarise(N = n()) %>%
  ungroup() %>%
  complete(year = c(min(year):max(year)),
           fill = list(N = 0))

ggplot(data = N.year) +
  geom_line(aes(x = year, y = N)) +
  theme_bw()

MAP <- raw.monthly %>%
  filter(abs(lat) <= 25) %>%
  group_by(station_id, year) %>%
  summarise(N = n(),
            lat = unique(lat),
            lon = unique(lon),
            MAP = sum(value_mm,na.rm = TRUE),
            .groups = 'keep')

ggplot() +
  guides(fill = "none") +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  geom_point(data = MAP %>%
               dplyr::filter(year %in% c(1960:1985)) %>%
               ungroup() %>%
               dplyr::select(station_id,lon,lat) %>%
               distinct(),
             aes(x = lon, y = lat), size = 0.6, alpha = 0.5) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  theme_bw()

MAP.sum <- MAP %>%
  group_by(year) %>%
  summarise(MAP.m = mean(MAP,na.rm = TRUE),
            N = n(),
            .groups = "keep") %>%
  dplyr::filter(year %in% c(1960:1985))

scale_factor <- max(MAP.sum$MAP.m, na.rm = TRUE) / max(MAP.sum$N, na.rm = TRUE)

ggplot(MAP.sum,
       aes(x = year)) +
  geom_line(aes(y = MAP.m, color = "MAP.m"), size = 1) +
  stat_smooth(aes(y = MAP.m),
              method = "lm", color = "grey17") +
  geom_line(aes(y = N * scale_factor, color = "N"), size = 1) +
  scale_y_continuous(
    name = "MAP.m",
    sec.axis = sec_axis(~ . / scale_factor, name = "N")
  ) +
  # scale_color_manual(values = c("MAP.m" = "blue", "N" = "red")) +
  labs(x = "Year", color = "") +
  theme_bw()

ggplot() +
  stat_bin2d(data = raw.monthly %>%
               dplyr::select(station_id,lon,lat) %>%
               distinct(), alpha = 0.8,
             aes(lon, lat),
             binwidth = c(1, 1)) +
  labs(fill = "Obs count", x = "Longitude", y = "Latitude") +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  scale_fill_gradient2(mid = "grey") +
  scale_x_continuous(limits = range(raw.monthly$lon,
                                    na.rm = TRUE) + c(-1, 1)) +
  scale_y_continuous(limits = range(raw.monthly$lat,
                                    na.rm = TRUE) + c(-1, 1)) +
  theme_minimal()


monthly.data <- raw.monthly  %>%
  filter(station_id %in% c(sites.selected)) %>%
  filter(!is.na(value_mm))

monthly.data.selected <- monthly.data %>%
  pull(station_id)

ggplot() +
  guides(fill = "none") +
  geom_sf(data = world, fill = NA, color = "grey30", linewidth = 0.2) +
  geom_point(data = metadata,
             aes(x = lon, y = lat), size = 0.6, alpha = 0.5) +
  coord_sf(xlim = c(-25, 65), ylim = c(-25, 25), expand = FALSE) +
  theme_bw()

saveRDS(monthly.data %>%
          filter(station_id %in% monthly.data.selected) %>%
          dplyr::select(-c(annual_total,lon,lat)) %>%
          rename(pre = value_mm),
        "./data/station_data/Nicholson/data.Nicholson.RDS")


raw.monthly %>%
  filter(station_id == "GUINBBISS2")

MD <- monthly.data %>%
  dplyr::select(c(station_id,lon,lat,elevation)) %>%
  distinct()


saveRDS(MD,
        "./data/station_data/Nicholson/metadata.Nicholson.RDS")


hist(monthly.data$elevation)

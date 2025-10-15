rm(list = ls())

library(dplyr)
library(ggplot2)
library(tidyr)
library(units)

all.metadata <-
  bind_rows(
    readRDS("./data/station_data/Nicholson/metadata.Nicholson.RDS") %>%
      mutate(source = 'Nicholson'),
    readRDS("./data/station_data/CRU/metadata.CRU.RDS") %>%
      mutate(station_id = as.character(station_id)) %>%
      mutate(source = 'CRU'),
    readRDS("./data/station_data/Castellanos/metadata.Castellanos.RDS") %>%
      mutate(source = 'Castellanos'),
    readRDS("./data/station_data/GHCN/metadata.GHCN.RDS") %>%
      mutate(source = 'GHCN'),
    readRDS("./data/station_data/Meteostat/metadata.Meteostat.RDS") %>%
      mutate(source = 'Meteostat') %>%
      mutate(station_id = as.character(station_id)),
    readRDS("./data/station_data/TAHMO/metadata.TAHMO.RDS") %>%
      mutate(source = 'TAHMO'),
    readRDS("./data/station_data/GISS/metadata.GISS.RDS") %>%
      mutate(source = "GISS"),
    readRDS("./data/station_data/GSOD/metadata.GSOD.RDS") %>%
      mutate(source = 'GSOD') %>%
      mutate(station_id = as.character(station_id))) %>%
  mutate(station_id = case_when(station_id == "0" ~ NA_character_,
                                TRUE ~ station_id))


######################################################################

all.data <-
  bind_rows(
    readRDS("./data/station_data/Nicholson/data.Nicholson.RDS") %>%
      mutate(source = 'Nicholson'),
    readRDS("./data/station_data/CRU/data.CRU.RDS") %>%
      mutate(station_id = as.character(station_id)) %>%
      mutate(source = 'CRU'),
    readRDS("./data/station_data/Castellanos/data.Castellanos.RDS") %>%
      mutate(source = 'Castellanos'),
    readRDS("./data/station_data/GHCN/data.GHCN.RDS") %>%
      mutate(source = 'GHCN'),
    readRDS("./data/station_data/Meteostat/data.Meteostat.RDS") %>%
      mutate(source = 'Meteostat') %>%
      mutate(station_id = as.character(station_id)),
    readRDS("./data/station_data/TAHMO/data.TAHMO.RDS") %>%
      mutate(source = 'TAHMO'),
    readRDS("./data/station_data/GISS/data.GISS.RDS") %>%
      mutate(source = "GISS"),
    readRDS("./data/station_data/GSOD/data.GSOD.RDS") %>%
      mutate(source = 'GSOD') %>%
      mutate(station_id = as.character(station_id))) %>%
  filter(year <= 2024) %>%
  ungroup() %>%
  mutate(pre = case_when(pre > 2000 ~ NA_real_,
                         TRUE ~ pre),
         tasmin = case_when(tasmin %in% c(-99,-9999,99) ~ NA_real_,
                            TRUE ~ tasmin),
         tas = case_when(tas %in% c(-99,-9999,99) ~ NA_real_,
                         tasmax > 60 ~ NA_real_,
                         tas < -39 ~ NA_real_,
                         TRUE ~ tas),
         tasmax = case_when(tasmax %in% c(-99,-9999,99) ~ NA_real_,
                            tasmax > 200 ~ tasmax/10,
                            tasmax > 60 ~ NA_real_,
                            TRUE ~ tasmax))




# Got it — the safest way to “remove likely duplicates but keep all data” is to (1) cluster station metadata that refer to the same physical site (by distance + elevation tolerance, and guarding against same station_id used for different places), then (2) attach a stable unified_id to every row of your time-series table via a crosswalk instead of dropping anything.
# Here’s a ready-to-use R function that does exactly that. It:
#   builds proximity links with great-circle distance (meters) and optional elevation checks,
# handles missing elevations gracefully,
# clusters with graph connected components,
# produces: (a) a canonical station table, (b) a crosswalk from each source record to a unified_id, and (c) your data with unified_id added.
# Notes / choices
# Same station_id ≠ auto-duplicate: we still require proximity within max_dist_m (and elevation if available), so two different sites reusing the same ID won’t be merged.
# Missing elevation: if either record lacks elevation, we rely on distance alone; if both have elevation, we also require |Δelev| ≤ max_elev_diff.
# No data loss: the function never drops time-series rows; it only annotates them with unified_id. You can later aggregate by unified_id if you want a single series per physical site.
# Priority source: set prefer_source_order to pick which source defines the canonical name/coords when multiple overlap.

dedupe_weather_stations <- function(
    meta,
    data,
    max_dist_m = 1000,
    max_elev_diff = 50,
    prefer_source_order = NULL,
    # NEW: handle duplicates inside the same source
    within_source_tol_m = 200,          # treat coords within 200 m as identical
    within_source_elev_diff = 30,       # and |Δelev| ≤ 30 m
    on_within_source_conflict = c("error","warn_pick","pick_first") # what if too far apart?
){
  on_within_source_conflict <- match.arg(on_within_source_conflict)
  stopifnot(all(c("station_id","lon","lat","source") %in% names(meta)))
  if (!"elevation" %in% names(meta)) meta$elevation <- NA_real_

  suppressPackageStartupMessages({library(sf); library(dplyr); library(igraph)})

  meta <- meta %>%
    filter(is.finite(lon), is.finite(lat)) %>%
    mutate(
      row_id    = row_number(),
      id_key    = paste0(source, "::", station_id),
      elevation = suppressWarnings(as.numeric(elevation))
    )

  # ---------- 0) DEDUPE *WITHIN* EACH SOURCE ----------
  # One row per (source, station_id). If multiple:
  # - If all members are within 'within_source_tol_m' and 'within_source_elev_diff',
  #   collapse to a canonical record (prefer non-NA elevation, then earliest row_id).
  # - Else: depending on 'on_within_source_conflict' -> error / warn and pick best / pick first.
  sf::sf_use_s2(TRUE)
  pts0 <- st_as_sf(meta, coords = c("lon","lat"), crs = 4326, remove = FALSE)

  dup_groups <- meta %>% count(source, station_id, name = "n") %>% filter(!is.na(station_id), n > 1L)
  handled <- list()

  if (nrow(dup_groups)) {
    for (i in seq_len(nrow(dup_groups))) {
      g <- dup_groups[i,]
      idx <- which(meta$source == g$source & meta$station_id == g$station_id)
      # distances within the group
      dmat <- as.matrix(st_distance(pts0[idx, ], pts0[idx, ])) # meters
      elev <- meta$elevation[idx]
      # tolerances

      tol_m <- set_units(within_source_tol_m, "m")

      dmat  <- sf::st_distance(pts0[idx, ], pts0[idx, ])         # units "m"
      max_d <- max(dmat, na.rm = TRUE)                           # still units "m"

      elev  <- meta$elevation[idx]
      de_mat <- abs(outer(elev, elev, "-"))
      max_de <- if (all(is.na(de_mat))) NA_real_ else max(de_mat, na.rm = TRUE)

      close_enough <- (is.finite(as.numeric(max_d)) && max_d <= tol_m) &&
        (is.na(max_de) || max_de <= within_source_elev_diff)




      if (!close_enough && on_within_source_conflict == "error") {
        stop(sprintf("Within-source conflict for %s::%s (max Δdist=%.1fm, Δelev=%s).",
                     g$source, g$station_id, max_d,
                     if (is.finite(max_de)) sprintf("%.1fm", max_de) else "NA"))
      }
      # pick canonical row
      pick <- meta[idx, ] %>%
        arrange(is.na(elevation), row_id) %>%
        slice(1)
      # replace all members by the canonical values (so (source,station_id) becomes unique)
      meta[idx, c("lon","lat","elevation")] <- pick[ , c("lon","lat","elevation")]
      handled[[length(handled)+1]] <- tibble::tibble(
        source = g$source, station_id = g$station_id,
        n_members = g$n, max_dist_m = as.numeric(max_d),
        max_elev_diff = if (is.finite(max_de)) as.numeric(max_de) else NA_real_,
        action = if (close_enough) "collapsed" else
          switch(on_within_source_conflict,
                 warn_pick = "warn_pick", pick_first = "pick_first", error = "error")
      )
      if (!close_enough && on_within_source_conflict == "warn_pick") {
        warning(sprintf("Within-source conflict for %s::%s (Δdist=%.1fm, Δelev=%s). Picking first.",
                        g$source, g$station_id, max_d,
                        if (is.finite(max_de)) sprintf("%.1fm", max_de) else "NA"))
      }
    }
  }
  within_source_report <- dplyr::bind_rows(handled)

  # Now enforce uniqueness explicitly
  meta <- meta %>%
    arrange(source, station_id, is.na(elevation), row_id) %>%
    distinct(source, station_id, .keep_all = TRUE)

  # ---------- 1) SPATIAL CLUSTERING *ACROSS* SOURCES ----------
  pts <- st_as_sf(meta, coords = c("lon","lat"), crs = 4326, remove = FALSE)
  nb <- st_is_within_distance(pts, pts, dist = max_dist_m)

  i_idx <- rep.int(seq_along(nb), lengths(nb))
  j_idx <- unlist(nb, use.names = FALSE)
  keep <- i_idx < j_idx
  i_idx <- i_idx[keep]; j_idx <- j_idx[keep]

  elev <- meta$elevation
  both_elev <- !is.na(elev[i_idx]) & !is.na(elev[j_idx])
  ok_elev   <- !both_elev | (abs(elev[i_idx] - elev[j_idx]) <= max_elev_diff)

  edges <- cbind(i_idx[ok_elev], j_idx[ok_elev, drop = FALSE])
  g <- igraph::graph.empty(n = nrow(meta), directed = FALSE)
  if (nrow(edges)) g <- igraph::add_edges(g, as.vector(t(edges)))
  comp <- igraph::components(g)$membership
  meta$unified_id <- sprintf("STN_%04d", as.integer(factor(comp)))

  if (is.null(prefer_source_order)) prefer_source_order <- character(0)
  src_rank <- function(src) {
    if (!length(prefer_source_order)) return(rep(.Machine$integer.max, length(src)))
    r <- match(src, prefer_source_order, nomatch = NA_integer_)
    ifelse(is.na(r), .Machine$integer.max, r)
  }

  meta_canon <- meta %>%
    group_by(unified_id) %>%
    arrange(src_rank(source), is.na(elevation), row_id, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      member_keys = {
        tmp <- meta %>% group_by(unified_id) %>%
          summarise(keys = paste(unique(paste0(source,"::",station_id)), collapse="; "),
                    .groups="drop")
        tmp$keys[match(unified_id, tmp$unified_id)]
      }
    ) %>%
    dplyr::select(
      unified_id,
      station_id_canonical = station_id,
      lon, lat, elevation, source_canonical = source,
      member_keys
    )

  # ---------- 2) CROSSWALK (unique by construction now) ----------
  crosswalk <- meta %>%
    dplyr::select(source, station_id, id_key, unified_id, lon, lat, elevation)

  # ---------- 3) JOIN unified_id TO TIME SERIES ----------
  data2 <- data %>%
    mutate(station_id = as.character(station_id),
           source = as.character(source)) %>%
    left_join(crosswalk %>% dplyr::select(source, station_id, unified_id),
              by = c("source","station_id"),
              relationship = "many-to-one")  # assert no dup in crosswalk

  list(
    stations              = meta_canon,
    crosswalk             = crosswalk,
    data_with_unified     = data2,
    within_source_report  = within_source_report   # useful audit table (possibly empty)
  )
}



collapse_station_observations <- function(
    data_with_unified,
    prefer_source_order = c("Castellanos","Nicholson","GHCN","CRU","Meteostat","TAHMO","GISS","GSOD"),
    combine = c("first_non_na","mean"),
    variables = c("pre","tas","tasmin","tasmax","dewpoint")){
  combine <- match.arg(combine)

  # --- setup
  stopifnot(all(c("unified_id","year","month","source") %in% names(data_with_unified)))
  vars <- intersect(variables, names(data_with_unified))
  if (length(vars) == 0) stop("None of 'variables' found in data_with_unified.")

  # data.table
  requireNamespace("data.table")
  dt <- data.table::as.data.table(data_with_unified)
  dt <- dt[!is.na(unified_id)]  # ignore rows we can’t map

  # enforce simple types for speed
  dt[, `:=`(
    unified_id = as.character(unified_id),
    year = as.integer(year),
    month = as.integer(month),
    source = as.character(source)
  )]

  # keep only what we need
  keep_cols <- c("unified_id","year","month","source", vars)
  dt <- dt[, ..keep_cols]

  # result skeleton
  out <- unique(dt[, .(unified_id, year, month)])
  data.table::setkey(out, unified_id, year, month)

  # function to wide-cast per variable, then collapse across sources
  collapse_one_var <- function(v) {
    tmp <- dt[, .(unified_id, year, month, source, value = .SD[[1]]), .SDcols = v]

    # wide by source; if duplicates within (uid, y, m, source) exist, take first (not mean!)
    wide <- data.table::dcast(
      tmp,
      unified_id + year + month ~ source,
      value.var = "value",
      fun.aggregate = function(x) x[1]  # first non-NA per source/timestamp
    )

    # order columns by preference
    src_cols <- intersect(prefer_source_order, names(wide))
    # also keep any “other” sources (rare), appended after preferred
    other_cols <- setdiff(setdiff(names(wide), c("unified_id","year","month")), src_cols)
    src_cols_all <- c(src_cols, sort(other_cols))

    if (combine == "first_non_na") {
      # pick first non-NA across ordered sources
      # data.table’s fcoalesce works like dplyr::coalesce but for many columns
      value_vec <- do.call(data.table::fcoalesce, c(wide[, ..src_cols_all], list(NA_real_)))

      # provenance: which source supplied the value?
      prov <- rep(NA_character_, nrow(wide))
      if (length(src_cols_all)) {
        m <- as.matrix(!is.na(wide[, ..src_cols_all]))
        # find first TRUE in each row
        first_idx <- max.col(m, ties.method = "first")
        prov[ rowSums(m) > 0 ] <- src_cols_all[first_idx[rowSums(m) > 0]]
      }

      ans <- wide[, .(unified_id, year, month)]
      ans[, (v) := value_vec]
      ans[, paste0(v, "_source") := prov]
      ans
    } else { # combine == "mean"
      if (length(src_cols_all) == 0) {
        ans <- wide[, .(unified_id, year, month)]
        ans[, (v) := NA_real_]
        ans[, paste0(v, "_source") := NA_character_]
        return(ans)
      }
      mat <- as.matrix(wide[, ..src_cols_all])
      mean_vec <- rowMeans(mat, na.rm = TRUE)
      mean_vec[is.nan(mean_vec)] <- NA_real_

      # provenance: list sources that contributed (optional; comment out if not needed)
      prov <- apply(mat, 1L, function(x) {
        s <- src_cols_all[!is.na(x)]
        if (length(s)) paste(s, collapse = ";") else NA_character_
      })

      ans <- wide[, .(unified_id, year, month)]
      ans[, (v) := mean_vec]
      ans[, paste0(v, "_source") := prov]
      ans
    }
  }

  # collapse each variable and merge into output
  for (v in vars) {
    ans_v <- collapse_one_var(v)
    data.table::setkey(ans_v, unified_id, year, month)
    out <- ans_v[out]  # fast keyed join (keeps all rows in 'out')
  }

  # return as data.frame
  data.table::setcolorder(out, c("unified_id","year","month",
                                 as.vector(rbind(vars, paste0(vars, "_source")))))
  as.data.frame(out)
}


out <- dedupe_weather_stations(
  meta = all.metadata,
  data = all.data,
  max_dist_m = 2000,
  max_elev_diff = 100,
  prefer_source_order = c("Castellanos","Nicholson","GHCN","CRU","Meteostat","TAHMO","GISS","GSOD"),
  within_source_tol_m = 2000,
  within_source_elev_diff = 100,
  on_within_source_conflict = "warn_pick"  # set to "warn_pick" if you prefer auto-resolve
)



# Canonical stations (one per physical site)
out$stations

# Crosswalk from (source, station_id) to unified_id
out$crosswalk

# Your original data, unchanged but with unified_id attached
out$data_with_unified


dedup_monthly_fast <- collapse_station_observations(
  out$data_with_unified,
  prefer_source_order = c("Castellanos","Nicholson","GHCN","CRU","Meteostat","TAHMO","GISS","GSOD"),
  combine = "mean"  # or "mean"
)


metadata2keep <- out$stations
all.data2keep <- dedup_monthly_fast


saveRDS(metadata2keep,
        "./outputs/all.metadata.RDS")
saveRDS(all.data2keep,
        "./outputs/all.data.RDS")
keys <- c("unified_id","year","month")
sources <- all.data2keep %>%
  dplyr::select(all_of(keys), ends_with("_source")) %>%
  rename_with(~ sub("_source$", "", .x), ends_with("_source")) %>%
  pivot_longer(
    cols = -all_of(keys),
    names_to = "var",
    values_to = "source",
    values_drop_na = TRUE
  ) %>%
  mutate(source = case_when(grepl(";",source) ~ "Multiple",
                            TRUE ~ source))

source.sum <- sources %>%
  group_by(year,source,var) %>%
  summarise(Nstation = length(unique(unified_id)),
            .groups = "keep")

ggplot(data = source.sum) +
  geom_line(aes(x = year,
                y = Nstation,
                color = source)) +
  theme_bw() +
  facet_wrap(~var)

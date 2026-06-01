rm(list = ls())

all.data <- read.csv("/home/femeunier/Documents/projects/Drying.CB/data/station_data/Castellanos/Monthly all station records and predictions (1901-2010).csv")
meta.data <- read.csv("/home/femeunier/Documents/projects/Drying.CB/data/station_data/Castellanos/Station Statistics.csv")

selected.stations <-
  meta.data %>%
  filter(abs(LAT) <= 25,
         LONG >= -25,LONG <= 65)

all.data.selected <- all.data %>%
  dplyr::filter(ID %in% (selected.stations %>% pull(ID)))


long_df <- all.data.selected %>%
  pivot_longer(
    cols = matches("^PRE|^TYPE"),                  # All PREXX and TYPEXX columns
    names_to = c(".value", "MONTH"),               # .value = uses PRE/TYPE as column names
    names_pattern = "(PRE|TYPE)(\\d{2})"           # Extract base name and 2-digit month
  ) %>%
  mutate(MONTH = as.integer(MONTH)) %>%
  dplyr::select(-starts_with("TYPE")) %>%
  ungroup() %>%
  mutate(
    source    = sub("_.*$", "", ID),          # everything before first "_"
    station_id = sub("^.*?_", "", ID)          # everything after first "_"
  ) %>%
  dplyr::select(-ID) %>%
  rename(year = YEAR,
         month = MONTH,
         pre = PRE) %>%
  na.omit()

summary(long_df$pre)


df2save <- selected.stations %>%
  mutate(
    elevation = case_when(ELEV_REP==9999 ~ NA,
                          TRUE ~ ELEV_REP),
    source    = sub("_.*$", "", ID),          # everything before first "_"
    station_id = sub("^.*?_", "", ID)          # everything after first "_"
  ) %>%
  rename(lat = LAT,
         lon = LONG) %>%
  dplyr::select(source,station_id,lat,lon,elevation)
table(df2save$source)
hist(df2save$elevation)

saveRDS(df2save,
        "./data/station_data/Castellanos/metadata.Castellanos.RDS")

saveRDS(long_df,
        "./data/station_data/Castellanos/data.Castellanos.RDS")

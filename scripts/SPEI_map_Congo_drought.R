library(sf)
library(terra)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(geodata)
library(dplyr)
library(ggpattern)

# -----------------------------
# 1. Study region
# -----------------------------

xlim <- c(-5, 55)
ylim <- c(-15, 10)

world <- ne_countries(scale = "medium", returnclass = "sf")

countries <- world %>%
  filter(admin %in% c(
    "Democratic Republic of the Congo",
    "Republic of Congo",
    "Gabon",
    "Cameroon",
    "Central African Republic",
    "Angola",
    "Nigeria",
    "South Sudan",
    "Uganda",
    "Rwanda",
    "Burundi",
    "Tanzania",
    "Zambia"
  ))

label_pts <- countries %>%
  mutate(geometry = st_point_on_surface(geometry)) %>%
  mutate(
    label = case_when(
      admin == "Democratic Republic of the Congo" ~ "Democratic Republic\nof the Congo",
      admin == "Republic of the Congo" ~ "Republic of\nthe Congo",
      admin == "Central African Republic" ~ "Central African\nRepublic",
      TRUE ~ admin
    ),
    geometry = st_point_on_surface(geometry)
  )

africa <- world %>% filter(continent == "Africa")

# -----------------------------
# 2. Major lakes
# -----------------------------

lakes <- ne_download(
  scale = 10,
  type = "lakes",
  category = "physical",
  returnclass = "sf"
)

big_lakes <- lakes %>%
  st_make_valid() %>%
  st_crop(xmin = xlim[1], xmax = xlim[2],
          ymin = ylim[1], ymax = ylim[2])

# -----------------------------
# 3. Major rivers
# -----------------------------

rivers <- ne_download(
  scale = 10,
  type = "rivers_lake_centerlines",
  category = "physical",
  returnclass = "sf"
)

major_rivers <- rivers %>%
  st_make_valid() %>%
  st_crop(xmin = xlim[1], xmax = xlim[2],
          ymin = ylim[1], ymax = ylim[2])

# Optional: keep only larger rivers if scalerank exists
if ("scalerank" %in% names(major_rivers)) {
  major_rivers <- major_rivers %>% filter(scalerank <= 6)
}

# -----------------------------
# 4. Peatland extent
# -----------------------------
# Example: your peatland polygon file
# Replace by your own path, e.g. Xu et al. peatland extent shapefile

peat <- st_read("./data/shapefiles/Peatland.shp", quiet = TRUE) %>%
  st_transform(4326) %>%
  st_make_valid() %>%
  st_crop(xmin = xlim[1], xmax = xlim[2],
          ymin = ylim[1], ymax = ylim[2])

# -----------------------------
# 5. Prepare drought dataframe
# -----------------------------

E <- ext(-180,180,-25,25)
Mask <- read_sf("./data/shapefiles/Rainforests.shp")

system2("scp",c("hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/all.climate/MEM_SPEI6_all.years.tif","./outputs/"))

Zscores <- crop(mask(rast("./outputs/MEM_SPEI6_all.years.tif"),Mask),
                E)
dates <- time(Zscores)

df <- as.data.frame(Zscores,xy = TRUE) %>%
  pivot_longer(cols = -c(x,y))
df$date <- dates[match(df$name, names(Zscores))]


df_duration <- df %>%
  rename(lon = x, lat =y) %>%
  mutate(
    time = as.Date(date),
    stress = (value <= -1)
  ) %>%
  arrange(lon, lat, time) %>%
  group_by(lon, lat) %>%
  mutate(
    spell_id = cumsum(stress != lag(stress, default = first(stress)))
  ) %>%
  group_by(lon, lat, spell_id) %>%
  mutate(
    spell_duration = ifelse(stress, n(), 0)
  ) %>%
  ungroup()


df_duration_2024 <- df %>%
  rename(lon = x, lat =y) %>%
  mutate(
    time = as.Date(date),
    stress = value <= -1
  ) %>%
  filter(format(time, "%Y") == "2024") %>%
  arrange(lon, lat, time) %>%
  group_by(lon, lat) %>%
  summarise(
    duration_2024_months = sum(stress, na.rm = TRUE),
    .groups = "drop"
  )

summary(df_duration_2024$duration_2024_months)


df_map <- as.data.frame(mean(Zscores[[which(year(time(Zscores)) == 2024 &
                     month(time(Zscores)) >= 6)]]),
                     xy = TRUE) %>%
  rename(lon = x, lat = y) %>%
  filter(
    lon >= xlim[1], lon <= xlim[2],
    lat >= ylim[1], lat <= ylim[2]
  ) %>%
  mutate(class = cut(
    mean,
    breaks = c(-Inf, -2, -1.5, -1, -0.5, Inf),
    labels = c(
      "≤ -2.0",
      "-2.0 to -1.5",
      "-1.5 to -1.0",
      "-1.0 to -0.5",
      "> -0.5"
    )
  ))

# -----------------------------
# 6. Plot
# -----------------------------


ggplot() +
  geom_raster(
    data = df_map,
    aes(x = lon, y = lat, fill = mean)
  ) +
  geom_sf_pattern(
    data = peat,
    aes(),
    fill = NA,
    color = "darkolivegreen",
    linewidth = 0.1,

    pattern = "stripe",
    pattern_fill = "darkolivegreen",
    pattern_colour = "darkolivegreen",

    pattern_angle = 45,
    pattern_density = 0.15,
    pattern_spacing = 0.01,
    pattern_size = 0.4
  ) +
  geom_sf(
    data = big_lakes,
    fill = "lightblue",
    color = "lightblue",
    linewidth = 0.2
  ) +
  geom_sf(
    data = major_rivers,
    color = "steelblue",
    linewidth = 0.25
  ) +
  geom_sf(
    data = africa,
    fill = NA,
    color = "grey25",
    linewidth = 0.5
  ) +
  geom_sf_text(
    data = label_pts,
    aes(label = label),
    size = 3.2,
    family = "sans",  fontface = "bold",
    lineheight = 0.9,
    color = "grey15"
  ) +
  scale_fill_gradient2(limits = c(-3,3),
                       low = "darkred",high = 'darkblue',
                       oob = scales::squish) +
  # scale_fill_manual(
  #   values = c(
  #     "≤ -2.0" = "#7f0000",
  #     "-2.0 to -1.5" = "#d7301f",
  #     "-1.5 to -1.0" = "#fc8d59",
  #     "-1.0 to -0.5" = "#fee8c8",
  #     "> -0.5" = "grey"
  #   ),
  #   drop = FALSE,
  #   name = "Drought severity"
  # ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  scale_y_continuous(limits = c(-7,8)) +
  scale_x_continuous(limits = c(5,35)) +
  theme_void(base_size = 11) +
  theme(
    legend.position = "left",
    plot.title = element_text(size = 15, face = "bold"),
    plot.subtitle = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
  )



ggplot() +
  geom_raster(
    data = df_map,
    aes(x = lon, y = lat, fill = class)
  ) +
  geom_sf_pattern(
    data = peat,
    aes(),
    fill = NA,
    color = "darkolivegreen",
    linewidth = 0.1,

    pattern = "stripe",
    pattern_fill = "darkolivegreen",
    pattern_colour = "darkolivegreen",

    pattern_angle = 45,
    pattern_density = 0.15,
    pattern_spacing = 0.01,
    pattern_size = 0.4
  ) +
  geom_sf(
    data = big_lakes,
    fill = "lightblue",
    color = "lightblue",
    linewidth = 0.2
  ) +
  geom_sf(
    data = major_rivers,
    color = "steelblue",
    linewidth = 0.25
  ) +
  geom_sf(
    data = africa,
    fill = NA,
    color = "grey25",
    linewidth = 0.5
  ) +
  geom_sf_text(
    data = label_pts,
    aes(label = label),
    size = 3.2,
    family = "sans",  fontface = "bold",
    lineheight = 0.9,
    color = "grey15"
  ) +
  scale_fill_manual(
    values = c(
      "≤ -2.0" = "#7f0000",
      "-2.0 to -1.5" = "#d7301f",
      "-1.5 to -1.0" = "#fc8d59",
      "-1.0 to -0.5" = "#fee8c8",
      "> -0.5" = "grey"
    ),
    drop = FALSE,
    name = "Drought severity"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  scale_y_continuous(limits = c(-7,8)) +
  scale_x_continuous(limits = c(5,35)) +
  theme_void(base_size = 11) +
  theme(
    legend.position = "left",
    plot.title = element_text(size = 15, face = "bold"),
    plot.subtitle = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
  )



df_donut <- df_map %>%
  mutate(
    drought_class = case_when(
      mean <= -2.0 ~ "Extreme",
      mean <= -1.5 ~ "Severe",
      mean <= -1.0 ~ "Moderate",
      TRUE ~ "No drought"
    )
  ) %>%
  count(drought_class) %>%
  mutate(
    pct = 100 * n / sum(n),
    ymax = cumsum(pct),
    ymin = lag(ymax, default = 0)
  )

affected <- round(sum(df_donut$pct[df_donut$drought_class != "No drought"]), 1)

ggplot(df_donut,
       aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = drought_class)) +
  geom_rect(color = "white", linewidth = 0.5) +
  coord_polar(theta = "y") +
  xlim(c(0, 4)) +

  scale_fill_manual(
    values = c(
      "Extreme"   = "#7f0000",
      "Severe"    = "#d7301f",
      "Moderate"  = "#fc8d59",
      "No drought" = "grey92"
    ),
    breaks = c("Extreme", "Severe", "Moderate", "No drought"),
    name = NULL
  ) +

  labs(title = "Drought extent") +

  theme_void(base_size = 11) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 12)
  )

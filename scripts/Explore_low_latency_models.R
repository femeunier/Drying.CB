rm(list = ls())

library(terra)
library(lubridate)
library(dplyr)
library(ggplot2)
library(sf)
library(raster)
library(tidyr)

Ext <- ext(-180,180,-25,25)

IF <- crop(rast("./data/shapefiles/Intact_mask_mod.tif"),
           Ext)


# Radar.data <- rast("./outputs/Radar_all.years.tif")
# # IF <- rasterize(vect(read_sf("~/Documents/projects/Drying.CB/data/shapefiles/Rainforests.shp")),
# #                 Radar.data, field = "LC")
# IF <- rasterize(vect(read_sf("~/Documents/projects/Drying.CB/data/shapefiles/Rainforests.shp")), Radar.data, field = "LC")
# plot(IF)

df.models <- data.frame()
vars <- c("gpp","nbp","ter")

all.rast <- all.rast.anom <-
  all.rast.norm.anom <- list()

for (cvar in vars){
  files <- list.files("~/Downloads/data/",recursive = TRUE,full.names = TRUE,
                      pattern = paste0(".*", cvar, ".*\\.nc$"))
  files <- files[!grepl("Prog_2025",files)]

  all.rast[[cvar]] <- list()

  for (ifile in seq(1,length(files))){
    cfile <- files[ifile]

    cmodel <- strsplit(basename(cfile),"_")[[1]][1]

    print(paste0(cvar,"-",cmodel))

    mult <- ifelse(cmodel == "LPJ",
                   12,86400*365)
    crast <- crop(rast(cfile)*mult,
                  Ext)

    crast.r <- resample(crast, IF, method = "bilinear")
    ccrast.r <- mask(crast.r,IF,maskvalues=c(0,NA))

    if (cmodel == "ORCHIDEE" & grepl("2025",cfile)){
      ctimes <- seq(from = as.Date("2025-01-01"),
                    by = "1 month",
                    length.out = nlyr(ccrast.r))
    } else if (cmodel == "ORCHIDEE"){
      ctimes <- rev(seq(from = as.Date("2024-12-01"),
                        by = "-1 month",
                        length.out = nlyr(ccrast.r)))
    } else {
      ctimes <- time(ccrast.r)
    }

    print("===========")
    print(ctimes[length(ctimes)])
    print(time(ccrast.r[[length(ctimes)]]))
    print("===========")

    mths <- month(ctimes)
    clim <- tapp(ccrast.r, index = mths, fun = mean, na.rm = TRUE)
    clim.sd <- tapp(ccrast.r, index = mths, fun = sd, na.rm = TRUE)

    ccrast.r.anom <- ccrast.r - clim[[mths]]
    ccrast.r.norm.anom <- ccrast.r.anom/clim.sd[[mths]]

    terra::time(ccrast.r) <- terra::time(ccrast.r.anom) <-
      terra::time(ccrast.r.norm.anom) <- ctimes

    df.models <- bind_rows(df.models,
                           data.frame(model = cmodel,
                                      time = ctimes,
                                      var = cvar,
                                      value = global(ccrast.r,mean,na.rm = TRUE)[["mean"]],
                                      anom = global(ccrast.r.anom,mean,na.rm = TRUE)[["mean"]],
                                      anom.norm = global(ccrast.r.norm.anom,mean,na.rm = TRUE)[["mean"]]))

    all.rast[[cvar]][[cmodel]] <-
      c(all.rast[[cvar]][[cmodel]],
        ccrast.r)

    all.rast.norm.anom[[cvar]][[cmodel]] <-
      c(all.rast.norm.anom[[cvar]][[cmodel]],
        ccrast.r.norm.anom)

    all.rast.anom[[cvar]][[cmodel]] <-
      c(all.rast.anom[[cvar]][[cmodel]],
        ccrast.r.anom)

  }
}


# We add CTEHR
files <- list.files("~/Downloads/data/CTE-HR-Global_202505-05x05/",full.names = TRUE)

rast.gpp <- rast.nep <- rast.ter <-
  list()
for (ifile in seq(1,length(files))){
  cfile <- files[ifile]
  cdate <- strsplit(tools::file_path_sans_ext(basename(cfile)),"_")[[1]][5]
  ctime <- paste0(substr(cdate,1,4),"/",substr(cdate,5,6),"/01")

  print(ctime)
  crast <- crop(rast(cfile),Ext)

  crast.r <- resample(crast, IF, method = "bilinear")
  ccrast.r <- mask(crast.r,IF,maskvalues=c(0,NA))

  Ndays <- (nlyr(crast)-1)/2
  pos.gpp <- 2:(Ndays+1)
  pos.ter <- (Ndays+2):nlyr(crast)
  gpp <- mean(-ccrast.r[[pos.gpp]])*365 *86400 *0.012
  ter <- mean(ccrast.r[[pos.ter]])*365 *86400 *0.012
  nep <- gpp - ter

  if (ifile == 1){
    pos.mask <- (gpp == 0) & (ter == 0)
  }

  gpp[pos.mask] <- ter[pos.mask] <-
    nep[pos.mask] <- NA

  final.rast.gpp <- gpp
  final.rast.nep <- nep
  final.rast.ter <- ter

  terra::time(final.rast.gpp) <-
    terra::time(final.rast.nep) <-
    terra::time(final.rast.ter) <-
    as.Date(ctime)

  rast.gpp[[ifile]] <- final.rast.gpp
  rast.ter[[ifile]] <- final.rast.ter
  rast.nep[[ifile]] <- final.rast.nep

}

rast.gpp <- rast(rast.gpp)
rast.ter <- rast(rast.ter)
rast.nep <- rast(rast.nep)

ctimes <- time(rast.nep)

mths <- month(ctimes)
clim.gpp <- tapp(rast.gpp, index = mths, fun = mean, na.rm = TRUE)
clim.gpp.sd <- tapp(rast.gpp, index = mths, fun = sd, na.rm = TRUE)

clim.ter <- tapp(rast.ter, index = mths, fun = mean, na.rm = TRUE)
clim.ter.sd <- tapp(rast.ter, index = mths, fun = sd, na.rm = TRUE)

clim.nep <- tapp(rast.nep, index = mths, fun = mean, na.rm = TRUE)
clim.nep.sd <- tapp(rast.nep, index = mths, fun = sd, na.rm = TRUE)

anom.gpp <- rast.gpp - clim.gpp[[mths]]
anom.ter <- rast.ter - clim.ter[[mths]]
anom.nep <- rast.nep - clim.nep[[mths]]

anom.norm.gpp <- anom.gpp/clim.gpp.sd[[mths]]
anom.norm.nep <- anom.nep/clim.nep.sd[[mths]]
anom.norm.ter <- anom.ter/clim.ter.sd[[mths]]

all.rast[["Sib4"]][["gpp"]] <- rast.gpp
all.rast[["Sib4"]][["nbp"]] <- rast.nep
all.rast[["Sib4"]][["ter"]] <- rast.ter

all.rast.anom[["Sib4"]][["gpp"]] <- anom.gpp
all.rast.anom[["Sib4"]][["nbp"]] <- anom.nep
all.rast.anom[["Sib4"]][["ter"]] <- anom.ter

all.rast.norm.anom[["Sib4"]][["gpp"]] <- anom.norm.gpp
all.rast.norm.anom[["Sib4"]][["nbp"]] <- anom.norm.nep
all.rast.norm.anom[["Sib4"]][["ter"]] <- anom.norm.ter

df.models <- bind_rows(df.models,
                       data.frame(model = "Sib4",
                                  time = ctimes,
                                  var = "gpp",
                                  value = global(rast.gpp,mean,na.rm = TRUE)[["mean"]],
                                  anom = global(anom.gpp,mean,na.rm = TRUE)[["mean"]],
                                  anom.norm = global(anom.norm.gpp,mean,na.rm = TRUE)[["mean"]]),
                       data.frame(model = "Sib4",
                                  time = ctimes,
                                  var = "nbp",
                                  value = global(rast.nep,mean,na.rm = TRUE)[["mean"]],
                                  anom = global(anom.nep,mean,na.rm = TRUE)[["mean"]],
                                  anom.norm = global(anom.norm.nep,mean,na.rm = TRUE)[["mean"]]),
                       data.frame(model = "Sib4",
                                  time = ctimes,
                                  var = "ter",
                                  value = global(rast.ter,mean,na.rm = TRUE)[["mean"]],
                                  anom = global(anom.ter,mean,na.rm = TRUE)[["mean"]],
                                  anom.norm = global(anom.norm.ter,mean,na.rm = TRUE)[["mean"]]))


##################################################
# Merge all

all.rast.together <-
  all.rast.norm.anom.together <-
  all.rast.anom.together <- list()

compt <- 1
for (cvar in vars){

  models <- names(all.rast[[cvar]])

  for(cmodel in models){

    all.rast[[cvar]][[cmodel]] <- rast(all.rast[[cvar]][[cmodel]])
    all.rast.anom[[cvar]][[cmodel]] <- rast(all.rast.anom[[cvar]][[cmodel]])
    all.rast.norm.anom[[cvar]][[cmodel]] <- rast(all.rast.norm.anom[[cvar]][[cmodel]])


    crast <- all.rast[[cvar]][[cmodel]]
    crast2 <- all.rast.norm.anom[[cvar]][[cmodel]]
    crast3 <- all.rast.anom[[cvar]][[cmodel]]

    names(crast) <- paste0(cmodel,"_",cvar,"_",1:nlyr(crast))
    names(crast2) <- paste0(cmodel,"_",cvar,"_norm.anom_",1:nlyr(crast))
    names(crast3) <- paste0(cmodel,"_",cvar,"_anom_",1:nlyr(crast))


    all.rast.together[[compt]] <- crast
    all.rast.anom.together[[compt]] <- crast3
    all.rast.norm.anom.together[[compt]] <- crast2

    compt <- compt + 1
  }
}

all.rast.together <- rast(all.rast.together)
all.rast.anom.together <- rast(all.rast.anom.together)
all.rast.norm.anom.together <- rast(all.rast.norm.anom.together)

pos <- (grepl("gpp",names(all.rast.norm.anom.together)) &
               year(time(all.rast.norm.anom.together)) == 2024 &
               month(time(all.rast.norm.anom.together)) >= 6)


r.sel <- all.rast.norm.anom.together[[pos]]

model.names <- as.factor(sapply(strsplit(names(r.sel), "_"), "[[", 1))
mdls <- as.numeric(model.names)


r_value <- all.rast.together[[pos]]
r_anom  <- all.rast.anom.together[[pos]]

r_clim <- r_value - r_anom

model.names <- as.factor(sapply(strsplit(names(r_value), "_"), "[[", 1))
mdls <- as.numeric(model.names)

reduc.per.model <- tapp(
  r_anom / r_clim,
  index = mdls,
  fun = mean,
  na.rm = TRUE
)

names(reduc.per.model) <- levels(model.names)

summary(mean(reduc.per.model))

anom.per.model <- tapp(r_anom, index = mdls, fun = mean, na.rm = TRUE)
clim.per.model <- tapp(r_clim, index = mdls, fun = mean, na.rm = TRUE)

reduc.per.model <- anom.per.model / clim.per.model
names(reduc.per.model) <- levels(model.names)

plot(reduc.per.model)

anomalies.per.model <- tapp(r.sel, index = mdls, fun = mean, na.rm = TRUE)
names(anomalies.per.model) <- levels(model.names)

df.anomalies.per.model <- as.data.frame(anomalies.per.model, xy = TRUE) %>%
  rename(lon = x, lat = y) %>%
  pivot_longer(cols = -c(lon, lat),
               names_to = "model")

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  geom_tile(data = df.anomalies.per.model,
              aes(x = lon, y = lat, fill = value)) +
  geom_sf(data = world, fill = NA, color = "black", linewidth = 0.3) +
  scale_fill_gradient2(low = "darkred",mid = "white",high = "darkblue",
                       midpoint = 0, oob = scales::squish) +
  theme_bw() +
  facet_wrap(~model) +
  coord_sf(xlim = c(5, 35), ylim = c(-10, 5), expand = FALSE) +
  labs(
    fill = "Z-anomalies"
  )


ggplot(data = df.models %>%
         filter(year(time) >= 1990)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(aes(x = time, y = anom.norm,
                color = model)) +
  facet_wrap(~var,scales = "free") +
  theme_bw()

df.models %>%
  filter(year(time) >= 1990, year(time)<2025) %>%
  group_by(model,var) %>%
  filter(anom.norm == anom.norm[which.min(anom.norm)])

df.when <- df.models %>%
  filter(year(time) >= 1990, year(time)<2025) %>%
  group_by(model,var) %>%
  arrange(anom.norm) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  arrange(model,var)

ggplot(data = df.models %>%
         filter(year(time) >= 1990)) +
  geom_rect(inherit.aes = FALSE,
            aes(xmin = as.Date("2024-01-01"),xmax = as.Date("2025-01-01"),
                ymin = -Inf,ymax = Inf),
            fill = "lightgrey",
            alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(aes(x = time, y = anom.norm,
                color = model),
            size = 0.2) +
  facet_wrap(~var,scales = "free") +
  geom_point(data = df.when,inherit.aes = FALSE,
             aes(x = time, y = anom.norm, color = model)) +
  theme_bw()

Nsmooth=6


df.models2 <- df.models %>%
  group_by(time) %>%
  summarise(anom.norm = mean(anom.norm,na.rm = TRUE),
            .groups = "keep") %>%
  ungroup() %>%
  mutate(anom.norm_ma = slide_dbl(anom.norm, mean, .before = (Nsmooth-1), .complete = TRUE, na.rm = TRUE))

ggplot(data = df.models %>%
         filter(year(time) >= 1990, year(time) < 2025)) +
  geom_rect(inherit.aes = FALSE,
            aes(xmin = as.Date("2024-01-01"),xmax = as.Date("2025-01-01"),
                ymin = -Inf,ymax = Inf),
            fill = "lightgrey",
            alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(aes(x = time, y = anom.norm,
                group = model), color = "grey",
            size = 0.2) +
  geom_line(data = df.models2 %>%
              filter(year(time) >= 1990, year(time) < 2025),
            aes(x = time, y = anom.norm_ma), color = "black",
            size = 1) +
  facet_wrap(~var,scales = "free") +
  theme_bw()



df.models3 <- df.models %>%
  mutate(year = year(time),
         month = month(time)) %>%
  group_by(var,year,month) %>%
  summarise(anom.norm = mean(anom.norm,na.rm = TRUE),
            .groups = "keep") %>%
  ungroup() %>%
  mutate(anom.norm_ma = slide_dbl(anom.norm, mean, .before = (Nsmooth-1), .complete = TRUE, na.rm = TRUE))

ggplot(data = df.models %>%
         filter(var == "gpp")) +
  geom_line(aes(x = year(time) + (month(time) - 1/2)/12, y = anom.norm,
                group = model),color = "grey") +
  scale_x_continuous(limits = c(2020,2025)) +
  geom_line(data = df.models3 %>%
              filter(var == "gpp") %>%
              filter(year >= 1990, year < 2025),
            aes(x = year + (month - 1/2)/12, y = anom.norm_ma), color = "black",
            size = 1) +
  geom_hline(yintercept = 0,linetype = 2) +
  labs(x = "", y = "GPP anomaly") +
  theme_bw() +
  theme(text = element_text(size = 24))

ggplot(data = df.models) +
  geom_line(aes(x = year(time) + (month(time) - 1/2)/12, y = anom,
                color = model)) +
  scale_x_continuous(limits = c(2020,2025)) +
  geom_hline(yintercept = 0,linetype = 2) +
  facet_wrap(~var) +
  labs(x = "") +
  theme_bw()


ggplot(data = df.models) +
  geom_line(aes(x = year(time) + (month(time) - 1/2)/12, y = value,
                color = model)) +
  scale_x_continuous(limits = c(2020,2025)) +
  geom_hline(yintercept = 0,linetype = 2) +
  facet_wrap(~var) +
  labs(x = "") +
  theme_bw()

A <- df.models %>%
  filter(var == "nbp") %>%
  mutate(year = year(time)) %>%
  group_by(model,year) %>%
  summarise(m = mean(value,na.rm = TRUE))

B <- A %>%
  group_by(year) %>%
  summarise(MEM = mean(m))


ggplot(data = B %>%
         filter(year <= 2024),
       aes(x = year, y = MEM)) +
  geom_hline(yintercept = 0,linetype = 2) +
  geom_line() +
  theme_bw()

S <- (as.data.frame(IF) %>%
        filter(lccs_class==1) %>%
        nrow())*(0.25*0.25*111*111)*1e6


B %>%
  filter(year == 2024) %>% pull(MEM)*S/1e12

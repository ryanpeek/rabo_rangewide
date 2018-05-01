# map of samples
library(sf)
library(here)
library(dplyr)
library(mapview)


site <- "rub_rasi"
siterds <- readRDS(paste0(here(), "/data_output/bamlist_dat_", site, "_25k_thresh.rds"))

summary(siterds)

# filter NAs
df <- siterds %>% filter(!is.na(lat))

# make sf
df_locs<- st_as_sf(df, 
                   coords=c("lon", "lat"),
                   remove=F,
                   crs=4326)
mapview(df_locs)

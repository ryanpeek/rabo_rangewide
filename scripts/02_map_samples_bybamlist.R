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




# MAP HUC8 ----------------------------------------------------------------


h8 <- st_read(unzip("data/HUC8_named_westcoast.zip"))

# then remove raw files since file is added in memory
file.remove(list.files(pattern = "HUC8_named_westcoast*",recursive = F))

mapview(h8)

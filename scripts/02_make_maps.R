# map of samples

# Load Libraries ----------------------------------------------------------

library(sf)
library(here)
library(dplyr)
library(readr)
library(mapview)
options(scipen = 12)
options(tibble.print_min = Inf) # print all rows 
# Load metadat ------------------------------------------------------------

siterds <- read_rds("data_output/rapture_metadata_rabo_quant.rds")
summary(siterds)

# filter NAs
#df <- siterds %>% filter(!is.na(lat))

# MAKE SF and MAPVIEW -----------------------------------------------------

df_locs<- st_as_sf(siterds, 
                   coords=c("lon", "lat"),
                   remove=F,
                   crs=4326)
mapview(df_locs)

# Find Areas with Few Samples ---------------------------------------------

# total per HUC8
df_locs %>% group_by(HUC_8) %>% tally 

# HUC8's with fewer than 8 samples:
df_locs %>% group_by(HUC_8) %>% mutate(tot_n=n()) %>% 
  filter(tot_n < 3) %>% mapview()

# MAP HUC8 ----------------------------------------------------------------

h8 <- st_read(unzip("data/HUC8_named_westcoast.zip"))

# then remove raw files since file is added in memory
file.remove(list.files(pattern = "HUC8_named_westcoast*",recursive = F))

# add h6 field
h8$HUC_6 <- strtrim(as.character(h8$HUC_8), width = 6)

#mapview(h8) + mapview(df_locs)

head(h8)

# MAKE H6 -----------------------------------------------------------------

h6 <- h8 %>%  
  group_by(HUC_6) %>% summarize()
mapview(h6)

# JOIN W METADATA ---------------------------------------------------------

# HUC6's with fewer than 8 samples:
(df_h6_tally <- df_locs %>% group_by(HUC_6) %>% tally %>% 
   as.data.frame() %>% # remove geometry
   select(HUC_6, n) %>% 
   mutate(HUC_6=as.character(HUC_6)))

# now join with H6 spatial layer
h6 <- left_join(h6, df_h6_tally, by="HUC_6")

summary(h6)

mapview(h6, zcol="n") + mapview(df_locs)

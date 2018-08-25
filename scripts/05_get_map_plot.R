
# GET MAP AND PLOT ------------------------------------------------------

library(tidyverse)
library(mapview)
library(sf)

# GET BAMLIST AND METADATA ------------------------------------------------

bamfile <- "all_rabo_filt10_100k"

# Get ID and pop info for each individual from bamlists
bams <- read.table(paste0("data_output/bamlists/",bamfile, "_thresh.bamlist"),stringsAsFactors = F, header = F)
bams$V2 <- sub('\\..*$', '', basename(bams$V1)) # remove the path and file extension

# get metadat
metadat <- read_rds(path = "data_output/rapture_metadata_rabo_quant.rds")

# join together
annot <- left_join(bams, metadat, by=c("V2"="Seq")) %>% select(-V1) # join with the metadata

# check for duplicates?
annot[duplicated(annot$V2),]
metadat[duplicated(metadat$Seq),]

# change ecoregion to factor
annot$EcoRegion <- as.factor(annot$EcoRegion)

# make factor
length(annot[is.na(annot$EcoRegion),])
annot[is.na(annot$EcoRegion),] 
summary(annot$EcoRegion)

# recode levels based on Shaffer paper (E=Southern Sierra Nevada, NE=Northern Sierra Nevada, NW=North Coast, W=Central Coast, SW=South Coast, NEW::::Sierra/Basin Range)

# view groups with:
# annot %>% distinct(EcoRegion)
# annot %>% filter(EcoRegion=="Sierra Nevada") %>% group_by(River) %>% tally()
# annot %>% filter(grepl("Southern OR Coastal|North Coast|Cascades|Klamath North Coast|Northern CA Coastal Foothills", EcoRegion)) %>% group_by(River) %>% distinct(River)
# annot %>% filter(grepl("Sierra/Basin Range", EcoRegion)) %>% group_by(River) %>% distinct(River)

annot<- annot %>% 
  mutate(admix_groups = case_when(
    grepl("STAN|TUO|SFA|CALAV", River) ~ "East", # southern siera
    grepl("ANTV|BEAR|DEER|MFA|MFY|NFA|NFMFA|NFY|SFY|RUB", River) ~ "North-East", # northern sierra
    grepl("CHETCO|SFEEL|VANDZ|TRIN|MAT|KLAM|SSANTIAM|PUT|MAD|LAGUN|SUMPQUA|RUSS|SMITH|EEL", River) ~ "North-West", # North Coast
    grepl("NFF|FEA", River) ~ "Feather-North", # feather
    grepl("PAJ|ALA|DRY|SOQUEL", River) ~ "West", # Central Coast
    grepl("SANCARP", River) ~ "South-West") # South Coast
  )


# 05c. MAKE QUICK MAP ----------------------------------------------------------


# jitter coords slightly for viewing
annot_sf <- annot %>% filter(!is.na(lat)) %>%  
  distinct(lat,.keep_all = T) %>% 
  select(River, Site, SampleID, LabID, SPP_ID, admix_groups, lat, lon, elev_m:Locality,Locality_details, EcoRegion)

# make sf:
annot_sf <- st_as_sf(annot_sf, 
                     coords = c("lon", "lat"), 
                     remove = F, # don't remove these lat/lon cols from df
                     crs = 4326) # add projection

# check proj
st_crs(annot_sf)

# write out and delete there's a file with same name (will give warning if first time saving)
st_write(annot_sf, paste0("data_output/sites_",bamfile, ".shp"), delete_dsn = T)

# make a quick mapview map
mapview(annot_sf, zcol="admix_groups") %>% addMouseCoordinates()

# Get Counties  -----------------------------------------------------------

library(USAboundaries)
library(purrr)

# Pick a State
state_names <- c("california", "oregon")

# Download STATE data and add projection
CA<-us_states(resolution = "low", states = state_names) #%>%
#st_transform(crs = 32610)
st_crs(CA)

# get COUNTY data for a given state
counties <- us_counties(resolution = "low", states=state_names) %>% # use list of state(s) here
  mutate(lon=map_dbl(geometry, ~st_centroid(.x)[[1]]), # add centroid values for labels
         lat=map_dbl(geometry, ~st_centroid(.x)[[2]])) # add centroid values for labels

# add color palette: 
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# plot
ggplot() + geom_sf(data=CA, fill="gray80") + geom_sf(data=counties, col="gray20", alpha=0.8) +
  geom_sf(data=annot_sf, aes(color=admix_groups), size=3.5) + 
  scale_color_manual("Admix Groups", 
                     values = c("East"=cbbPalette[1], 
                                "North-East"=cbbPalette[2], 
                                "North-West"=cbbPalette[3],
                                "Feather-North"=cbbPalette[4],
                                "West"=cbbPalette[5], 
                                "South-West"=cbbPalette[6])) +
  theme_bw(base_family = "Roboto Condensed") 
  
  #scale_color_viridis(discrete = TRUE) 

ggsave(filename = paste0("figs/maps_", bamfile, "_admix_groups.png"), width = 8, height = 11, 
              units = "in", dpi = 300)




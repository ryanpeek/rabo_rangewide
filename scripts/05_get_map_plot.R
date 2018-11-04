
# LIBRARIES

library(tidyverse)
library(mapview)
library(sf)
library(USAboundaries)
library(purrr)
library(ggsn)
library(dplyr)
library(smoothr)

# GET BAMLIST AND METADATA ------------------------------------------------

bamfile <- "all_rabo_filt10_1_100k"
# all_rabo_filt_100k
# all_rabo_filt10_1_100k
# all_rabo_filt_1_100k

# Get ID and pop info for each individual from bamlists
bams <- read.table(paste0("data_output/bamlists/",bamfile, "_thresh.bamlist"),stringsAsFactors = F, header = F)
bams$V2 <- sub('\\..*$', '', basename(bams$V1)) # remove the path and file extension

# get metadat and fix sandy bar space and deer ck
metadat <- read_rds(path = "data_output/rapture_metadata_rabo_quant.rds")

# fix erroneous Lassen Locality
# metadat <- metadat %>% mutate(Locality_details=if_else(SampleID=="RAP-327", "Mineral Bar, Iowa Hill Rd, at NF American River, Placer County, el 2400 feet, 5 July, 1994", Locality_details),
#                               Locality = if_else(SampleID=="RAP-327", "NFA", Locality),
#                               River = if_else(SampleID=="RAP-327","NFA", River),
#                               Site = if_else(SampleID =="RAP-327", as.character(NA), Site),
#                               lat = if_else(SampleID=="RAP-327", as.numeric(39.11115), lat),
#                               lon = if_else(SampleID=="RAP-327", as.numeric(-120.9168), lon),
#                               elev_m = if_else(SampleID=="RAP-327", 386.4710, elev_m),
#                               HUC_6 = if_else(SampleID=="RAP-327", as.integer(180201), HUC_6),
#                               HUC_8 = if_else(SampleID=="RAP-327", as.integer(18020128), HUC_8),
#                               HUC_10 = if_else(SampleID=="RAP-327", as.integer(1802012801), HUC_10),
#                               HU_6_NAME = if_else(SampleID=="RAP-327", "Lower Sacramento", HU_6_NAME),
#                               HU_8_NAME = if_else(SampleID=="RAP-327", "North Fork American", HU_8_NAME),
#                               HU_10_NAME = if_else(SampleID=="RAP-327", "Upper North Fork American River", HU_10_NAME),
#                               HU_12_NAME = if_else(SampleID=="RAP-327", "Indian Creek-North Fork American River", HU_12_NAME),
#                               county = if_else(SampleID == "RAP-327", "Placer", county),
#                               dist2stream_m = if_else(SampleID == "RAP-327", as.numeric(37.22), dist2stream_m),
#                               COMID = if_else(SampleID == "RAP-327", as.integer(14992581), COMID),
#                               GNIS_ID = if_else(SampleID == "RAP-327", as.integer(267066), GNIS_ID),
#                               GNIS_NAME = if_else(SampleID == "RAP-327", "North Fork American River", GNIS_NAME),
#                               NHDlength_km = if_else(SampleID == "RAP-327", as.numeric(3.238), NHDlength_km),
#                               NHD_StreamLevel = if_else(SampleID == "RAP-327", as.integer(2), NHD_StreamLevel),
#                               NHD_StreamOrder = if_else(SampleID == "RAP-327", as.integer(4), NHD_StreamOrder),
#                               NHD_Tot_DA_sqkm = if_else(SampleID == "RAP-327", as.numeric(605.3175), NHD_Tot_DA_sqkm),
#                               EcoRegion = if_else(SampleID == "RAP-327", "Sierra Nevada", EcoRegion))
# 

# fix trinity spaces
# metadat$Locality<-gsub(pattern = "[[:space:]]", replacement = "-", x = metadat$Locality)

# fix deer-clearck/ deer-clec
# metadat$Locality <- gsub(pattern="deer-clearck", replacement = "deer-clec", x=metadat$Locality)

# save metadat back out:
# write_rds(metadat, path = "data_output/rapture_metadata_rabo_quant.rds")

# join together
annot <- left_join(bams, metadat, by=c("V2"="Seq")) %>% select(-V1) # join with the metadata

# check for duplicates?
annot[duplicated(annot$V2),]
metadat[duplicated(metadat$Seq),]

# add admix groupings
annot<- annot %>% 
  mutate(
    admix_orig = case_when(
      grepl("STAN|TUO|CALAV", River) ~ "East", # southern siera
      grepl("SFA", River) ~ "Unknown",
      grepl("NFF|FEA|ANTV|BEAR|DEER|MFA|MFY|NFA|NFMFA|NFY|SFY|RUB", River) ~ "North-East", # northern sierra
      grepl("CHETCO|SFEEL|COW|VANDZ|TRIN|MAT|KLAM|SSANTIAM|PUT|^MAD$|LAGUN|SUMPQUA|RUSS|SMITH|EEL", River) ~ "North-West", # North Coast
      #grepl("NFF|FEA", River) ~ "North-Feather", # feather
      grepl("PAJ|ALA|DRY|SOQUEL", River) ~ "West", # Central Coast
      grepl("SANCARP|SALIN", River) ~ "South-West"), # South Coast
    admix_groups = case_when(
      grepl("STAN|TUO|CALAV|SFA", River) ~ "East", # southern siera
      grepl("BEAR|DEER|MFA|MFY|NFA|NFMFA|NFY|SFY|RUB", River) ~ "North-East", # northern sierra
      grepl("CHETCO|SFEEL|COW|VANDZ|TRIN|MAT|KLAM|SSANTIAM|PUT|^MAD$|LAGUN|SUMPQUA|RUSS|SMITH|EEL", River) ~ "North-West", # North Coast
      grepl("NFF|FEA", River) ~ "North-Feather", # feather
      grepl("PAJ|ALA|DRY|SOQUEL", River) ~ "West", # Central Coast
      grepl("SANCARP|SALIN", River) ~ "South-West") # South Coast
  )

ords_admix_grps <- c("East", "North-East", "North-Feather", "North-West", "South-West", "West")

annot$admix_groups <- factor(annot$admix_groups, levels = ords_admix_grps)

# MAKE SPATIAL ------------------------------------------------------------

# jitter coords slightly for viewing
annot_sf <- annot %>% filter(!is.na(lat)) %>%  
  group_by(Locality) %>% add_tally() %>% 
  distinct(Locality,.keep_all = T) %>% 
  select(River, Site, SampleID, LabID, SPP_ID, admix_groups, admix_orig, lat, lon, elev_m:Locality,Locality_details, EcoRegion, n) %>% 
  ungroup() %>% 
  arrange(admix_groups, Locality) %>% 
  mutate(siteID=row_number()) # for matching labels, etc

# save out table of sites
annot_out <- select(annot_sf, siteID, Locality, River, Site, admix_groups, admix_orig, lat, lon, HUC_6, county, n) %>% 
  dplyr::rename("clade" = admix_groups, "n_samples"=n)

#write_csv(annot_out, path = paste0("data_output/table_site_localities_clades_", bamfile, ".csv"))


# make sf:
annot_sf <- st_as_sf(annot_sf, 
                     coords = c("lon", "lat"), 
                     remove = F, # don't remove these lat/lon cols from df
                     crs = 4326) # add projection

# check proj
st_crs(annot_sf)

# write out and delete there's a file with same name (will give warning if first time saving)
#st_write(annot_sf, paste0("data_output/sites_",bamfile, ".shp"), delete_dsn = T)
#save(annot_sf, file = paste0("data_output/sites_",bamfile, ".rda"))

# MAKE QUICK MAP ----------------------------------------------------------

# make a quick mapview map
mapview(annot_sf, zcol="admix_groups") %>% addMouseCoordinates()

# SIMPLIFY FIX RANGE -----------------------------------------------------------

# get shp of range
rb_range <- st_read(unzip("data/shps/Rb_Potential_Range_CAandOR.zip")) %>% 
   st_transform(crs = 4326)
file.remove(list.files(pattern="Rb_Potential_Range_CAandOR", recursive = F))
st_crs(rb_range)
plot(rb_range$geometry)

# simplify
rb_range_simple <- st_simplify(rb_range, dTolerance = .05)
plot(rb_range_simple$geometry, col="skyblue")
 
# smooth
rb_smooth <- smooth(rb_range_simple, method = "ksmooth") # or spline / ksmooth
plot(rb_smooth$geometry, col="maroon")
 
# rm small bits
area_thresh <- units::set_units(400, km^2)
rb_dropped <- drop_crumbs(rb_smooth, threshold = area_thresh)
plot(rb_dropped$geometry, col="orange")
 
# fill holes
area_thresh <- units::set_units(800, km^2)
rb_filled <- fill_holes(rb_dropped, threshold = area_thresh)
plot(rb_filled$geometry, col="purple")
 
# now add fields
rb_filled <- rb_filled %>% 
  mutate(state = case_when(
    SEASON=="Y" ~ "CA",
    is.na(SEASON) ~ "OR"),
    keeppoly = "Y") %>%
  select(-SEASON, -SHAPE_NAME, -Id)
 
# get OR
or_rb <- rb_filled %>% filter(state=="OR")
or_rb_buff <- st_buffer(or_rb, dist = .08)
 
plot(rb_filled$geometry, border="gray", lwd=1)
plot(or_rb_buff$geometry, border="blue", add=T)
plot(or_rb$geometry, border="red", add=T)
 
# get CA only
ca_rb <- rb_filled %>% filter(state=="CA")
 
# combine
rb_all <- st_union(or_rb_buff, ca_rb, by_feature = "state")
plot(rb_all$geometry, col="pink")
 
rb_diss <- rb_all %>%
  group_by(keeppoly) %>%
  summarize()

plot(rb_diss$geometry, col="purple")

# write out
st_write(rb_diss, "data_output/rabo_range_simple.shp", delete_dsn = T)

# rm interm files:
rm(ca_rb, or_rb, or_rb_buff, rb_all, rb_dropped, rb_filled, rb_range, rb_range_simple, rb_smooth)


# 01. GET SHAPES --------------------------------------------------------------

#bamfile <- "all_rabo_filt_100k"
bamfile <- "all_rabo_filt10_1_100k"

# get shp file with data
#annot_sf <- st_read(paste0("data_output/sites_", bamfile, ".shp"))
#st_crs(annot_sf)

# revised range:
rb_range <- st_read("data/rabo_range_simple.shp")
# revised range modified in GIS to include LASSEN
# rb_range <- st_read("data/shps/rabo_range_simple_v2.shp") # 


st_crs(rb_range)
plot(rb_range$geometry, border="orange", lwd=2) 

# load annot_sf data:
load(paste0("data_output/sites_",bamfile, ".rda"))

# lakes and streams
lakes <- st_read("data/shps/NHDWaterbody_lakes_greater_15sqkm.shp") 
lakes <- st_transform(lakes, crs = 4326)
st_crs(lakes)

rivers <- st_read("data/shps/CA_major_river_SJV_disconnected.shp")
rivers <- st_transform(rivers, crs = 4326)
st_crs(rivers)

# 03. GET SHAFFER SITES -------------------------------------------------------

# get shaffer sites:
shaff <- read_csv("data/mccartney_shaffer_samples.csv")

# make spatial
shaff_sf <- st_as_sf(shaff, 
                     coords = c("lon", "lat"), 
                     remove = F, # don't remove these lat/lon cols from df
                     crs = 4326) # add projection

# check proj
st_crs(shaff_sf)

# quick map
#mapview(shaff_sf) + mapview(annot_sf)

# GET STATE AND COUNTIES --------------------------------------------------

# Get states
state_names <- c("california", "oregon") # for RABO range
bg_states <- c("california", "oregon", "washington", "idaho","montana", "nevada","utah", "arizona")
# Download STATE data and add projection
CA<-us_states(resolution = "low", states = state_names)
bgSTs <- us_states(resolution = "low", states=bg_states)
#st_crs(CA)
#plot(st_geometry(bgSTs))

# get COUNTY data for a given state
counties <- us_counties(resolution = "low", states=state_names) %>% # use list of state(s) here
  mutate(lon=map_dbl(geometry, ~st_centroid(.x)[[1]]), # add centroid values for labels
         lat=map_dbl(geometry, ~st_centroid(.x)[[2]])) # add centroid values for labels

# SET UP COLOR AND RANGE --------------------------------------------------

# add color palette: 
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#CFCFCF")

# get range of lat/longs from for mapping
mapRange <- c(range(st_coordinates(counties)[,1]),range(st_coordinates(counties)[,2]))

# save data out for map:
save(bgSTs, CA, counties, rb_range, 
     #lakes, rivers,
     annot_sf, shaff_sf, mapRange, cbbPalette, file = "data_output/rangewide_map_data.rda")

# FIRST MAP OF SHAFF VS ANNOT SITES ---------------------------------------

load("data_output/rangewide_map_data.rda")
#rb_diss <- st_read("data_output/rabo_range_simple.shp")

# resave:
# save(bgSTs, CA, counties, lakes, rivers, rb_range, rb_diss,
#        annot_sf, shaff_sf, mapRange, cbbPalette, file = "data_output/rangewide_map_data.rda")


# plot
ggplot() + 
  geom_sf(data=bgSTs, fill="gray90", col="gray50", lty=1.2) +
  geom_sf(data=CA, fill="gray20") + xlab("")+
  geom_sf(data=counties, col="gray50", alpha=0.9) + ylab("") +
  #geom_sf(data=lakes, col="deepskyblue4", fill="deepskyblue4") +
  #geom_sf(data=rivers, col="darkblue", alpha=0.5) +
  geom_sf(data=rb_range, col="orange", fill="orange", alpha=0.7) +
  geom_sf(data=annot_sf, aes(fill="black"), alpha=0.7,size=2.6, pch=21, show.legend = 'point') +
  ggrepel::geom_text_repel(data=annot_sf, aes(x=lon, y=lat, label=siteID), size=1.7, segment.size = .25) +
  geom_sf(data=shaff_sf, aes(fill="white"), size=1.4, pch=21, alpha=0.8, show.legend = 'point') +
  #ggforce::geom_mark_ellipsis(data=annot_sf, aes(y=lat, x=lon, group=admix_orig), col="black") +
  scale_fill_manual(name = 'Localities', 
                    values =c('white'='white', 'black'='black'), 
                    labels = c( 'Peek et al.', 'McCartney-Melstad et al. 2018' )) +
  theme_bw(base_family = "Roboto Condensed") +
  # remove graticule and rotate x axis labels
  theme(panel.grid.major = element_line(color = 'transparent'),
        panel.background = element_rect(fill="darkslategray4"),
        axis.text.x = element_text(angle = 45, vjust = .75),
        legend.justification = c(0.1, 0.1),
        legend.position = c(0.55, 0.5),
        legend.title = element_blank(),
        legend.key.height=unit(1,"line"),
        legend.key.width = unit(1,"line"),
        legend.key = element_blank(),
        legend.background = element_rect(fill = "white")) +
  coord_sf(xlim = mapRange[c(1:2)], ylim = mapRange[c(3:4)]) +
  # add north arrow
  north(x.min = -124.2, x.max = -122.2,
        y.min = 33, y.max = 34.5,
        location = "bottomleft", scale = 0.5) +
  # add scale bar
  scalebar(location = "bottomleft", dist = 200,
           dd2km = TRUE, model = 'WGS84',           
           x.min = -124, x.max = -121,
           y.min = 32.4, y.max = 33.7, height = .15,
           st.size = 2.5, st.dist = .2)


ggsave(filename = paste0("figs/maps_", bamfile, "_revrange_localities.png"), width = 8, height = 11, 
              units = "in", dpi = 300)

# FIG OF ADMIX GROUPS -----------------------------------------------------

# plot
ggplot() + 
  geom_sf(data=bgSTs, fill="gray90", col="gray50", lty=1.2) +
  geom_sf(data=CA, fill="gray20") + xlab("") +
  geom_sf(data=counties, col="gray50", alpha=0.9) + ylab("") +
  geom_sf(data=rb_range, col="gray40", fill="gray50", alpha=0.6) +
  #geom_sf(data=shaff_sf, fill="white", col="gray30", size=1.2, pch=21, show.legend = 'point') +
  ggrepel::geom_text_repel(data=annot_sf, aes(x=lon, y=lat, label=siteID), size=1.7, segment.size = .25) +
  geom_sf(data=annot_sf, aes(fill=admix_groups), size=2, pch=21, show.legend = 'point') +
  #coord_sf(datum=sf::st_crs(4326), ndiscr = 5) + # include for graticule
  scale_fill_manual("Groups",
                     values = c("East"=cbbPalette[1], # E
                                "North-East"=cbbPalette[2], # NE
                                "North-West"=cbbPalette[3], # NW
                                "North-Feather"=cbbPalette[4], #N Feather
                                "West"=cbbPalette[5],  # W
                                "South-West"=cbbPalette[6]),# SW
                     #"Unknown"=cbbPalette[9]), # SFAmerican
                     labels = c("S. Sierra (E)", # E
                                "N. Sierra (NE)", # NE
                                "N. Sierra-Feather", #N Feather
                                "N. Coast (NW)", # NW
                                "S. Coast (W)",  # W
                                #"Unknown",
                                "C. Coast (SW)")) +
  # scale_fill_manual(values = c("East"=cbbPalette[1], 
  #                               "North-East"=cbbPalette[2], 
  #                               "North-West"=cbbPalette[3],
  #                               "North-Feather"=cbbPalette[4],
  #                               "West"=cbbPalette[5], 
  #                               "South-West"=cbbPalette[6])) +
  theme_bw(base_family = "Roboto Condensed") + 
  # remove graticule and rotate x axis labels
  theme(panel.grid.major = element_line(color = 'transparent'),
        panel.background = element_rect(fill="darkslategray4"),
        axis.text.x = element_text(angle = 45, vjust = .75),
        legend.justification = c(0.1, 0.1),
        legend.position = c(0.6, 0.5),
        legend.title = element_blank(),
        legend.key.height=unit(1,"line"),
        legend.key.width = unit(1,"line"),
        legend.key = element_blank(),
        legend.background = element_rect(fill = "white")) +
  coord_sf(xlim = mapRange[c(1:2)], ylim = mapRange[c(3:4)]) +
  # add north arrow
  north(x.min = -124.2, x.max = -122.2,
        y.min = 33, y.max = 34.5,
        location = "bottomleft", scale = 0.5) +
  # add scale bar
  scalebar(location = "bottomleft", dist = 200,
           dd2km = TRUE, model = 'WGS84',           
           x.min = -124, x.max = -121,
           y.min = 32.4, y.max = 33.7, height = .15,
           st.size = 2.5, st.dist = .2)

# with site ID's on map
ggsave(filename = paste0("figs/maps_", bamfile, "_clades.png"), width = 8, height = 11, 
       units = "in", dpi = 300)
# without site ID's on map
ggsave(filename = paste0("figs/maps_", bamfile, "_clades_localities.png"), width = 8, height = 11, 
       units = "in", dpi = 300)


mapview(annot_sf, zcol="admix_groups", layer.name=c("Genetic Groups")) + 
  mapview(rb_range, col.regions="orange", alpha.regions=0.7) + 
  mapview(counties, col.regions="gray")# %>% addMouseCoordinates()

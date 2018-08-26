# Figure 1 Map.
library(tidyverse)
library(USAboundaries)
library(purrr)
library(sf)
library(ggsn)
library(ggrepel)

# get name
bamfile <- "all_rabo_filt_100k"

# get shp file with data
annot_sf <- st_read(paste0("data_output/sites_", bamfile, ".shp"))
st_crs(annot_sf)

# filter to specific localities and add siteID number
annot_sf <- annot_sf %>% distinct(Localty, .keep_all = T) %>% 
  arrange(admx_gr, EcoRegn, River) %>% 
  mutate(localityID=row_number())

# get shp of range
rb_range <- st_read(unzip("data/Rb_Potential_Range_CAandOR.zip")) 
file.remove(list.files(pattern="Rb_Potential_Range_CAandOR", recursive = F))
rb_range <- rb_range %>% st_transform(crs = 4326)
st_crs(rb_range)
rb_range_simple <- st_simplify(rb_range, dTolerance = 0.05)
#plot(rb_range_simple$geom, col = "orange")

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

# add color palette: 
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# get range of lat/longs from for mapping
mapRange <- c(range(st_coordinates(counties)[,1]),range(st_coordinates(counties)[,2]))

# plot
ggplot() + 
  geom_sf(data=bgSTs, fill="gray90", col="gray50", lty=1.2) +
  geom_sf(data=CA, fill="gray20") + 
  geom_sf(data=counties, col="gray50", alpha=0.9) + ylab("") +
  geom_sf(data=rb_range_simple, fill="orange", alpha=0.5)+
  geom_sf(data=annot_sf, col="gray50", fill="black", size=1.5, pch=21, show.legend = 'point') +
  geom_text_repel(data=annot_sf, aes(x = lon, y=lat, label=localityID), segment.size = 0.2, size=2)+
  #coord_sf(datum=sf::st_crs(4326), ndiscr = 5) + # include for graticule
  
  theme_bw(base_family = "Roboto Condensed") +
  # remove graticule and rotate x axis labels
  theme(panel.grid.major = element_line(color = 'transparent'),
        panel.background = element_rect(fill="darkslategray4"),
        axis.text.x = element_text(angle = 45, vjust = .75)) +
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


ggsave(filename = paste0("figs/fig1_maps_", bamfile, "_sites_range.png"), width = 8, height = 11, 
       units = "in", dpi = 300)




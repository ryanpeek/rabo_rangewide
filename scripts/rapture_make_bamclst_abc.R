# RAPTURE bamlist & bamclst creation for different groups/sites
# 2017-Jun

# This script can be used to generate bamlists and bamlist_clst files for different subsampled groups.
# The following pieces should have already occurred and the subsampled bamlists should already be locally 
# available in a "data_output" folder within your R Project.

# 1. USE SFTP TO GRAB FILE ------------------------------------------------

# cd to location where your RProject/data_output folder lives
# sftp -P 2022 USERNAME@farm.cse.ucdavis.edu
# cd to location where the subsampled bamlist lieves (e.g., cd projects/rana_rapture/fastq)
# use GET to pull files from cluster/farm to your local drive and PUT for opposite
# get bamlist_flt_100k

# 2. LOAD LIBRARIES ----------------------------------------------------------

suppressMessages({
  library(tidyverse);
  library(lubridate);
  library(magrittr)
})

options(scipen = 12) # to avoid issues with paste functions, & joining on subsample number

# 3. GET DATA ----------------------------------------------------------------

# set subsample number (so need corresponding subsampled bamlist locally, see Steps 1A-1C)
bamNo<-25

# set site name (will be appended into filename)
#site <- "fea_rasi"
#site <- "fea_rabo"
site <- "rub_rasi"
#site <- "rub_rabo"
#site <- "yub_rasi"
#site <- "yub_rabo"
#site <- "yub_rabo_all"
#site <- "rubyubfea"


# FIX METADATA ------------------------------------------------------------

# METADATA
metadat<- read_csv("data_output/rapture06_metadata_revised.csv") %>% arrange(Seq)

# need to change Seq column to match the new merged one:
#metadat$Seq <- gsub(pattern = "SOMM165", replacement = "SOMM163", x = metadat$Seq)

# drop unneeded cols:
metadat <- select(metadat, -c(Regulated:HUC_12))

# filter out WESTUS Samples and all non-Sierra Samples:
metadat <- filter(metadat, EcoRegion=="Sierra Nevada" | EcoRegion=="Sierra/Basin Range") %>%
  filter(!HU_6_NAME=="San Joaquin")

# split out RIVER as separate col:
metadat <- metadat %>%
  separate(Locality, c("River", "Site"), "-", remove = FALSE)

# add combined SPP col:
metadat <- metadat %>%
  mutate(SPP_ID2 = case_when(
    SPP_ID=="RABO" ~ "RABO",
    SPP_pc1 == "RABO" ~ "RABO",
    SPP_ID == "RASI" ~ "RASI",
    SPP_pc1 == "RASI" ~ "RASI",
    TRUE ~ "RANA"
  ))

metadat %>% group_by(River) %>% tally

metadat <- metadat %>%
  mutate(watershed= case_when(
    River == "BEAR" ~ "Bear",
    grepl("FEA|NFF|MFF", River)  ~ "Feather",
    grepl("NFY|MFY|SFY|DEER|FORD", River)  ~ "Yuba",
    grepl("NFA|MFA|RUB|NFMFA|SFA", River) ~ "American")
  )

metadat %>% group_by(watershed, SPP_ID2) %>% tally

#save(metadat, file = "data_output/rapture06_metadata_sierra.rda")


# 3b. SUBSAMPLED BAMLIST --------------------------------------------------

# this is the subsampled bamlist: e.g., "bamlist_flt_50k"
# modify per user's naming convention

bams <- read_tsv(paste0("data_output/bamlists/bamlist_flt_mrg_",bamNo,"k"), col_names = F)

# remove the *000 component for join, requires fixing scipen for digits
subsamp<-bamNo*1000

# bamlist
bams$X1<-gsub(pattern = paste0(".sortflt.mrg_",subsamp,".bam"), replacement = "", bams$X1)

# 4. FILTER BY: -----------------------------------------------------------

set.seed(111)

# filter out hybrids
datHYB <- metadat %>% filter(grepl("RAP1745|RAP1587", SampleID))

# any number of filtering options can occur...below a few different versions
#datRANA <- metadat %>% filter(SPP_ID2=="RANA")

metadat %>% filter(!SPP_ID2=="RANA") %>% 
  group_by(watershed, SPP_ID2) %>% 
  tally

# general list of sites
metadat %>% filter(!(SPP_ID2=="RANA"), !watershed=="Bear") %>% 
  filter(!(grepl("^MFA-|^NFA|SFA-CAMI|SFY-HUMBUG", Locality)), 
         !(grepl("FEA-SFRockCk|^NFF-Poe|^NFF-EBNFF|^FEA-EBNFF", Locality)),
         !(grepl("RAP1745|RAP1587", SampleID)),
         !(grepl("RAP-122", SampleID)), # weird older sample
         !(grepl("Bean45", LabID)),
         !(grepl("RUB-Zitella", Locality)), # Rub RASI
         !(grepl("^NFY|^MFY|^DEER|^SFY-Scotc",Locality)))  %>% 
  group_by(watershed, SPP_ID2) %>%
  tally

# ALL SIERRAS: sample even number in each watershed and of each species:
# dat <- metadat %>% filter(!(SPP_ID2=="RANA"), !watershed=="Bear") %>% 
#   filter(!(grepl("^MFA-|^NFA|SFA-CAMI|SFY-HUMBUG", Locality)), # american
#          !(grepl("FEA-SFRockCk|^NFF-Poe|^NFF-EBNFF|^FEA-EBNFF", Locality)), # Feather
#          !(grepl("RAP1745|RAP1587", SampleID)), # Hybrids
#          !(grepl("RAP-122", SampleID)), # weird older sample
#          !(grepl("RUB-Zitella", Locality)), # Yuba RASI
#          !(grepl("FORD-NorthCkTrib|FORD-Mossy-P1|FORD-Mossy-P2|FORD-Mossy-P3", Locality)), # Yuba RASI
#          !(grepl("Bean45", LabID)), # Hybrids
#          !(grepl("^NFY|^MFY|^DEER|^SFY-Scotc",Locality))) %>%  # all other sites
#   group_by(watershed, SPP_ID2) %>% 
#   sample_n(size = 15)

# RUB_RASI
dat <- metadat %>% filter(!(SPP_ID2=="RANA"), !watershed=="Bear") %>% 
  filter(!(grepl("^MFA-|^NFA|SFA-CAMI|SFY-HUMBUG", Locality)), # american
         !(grepl("FEA-SFRockCk|^NFF-Poe|^NFF-EBNFF|^FEA-EBNFF", Locality)), # Feather
         !(grepl("RAP1745|RAP1587", SampleID)), # Hybrids
         !(grepl("RAP-122|RAP678", SampleID)), # weird older sample
         !(grepl("RUB-Zitella", Locality)), # Yuba RASI
         !(grepl("FORD-NorthCkTrib|FORD-Mossy-P1|FORD-Mossy-P2|FORD-Mossy-P3", Locality)), # Yuba RASI
         !(grepl("Bean45", LabID)), # Hybrids
         !(grepl("^NFY|^MFY|^DEER|^SFY-Scotc|^SFY-ShadyCk|^SFY-Spring",Locality))) %>%  # all other sites
  filter(watershed=="American", SPP_ID2=="RASI") #%>%


# YUB_RABO_v4
dat <- metadat %>% filter(!(SPP_ID2=="RANA"), !watershed=="Bear") %>% 
  filter(!(grepl("^MFA-|^NFA|SFA-CAMI|SFY-HUMBUG", Locality)), # american
         !(grepl("FEA-SFRockCk|^NFF-Poe|^NFF-EBNFF|^FEA-EBNFF", Locality)), # Feather
         !(grepl("RAP1745|RAP1587", SampleID)), # Hybrids
         !(grepl("RAP-122|RAP678", SampleID)), # weird older sample
         !(grepl("RUB-Zitella", Locality)), # Yuba RASI
         !(grepl("FORD-NorthCkTrib|FORD-Mossy-P1|FORD-Mossy-P2|FORD-Mossy-P3", Locality)), # Yuba RASI
         !(grepl("Bean45", LabID)), # Hybrids
         !(grepl("^NFY|^MFY|^DEER|^SFY-Scotc|^SFY-ShadyCk|^SFY-Spring",Locality))) %>%  # all other sites
  filter(watershed=="Yuba", SPP_ID2=="RABO") %>%
  filter(SampleID %in% c('RAP644','RAP643','RAP641','RAP640','RAP647',
                         'AAA824', 'AAA825', 'AAA826','AAA827','AAA828',
                         'RAP649', 'RAP676', 'RAP677', 'RAP679', 'RAP680'))

dat %>% group_by(Locality, SPP_ID2) %>% tally

#dat <- dat %>% sample_n(size=16)

# ALL_YUB_RABO ONLY
# dat <- metadat %>% filter(!(SPP_ID2=="RANA"), !watershed=="Bear") %>% 
#   filter(!(grepl("^MFA-|^NFA|SFA-CAMI|SFY-HUMBUG", Locality)), # american
#          !(grepl("FEA-SFRockCk|^NFF-Poe|^NFF-EBNFF|^FEA-EBNFF", Locality)), # Feather
#          !(grepl("RAP1745|RAP1587", SampleID)), # Hybrids
#          !(grepl("RAP-122", SampleID)), # weird older sample
#          !(grepl("RUB-Zitella", Locality)), # Yuba RASI
#          !(grepl("FORD-NorthCkTrib|FORD-Mossy-P1|FORD-Mossy-P2|FORD-Mossy-P3", Locality)), # Yuba RASI
#          !(grepl("Bean45", LabID)), # Hybrids
#          !(grepl("^DEER|^SFY-Scotc|^SFY-ShadyCk",Locality))) %>%  # all other sites
#   filter(watershed=="Yuba", SPP_ID2=="RABO") #%>% 
  #group_by(Locality, SPP_ID2) %>% tally

# add the RANA back in:
#dat <- bind_rows(dat, datRANA)

# # double check:
# dat %>% group_by(watershed, SPP_ID2) %>% tally
# 
# # split out by watershed
# datRUBrabo <- dat %>% filter(watershed=="American", SPP_ID2=="RABO")
# datRUBrasi <- dat %>% filter(watershed=="American", SPP_ID2=="RASI")
# datYUBrabo <- dat %>% filter(watershed=="Yuba", SPP_ID2=="RABO")
# datYUBrasi <- dat %>% filter(watershed=="Yuba", SPP_ID2=="RASI")
# datFEArabo <- dat %>% filter(watershed=="Feather", SPP_ID2=="RABO")
# datFEArasi <- dat %>% filter(watershed=="Feather", SPP_ID2=="RASI")

# 5. JOIN WITH FLT SUBSAMPLE LIST -----------------------------------------

# check col names for join...should be BAMFILE name with plate/well ID
dfout <- inner_join(dat, bams, by=c("Seq"="X1")) %>% arrange(Seq)

# JOIN WITH BAMLIST TO DOUBLE CHECK
bams2 <- read.table("./data_output/bamlists/rub_rasi_25k_thresh.bamlist", stringsAsFactors = F, header = F) # read in the bamlist
bams2$V2 <- sub('\\..*$', '', basename(bams2$V1)) # remove the path and file extension
dfout2 <- inner_join(dfout, bams2, by=c("Seq"="V2")) %>% select(-V1) # join with the

# add MH014
dfout2 <- dfout2 %>% bind_rows(., filter(dfout, SampleID=="MH014")) %>% arrange(Seq)

# dfout %>% group_by(watershed, SPP_ID2) %>% tally
# dfout %>% group_by(Locality, SPP_ID2) %>% tally

# write to file for later use/filtering
#write_rds(x = dfout, path = paste0("data_output/bamlist_dat_", site,"_",bamNo, "k_thresh.rds"))

# individual subpops

# # make one for RUB_rabo
# dfout1 <- inner_join(datRUBrabo, bams, by=c("Seq"="X1")) %>% arrange(Seq)
# # make one for RUB_rasi
# dfout2 <- inner_join(datRUBrasi, bams, by=c("Seq"="X1")) %>% arrange(Seq)
# # make one for YUB_rabo
# dfout3 <- inner_join(datYUBrabo, bams, by=c("Seq"="X1")) %>% arrange(Seq)
# # make one for YUB_raso
# dfout4 <- inner_join(datYUBrasi, bams, by=c("Seq"="X1")) %>% arrange(Seq)
# # make one for FEA_rabo
# dfout5 <- inner_join(datFEArabo, bams, by=c("Seq"="X1")) %>% arrange(Seq)
# # make one for FEA_rasi
# dfout6 <- inner_join(datFEArasi, bams, by=c("Seq"="X1")) %>% arrange(Seq)

# 6. WRITE TO BAMLIST PURRR -----------------------------------------------------

site <- "rubyubfea"

site1 <- "rub_rabo"
site2 <- "rub_rasi"
site3 <- "yub_rabo"
site4 <- "yub_rasi"
site5 <- "fea_rabo"
site6 <- "fea_rasi"

sites <- list(site1,site2,site3,site4,site5,site6)
dfouts <- list(dfout1,dfout2, dfout3, dfout4, dfout5, dfout6)

# use purrr package and map2 to work over two lists for final bamlist
map2(dfouts, sites, ~ write_delim(as.data.frame(
  paste0("/home/rapeek/projects/rasi_hybrid/alignments/", .x$Seq, ".sortflt.mrg.bam")),
  path = paste0("data_output/bamlists/",.y,"_",bamNo,"k_thresh.bamlist"), col_names = F))


# 6. WRITE TO BAMLIST SINGLE -----------------------------------------------------

# Write to bamlist for angsd call (w subsample appended)
#write_delim(as.data.frame(paste0("/home/rapeek/projects/rasi_hybrid/alignments/", dfout$Seq, ".sortflt.mrg_",subsamp,".bam")),path = paste0("data_output/bamlists/",site,"_",bamNo,"k.bamlist"), col_names = F)

# Write to bamlist for angsd call (for IBS, no bamNo)
write_delim(as.data.frame(paste0("/home/rapeek/projects/rasi_hybrid/alignments/", dfout2$Seq, ".sortflt.mrg.bam")), path = paste0("data_output/bamlists/",site,"_",bamNo,"k_thresh.bamlist"), col_names = F)

# write to file for later use/filtering
#write_rds(x = dfout, path = paste0("data_output/bamlist_dat_",site,"_",bamNo, "k_thresh.rds"))

# 7. WRITE CLST FILE ---------------------------------------------------------

names(dfout)

clst_outs <- map(dfouts, ~ .x %>% as.data.frame() %>%   
                   dplyr::select(watershed, SampleID, Locality) %>%
                   dplyr::rename(FID=watershed, IID=SampleID, CLUSTER=Locality))
map2(clst_outs, sites, ~ write_delim(.x, path=paste0("data_output/bamlists/",.y,"_",bamNo,"k.clst")))
       

# FINE SCALE
# clst_out<-dfout %>%
#   dplyr::select(watershed, SampleID, Locality) %>%
#   dplyr::rename(FID=watershed, IID=SampleID, CLUSTER=Locality)
# length(unique(clst_out$CLUSTER))
# length(unique(clst_out$FID))
# 
# # count NAs in your clst file
# clst_out %>% filter(is.na(FID)) %>% tally
# clst_out %>% filter(is.na(CLUSTER)) %>% tally
# 
# # view list
# clst_out
# 
# write_delim(clst_out, path=paste0("data_output/bamlists/bamlist_mrg_",site1,"_",bamNo,"k_clst"))

# 8. TERMINAL SFTP --------------------------------------------------------

# farmer (sftp)
# cd projects/rasi_hybrid/bamlists/subpops
# cd projects/rasi_hybrid/pop_gen
paste0("put ",site,"_*",bamNo,"k*.bamlist") # (this goes from local to cluster)

# 9. BASH: PCA CALC SITES --------------------------------------------------

# Use angsd to run pca_calc_sites script: 

# create command:
lsite<- tolower(site)

# NEW IBS METHOD

# less bamlist_mrg_RASI_all_50k | sed 's/\_50000//g' > bamlist_mrg_RASI_all_50k_thresh

paste0("sbatch -p high -t 48:00:00 --mail-type ALL --mail-user rapeek@ucdavis.edu scripts/03a_pca_new.sh bamlist_mrg_",site,"_",bamNo,"k_thresh", " ", lsite, "_",bamNo,"k_thresh", " bamlists/subpops")


# 10. BASH: MAKE PCA PLOTS --------------------------------------------------

# cd projects/rana_rapture/fastq
## (sbatch or srun)

lsite<- tolower(site)
paste0("srun -p med -t 24:00:00 04_pca_plot.sh ", lsite, "_",bamNo,"k", " ","bamlist_",site,"_",bamNo,"k_clst")

# sbatch -p med -t 12:00:00 04_pca_plot.sh RABO_nfa_75k bamlist_RABO_nfa_75k_clst

# can check datestamp & sort by time: ls -lt FILE*


# 11. SFTP and GET PDFS ---------------------------------------------------

# then sftp and "get" to get pdfs back to your computer (see step 1C).

# cd Documents/github/rana_rapture/data_output/bamlists
# farmer (sftp)
# cd projects/rana_rapture/fastq/pca_pdfs

# get rabo_nfa_75k* # (this goes from cluster to local)
paste0(lsite,"_",bamNo,"k*")


# MAPS: STATIC ----------------------------------------------------

# map_data("county") %>% 
#   filter(region %in% c("california","oregon")) %>%
#   ggplot() +
#   geom_polygon(aes(x = long, y = lat, group = group)) + 
#   coord_map("gilbert") +
#   geom_point(data=dfout, aes(x=lon, y=lat, fill=HU_6_NAME), size=4, pch=21)+
#   scale_fill_viridis_d("HUC") + 
#   ggtitle(paste0(site, ": ", bamNo,"k (n=",nrow(dfout),")"))

# MAPS: LEAFLET ----------------------------------------

## if you have X/Y for data, look at a map of samples from subsample

# library(sf)
# library(leaflet)
# 
# # hydrology (HUCS)
# h8 <- st_read("data_output/shps/HUC8_named_westcoast.shp", quiet = T)
# h8 <- st_transform(h8, crs = 4326)
# h6 <- st_read("data_output/shps/WBD_HUC6.shp", quiet = T)
# h6 <- st_transform(h6, crs = 4326)
# 
# # make map
# subMAP <- leaflet() %>%
#   addTiles(group = "Basemap") %>%
#   #setView(lng = -120.8, lat = 39, zoom = 5) %>%  # set to Auburn/Colfax, zoom 5 for CA
#   addProviderTiles("Stamen.TopOSMFeatures", group = "OSM Features") %>%
#   addProviderTiles("Esri.WorldTopoMap", group = "Topo") %>%
#   addProviderTiles("Esri.WorldImagery", group = "ESRI Aerial") %>%
# 
#   # huc6
#   addPolygons(data=h6, group="HUC6", color="#473C8B", weight = 1.5,
#               fillColor = "transparent", label = ~paste0(HU_6_Name, ": " ,HUC_6)) %>%
#   hideGroup("HUC6") %>%
#   # huc8
#   addPolygons(data=h8, group="HUC8", color="darkblue", weight = 1.3,
#               fillColor = "transparent", label = ~HU_8_NAME) %>%
#   hideGroup("HUC8") %>%
# 
#   # add scale bar
#   addMeasure(position = "topright",
#              primaryLengthUnit = "kilometers",
#              primaryAreaUnit = "hectares",
#              secondaryAreaUnit = "sqmiles",
#              activeColor = "#3D535D",
#              completedColor = "#7D4479") %>%
# 
#   # add fancy zoom button to recenter
#   # addEasyButton(easyButton(
#   #   icon="fa-globe", title="Zoom to Level 5",
#   #   onClick=JS("function(btn, map){ map.setZoom(5); }"))) %>%
# 
#   # add subsamples
#   addCircleMarkers(data=dfout, group="Samples",
#                    clusterOptions = markerClusterOptions(),
#                    lng = ~lon, lat=~lat,opacity = 0.5,
#                    popup=paste0("Locality: ", dfout$Locality, "<br>",
#                                 "HUC12Name: ", dfout$HU_12_NAME,
#                                 "<br>","SampleID: ",dfout$SampleID,
#                                 "<br>", "SPP_ID: ",dfout$SPP_ID,
#                                 "<br>", "Elev (m): ", dfout$elev_m, "<br>",
#                                 "StreamName: ", dfout$GNIS_NAME),
#                    weight=0.6,radius=10, stroke=TRUE,
#                    fillColor = ifelse(dfout$SPP_ID=="RABO" |
#                                         dfout$SPP_pc1=="RABO", "yellow", "dodgerblue")) %>%
# 
#   # add layer/legend control
#   addLayersControl(
#     baseGroups = c("Basemap", "Topo", "ESRI Aerial", "OSM Features"),
#     overlayGroups = c("HUC6", "HUC8", "Samples"),
#     options = layersControlOptions(collapsed = T))
# 
# subMAP





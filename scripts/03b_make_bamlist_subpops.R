# RAPTURE bamlist creation for different groups/sites
# 2018-May

# This script can be used to generate bamlists for different analyses. Typically this requires a few components:
## A data_output folder
## Raw bamlist of all bams from cluster
## Can use subsampled data as a threshold, here use 30000 reads as threshold:
## e.g., 'sbatch -p high -t 24:00:00 run_subsample.sh bam_sort_list 30000'
## Then create bamlist based on that subsample value: 'ls *flt_100000* > bamlist_flt_100k'

# 01. LOAD LIBRARIES & Global Settings ---------------------------

suppressPackageStartupMessages(library(tidyverse))

# Global Setting: set digits to 12 to avoid paste function/joining issues
options(scipen = 12)
set.seed(111) # for repeatable random sampling

# set the reads threshold (number of minimum reads subsampled
bamNo<-25

# set site name (will be appended into filename)
#site <- "ncoast_rabo_n9" # all_rabo, all_rabo_filt01

# 02. LOAD METADATA --------------------------------------------

metadat <- read_rds(path = "data_output/rapture_metadata_rabo_quant.rds")

# view summary of data:
metadat %>% group_by(HUC_8) %>% tally %>% print(n=Inf)

# check for duplicates:
metadat[duplicated(metadat$Seq),] %>% arrange(SampleID) %>% tally()

# 03. READ FULL BAMLIST --------------------------------------------------

# this can be a full bamlist (all bams) or subsampled bamlist

bams <- read_tsv(paste0("data_output/bamlists/bamlist_flt_mrg_",bamNo,"k"), col_names = F)

# remove the *000 component for join, requires fixing scipen for digits
subsamp<-bamNo*1000

# Fix the name in bamlist (based on subsample output)
bams$X1<-gsub(pattern = paste0(".sortflt.mrg_",subsamp,".bam"), replacement = "", bams$X1)

# 04. RANGEWIDE Filter Outliers ---------------------------------------------------------

filt01 <- c("RAP-092","RAP-104","RAP-122", "RAP-226", "RAP-278", "RAP-335",
           "RAP-346", "RAP-347", "RAP-348", "RAP-1649",
           # hybrids:
           "RAP1745", "RAP1587")

# NO FILTER
# dat <- metadat # all samples, NO FILTER
# view and sort
# summary(as.factor(metadat$HU_8_NAME))
# dat %>% group_by(River) %>% tally %>% arrange(n) %>%  print(n=Inf)

# FILTER out outlier samples (after first run)
dat_filt <- filter(metadat, !SampleID %in% filt01)

# see data
dat_filt %>% group_by(River) %>% tally %>% arrange(n) %>%  print(n=Inf) # river
dat_filt %>% group_by(River,Site) %>% tally %>% arrange(n) %>%  View() # Site

# filter to sites with more than 1 sample
dat_filt <- dat_filt %>% group_by(River, Site) %>% add_count(Site) %>% rename(n_site=n) %>%
  arrange(n_site, River) %>% # View
  filter(n_site>1)  # keep only sites with more than X samples

dat_filt %>% group_by(River, EcoRegion) %>% tally %>% arrange(n) %>%  View()
dat_filt %>% group_by(EcoRegion, River, Site) %>% tally %>% filter(n > 10) %>% arrange(n) %>%  View()

# 05a. JOIN WITH FLT SUBSAMPLE LIST ----------------------------------------

# check col names for join...should be BAMFILE name with plate/well ID
dfout <- inner_join(dat_filt, bams, by=c("Seq"="X1")) %>% arrange(Seq) %>% 
  select(-n_site) %>% 
  group_by(River, Site) %>% 
  add_count(Locality) %>% 
  rename(n_site=n)

# look at counts by site and watershed
# dfout %>% filter(SPP_ID=="RABO") %>% group_by(Locality) %>% tally %>% as.data.frame

# filter localities with less than 5 samples
# dfout <- dfout %>% add_count(Locality) %>% rename(n_local=n) %>% 
#  filter(n_local>4) 

# dfout %>% filter(SPP_ID=="RABO") %>% group_by(Locality) %>% tally %>% 
#  arrange(n) %>% as.data.frame

# check for dups
dfout[duplicated(dfout$Seq),] %>% arrange(SampleID) %>% tally()

# 05b. SUBSAMPLE SITES AND SPLIT INTO LIST --------------------------------

# split out sites with less than 9 samples (to rejoin later):
dfout_low_n <- dfout %>% filter(n_site < 11, n_site>2)
dfout_low_n %>% filter(SPP_ID=="RABO") %>% group_by(Locality) %>% tally %>% 
  arrange(n) %>% as.data.frame

dfout_low_n %>% group_by(River, Site) %>% tally %>% View()

# now sample down to 10 per Locality :
dfout_sampled <- dfout %>% filter(n_site > 10) %>% group_by(Locality) %>% 
  sample_n(10, replace = F)

# rebind
dfout_bind <- bind_rows(dfout_low_n, dfout_sampled)

dfout_bind %>% group_by(EcoRegion, Locality) %>% tally %>% View()

# split into list of dataframes by Locality 
sites <- dfout_bind %>% split(.$Locality)

# 05c. MAKE QUICK MAP ----------------------------------------------------------

library(mapview)
library(sf)

# jitter coords slightly for viewing
dfout_sf <- dfout_bind
dfout_sf$lat_jitter <- jitter(dfout_sf$lat, factor = 0.001)
dfout_sf$lon_jitter <- jitter(dfout_sf$lon, factor = 0.01)
dfout_sf <- dfout_sf %>% filter(!is.na(lat)) %>%  
  distinct(lat,.keep_all = T) %>% 
  select(River, Site, SampleID, LabID, SPP_ID, lat, lon, elev_m:Locality,Locality_details, EcoRegion, n_site:lon_jitter)

# make sf:
dfout_sf <- st_as_sf(dfout_sf, 
                     coords = c("lon_jitter", "lat_jitter"), 
                     remove = F, # don't remove these lat/lon cols from df
                     crs = 4326) # add projection

st_write(dfout_sf, "data_output/sites_25k_n3.shp", delete_dsn = T)

mapview(dfout_sf) %>% addMouseCoordinates()


# 06. WRITE TO BAMLISTS WITH PURRR ---------------------------------------------

# use purrr package and map2 to work over two lists for final bamlist

# first make list of sites for use in file_names
sitenames <- as.list(c(tolower(names(sites))))
# write out list of sites
subpops_list_n2 <- tolower(names(sites))
write_lines(subpops_list_n2, path = "data_output/bamlists/subpops_list_n2")


map2(sites, sitenames, ~ write_delim(as.data.frame(
 paste0("/home/rapeek/projects/rangewide/alignments/", .x$Seq, ".sortflt.mrg.bam")),
 path = paste0("data_output/bamlists/subpops/",.y,"_",bamNo,"k_thresh.bamlist"), col_names = F))


# 07. TERMINAL SFTP --------------------------------------------------------

# farmer (sftp)
# cd projects/rangewide/pop_gen/bamlists/subpops

paste0("put ",site,"_*",bamNo,"k*.bamlist") # (this goes from local to cluster)

# 08. BASH: PCA CALC SITES --------------------------------------------------

# Use angsd to run pca_calc_sites script: 

# create command:
lsite<- tolower(site)

# NEW IBS METHOD
map(sitenames, ~ paste0("sbatch -t 2880 --mail-type ALL --mail-user rapeek@ucdavis.edu 03_pca_ibs.sh ",.x,"_",bamNo,"k_thresh.bamlist", " ", .x, "_",bamNo,"k"))









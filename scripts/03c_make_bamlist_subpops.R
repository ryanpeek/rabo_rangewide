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
site <- "ncoast_rabo_n9" # all_rabo, all_rabo_filt01

# 02. LOAD METADATA --------------------------------------------

metadat <- read_csv("data/rapture_metadata_rabo.csv") %>% arrange(Seq)

# view summary of data:
#metadat %>% group_by(HUC_8) %>% tally

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
#dat <- metadat # all samples, NO FILTER

# FILTER out outlier samples (after first run)
dat <- filter(metadat, !SampleID %in% filt01)

# view and sort
summary(as.factor(metadat$HU_8_NAME))
dat %>% group_by(River) %>% tally %>% arrange(n) %>%  print(n=Inf)

# tally by river 
dat_filt <- dat %>% add_count(River) %>% rename(n_river=n) %>% 
  #arrange(River, n_river) %>% 
  filter(n_river>8)  # keep only rivers with more than 8 samples

dat_filt %>% group_by(River) %>% tally %>% arrange(n) %>%  print(n=Inf)

# 05a. JOIN WITH FLT SUBSAMPLE LIST ----------------------------------------

# check col names for join...should be BAMFILE name with plate/well ID
dfout <- inner_join(dat_filt, bams, by=c("Seq"="X1")) %>% arrange(Seq)

# look at counts by site and watershed
dfout %>% filter(SPP_ID=="RABO") %>% group_by(Locality) %>% tally %>% as.data.frame

# filter localities with less than 5 samples
dfout <- dfout %>% add_count(Locality) %>% rename(n_local=n) %>% 
  filter(n_local>4) 

dfout %>% filter(SPP_ID=="RABO") %>% group_by(Locality) %>% tally %>% 
  arrange(n) %>% as.data.frame


# 05b. SUBSAMPLE SITES AND SPLIT INTO LIST --------------------------------

# split out sites with less than 9 samples (to rejoin later):
dfout_low_n <- dfout %>% filter(n_local < 11)
dfout_low_n %>% filter(SPP_ID=="RABO") %>% group_by(Locality) %>% tally %>% 
  arrange(n) %>% as.data.frame

# now sample down to 10 per Locality :
dfout_sampled <- dfout %>% filter(n_local > 10) %>% group_by(Locality) %>% 
  sample_n(10, replace = F)

# rebind
dfout_bind <- bind_rows(dfout_low_n, dfout_sampled)

dfout_bind %>% group_by(Locality) %>% tally %>% as.data.frame

# split into list of dataframes by Locality 
sites <- dfout_bind %>% split(.$Locality)

# 05c. MAKE QUICK MAP ----------------------------------------------------------

library(mapview)
library(sf)

# jitter coords slightly for viewing
dfout_sf <- dfout_bind
dfout_sf$lat_jitter <- jitter(dfout_sf$lat, factor = 0.001)
dfout_sf$lon_jitter <- jitter(dfout_sf$lon, factor = 0.01)
dfout_sf <- dfout_sf %>% filter(!is.na(lat)) %>%  select(-lat, -lon)

# make sf:
dfout_sf <- st_as_sf(dfout_sf, 
                     coords = c("lon_jitter", "lat_jitter"), 
                     remove = F, # don't remove these lat/lon cols from df
                     crs = 4326) # add projection


mapview(dfout_sf) %>% addMouseCoordinates()


# 06. WRITE TO BAMLISTS WITH PURRR ---------------------------------------------

# use purrr package and map2 to work over two lists for final bamlist

# first make list of sites for use in file_names
sitenames <- as.list(c(tolower(names(sites))))

map2(sites, sitenames, ~ write_delim(as.data.frame(
 paste0("/home/rapeek/projects/rangewide/alignments/", .x$Seq, ".sortflt.mrg.bam")),
 path = paste0("data_output/bamlists/",.y,"_",bamNo,"k_thresh.bamlist"), col_names = F))


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









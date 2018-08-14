## Create Bamlists for Analysis
## R. Peek 2018

# Fri May  4 21:30:29 2018 ------------------------------

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
site <- "ncoast_rabo_n4"

# 02. LOAD METADATA --------------------------------------------

metadat <- read_csv("data/rapture_metadata_rabo.csv") %>% arrange(Seq)

# view summary of data:
#metadat %>% group_by(HUC_8) %>% tally
#summary(as.factor(metadat$HUC_8))

# 03. READ FULL BAMLIST --------------------------------------------------

# this can be a full bamlist (all bams) or subsampled bamlist

bams <- read_tsv(paste0("data_output/bamlists/bamlist_flt_mrg_",bamNo,"k"), col_names = F)

# remove the *000 component for join, requires fixing scipen for digits
subsamp<-bamNo*1000

# Fix the name in bamlist (based on subsample output)
bams$X1<-gsub(pattern = paste0(".sortflt.mrg_",subsamp,".bam"), replacement = "", bams$X1)

# 04a. NORTHCOAST ---------------------------------------------------------

# By EcoRgions
dat <- filter(metadat, grepl("^North Coast|North Coast$|^Northern CA Coastal", EcoRegion))
summary(as.factor(dat$HU_8_NAME))
summary(as.factor(dat$EcoRegion))
summary(as.factor(dat$River))

# 05a. JOIN WITH RAW BAMLIST -----------------------------------------

# check col names for join...should be BAMFILE name with plate/well ID
dfout <- inner_join(dat, bams, by=c("Seq"="X1")) %>% arrange(Seq)

# check tally's of groups
dfout %>% group_by(River) %>% 
  arrange(River) %>% tally 
dfout %>% group_by(Locality) %>% tally %>% as.data.frame

# sample 10 from every River where possible, drop singles (GUAL) 
(dfout_subsampled <- dfout %>% group_by(River) %>% 
    arrange(River) %>% # order by river
    nest %>% # collapse data 
    mutate(n=unlist(map(data, tally))) %>% # count number per group
    filter(n >=4) %>%  # filter to a minimum of 4 samples or more
    mutate(samples=map2(data, 4, sample_n, replace = F)) %>% # sample 4 per group
    select(River, samples) %>% # pull out the grouping and sample data only
    unnest()) # IT WORKS!!!!

dfout_subsampled %>% group_by(River) %>% tally


# 05b. MAKE QUICK MAP ----------------------------------------------------------

# load function, requires a bamlist ending in ".bamlist"
#source("scripts/functions/f_map_from_bamlists.R")
#make_map_from_bamlist("fea_rasi", 25)

library(mapview)
library(sf)

# jitter coords slightly for viewing
dfout_sf <- dfout_subsampled
dfout_sf$lat_jitter <- jitter(dfout_sf$lat, factor = 0.001)
dfout_sf$lon_jitter <- jitter(dfout_sf$lon, factor = 0.01)
dfout_sf <- dfout_sf %>% select(-lat, -lon)

# make sf:
dfout_sf <- st_as_sf(dfout_sf, 
                   coords = c("lon_jitter", "lat_jitter"), 
                   remove = F, # don't remove these lat/lon cols from df
                   crs = 4326) # add projection


mapview(dfout_sf)

# 06. WRITE TO BAMLIST SINGLE ----------------------------------------

# Write to bamlist for angsd call (for IBS, no bamNo)
write_delim(as.data.frame(paste0("/home/rapeek/projects/rangewide/alignments/", dfout_subsampled$Seq, ".sortflt.mrg.bam")), path = paste0("data_output/bamlists/",site,"_",bamNo,"k_thresh.bamlist"), col_names = F)

# 07. PUT BAMLIST ON CLUSTER -----------------------------------------

# go to local dir with bamlists
# cd data_output/bamlists/

# farmer (sftp)
# cd projects/rangewide/pop_gen

# get file
paste0("put ",site,"_*",bamNo,"k*.bamlist") # (this goes from local to cluster)

# 08. ANGSD PCA (IBS) ------------------------------------------------

# make sure site is all lowercase
lsite<- tolower(site)

# NEW IBS METHOD
paste0("sbatch -t 2880 --mail-type ALL --mail-user rapeek@ucdavis.edu 03_pca_ibs.sh ",site,"_",bamNo,"k_thresh.bamlist", " ", lsite, "_",bamNo,"k")




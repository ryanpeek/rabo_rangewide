## Create Bamlists for Analysis
## R. Peek 2018

# Mon Aug 13 16:25:48 2018 ------------------------------

# This script can be used to generate bamlists for different analyses. Typically this requires a few components:
## A data_output folder
## Raw bamlist of all bams from cluster
## Can use subsampled data as a threshold, here use 30000 reads as threshold:
## e.g., 'sbatch -p high -t 24:00:00 run_subsample.sh bam_sort_list 30000'
## Then create bamlist based on that subsample value: 'ls *flt_100000* > bamlist_flt_100k'

# 01. LOAD LIBRARIES & Global Settings ---------------------------

suppressPackageStartupMessages(library(tidyverse))
options(tibble.print_min = Inf) # print all rows 
# Global Setting: set digits to 12 to avoid paste function/joining issues
options(scipen = 12)
set.seed(111) # for repeatable random sampling

# SET SITES AND SUBSAMPLE LEVEL -------------------------------------------

# set the reads threshold (number of minimum reads subsampled
bamNo<-100

site <- "all_rabo_filt_100k" 

# rabo_nosocoast_filt10_1_100k: no southwest coast samples, filtered outliers, samples per locality between 1-10
# rabo_nofeath_filt10_1_100k: no feather samples, filtered outliers, samples per locality between 1-10
# all_rabo_filt10_100k: filtered outliers, samples per locality between 2-10
# all_rabo_filt10_1_100k: filtered outliers, samples per locality between 1-10
# all_rabo_filt_100k: filtered outliers, samples per locality > 2
# all_rabo_100k: no filter, no restriction, all samples.



# 02. LOAD METADATA --------------------------------------------

metadat <- read_rds(path = "data_output/rapture_metadata_rabo_quant.rds")

# Fix trin-sftrinity sandybar
unique(metadat$Locality) %>% sort()
metadat$Locality<-tolower(gsub(pattern = "[[:space:]]", replacement = "-", x = metadat$Locality))

# fix deer-clearck/ deer-clec
metadat$Locality <- gsub(pattern="deer-clearck", replacement = "deer-clec", x=metadat$Locality)



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

# check for duplicates:
bams[duplicated(bams$X1),] # yay

# 04a. ALL SAMPLES --------------------------------------------------------

dat <- metadat # for all samples if no filter needed

# 04b. ECOREGIONS ---------------------------------------------------------

# By EcoRgions
# dat %>% group_by(EcoRegion) %>% tally
# 
# # summary(as.factor(metadat$EcoRegion))
# dat <- filter(metadat, !grepl("^Central CA Coastal Foothills|Central Coast$", EcoRegion))
# 
# dat %>% group_by(EcoRegion) %>% tally

# Filter out Feather Sierras/BASIN Range Samples
#dat <- filter(metadat, !EcoRegion=="Sierra/Basin Range")

# be aware outliers: 
## BEAR-MISS Canyon RAP-278
## SFY-SpringCk RAP 122

# 04d. FILTER BY SITE NAME -------------------------------------------------

# NFA
#dat <- filter(metadat, SPP_ID=="RABO", grepl('^NFA', Locality)) 

# MFA
# dat <- filter(metadat, SPP_ID=="RABO", grepl('^MFA|^RUB|^NFMFA|^SFA', Locality))

# SFY
#dat <- metadat[startsWith(x = metadat$Locality, "SFY") & metadat$SPP_ID=="RABO",]

# FEATHER
#dat <- filter(metadat, SPP_ID=="RABO", grepl('^FEA|^NFF', Locality)) 

# BEAR 
#dat <- filter(metadat, SPP_ID=="RABO", grepl('^BEAR', Locality)) 

# VAN DUZEN
#dat <- filter(metadat, SPP_ID=="RABO", grepl('^VANDZ', Locality)) 

#  EEL
#dat <- filter(metadat, SPP_ID=="RABO", grepl('^EEL|SFEEL', Locality)) 

#  TRIN
#dat <- filter(metadat, SPP_ID=="RABO", grepl('^TRIN|KLAM', Locality)) 

#  MFA_RUB
#dat <- metadat[startsWith(x = metadat$Locality, "MFA")| startsWith(x = metadat$Locality, "RUB"),]


# 05. JOIN WITH RAW BAMLIST -----------------------------------------

# check col names for join...should be BAMFILE name with plate/well ID
dfout <- inner_join(dat, bams, by=c("Seq"="X1")) %>% arrange(Seq)

# check for duplicates:
#dfout[duplicated(dfout$Seq),] %>% arrange(SampleID) %>% tally()

# check tally's of groups
#dfout %>% group_by(River) %>% tally
#dfout %>% group_by(Locality) %>% tally

# 05b. Filter out sites with low sample n ---------------------------------

dfout <- dfout %>% group_by(Locality) %>% add_tally() %>% 
  rename(n_locality=n)
#dfout %>% select(Locality, n_locality) %>% View

# filter down to all samples between 2 & 10 per locality:
#dfout <- bind_rows(filter(dfout, n_locality<10, n_locality>1), sample_n(dfout[dfout$n_locality>10,], 10))

# filter down to all samples between 1 & 10 per locality:
dfout <- bind_rows(filter(dfout, n_locality<10), sample_n(dfout[dfout$n_locality>10,], 10))

# filter out localities with less than 3 samples
dfout <- dfout %>% filter(n_locality>2)

# filter out outliers
outliers <- c("RAP-040","RAP-092", "RAP-097", "RAP-104","RAP-122",
            "RAP-177","RAP-226", "RAP-278", "RAP-335",
            "RAP-346", "RAP-347", "RAP-348", "RAP-1649", 
            # hybrids:
            "RAP1745", "RAP1587")
# possible outliers: RAP-357, RAP-039, RAP-420

# check tallys of NA
#dfout %>% filter(is.na(River)) %>% tally
#dfout %>% filter(is.na(Locality)) %>% tally
dfout <- filter(dfout, !SampleID %in% outliers)

# 06. WRITE TO BAMLIST SINGLE ----------------------------------------

# Write to bamlist for angsd call (for IBS, no bamNo)
write_delim(as.data.frame(paste0("/home/rapeek/projects/rangewide/alignments/", dfout$Seq, ".sortflt.mrg.bam")), path = paste0("data_output/bamlists/",site,"_",bamNo,"k_thresh.bamlist"), col_names = F)

# 07. PUT BAMLIST ON CLUSTER -----------------------------------------

# go to local dir with bamlists
# cd data_output/bamlists/

# farmer2 (sftp)
#cd projects/rangewide/pop_gen/bamlists
paste0("put ",site,"_*",bamNo,"k*.bamlist") # (this goes from local to cluster)

# 08. ANGSD PCA (IBS) ------------------------------------------------

# create command:
lsite<- tolower(site)

# NEW IBS METHOD
paste0("sbatch -p high -t 3600 --mem=32G --mail-type ALL --mail-user rapeek@ucdavis.edu 03_pca_ibs.sh ",site,"_",bamNo,"k")


# ADMIX k=12 (sbatch -t 3600 -p med --mem=60G 04_get_admix.sh all_rabo_filt10_100k 12)
paste0("sbatch -p med -t 3600 --mem=60G --mail-type ALL --mail-user rapeek@ucdavis.edu 04_get_admix.sh ", site,"_",bamNo,"k ", 12)



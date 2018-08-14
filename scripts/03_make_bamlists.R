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
bamNo<-25

site <- "all_rabo"

# 02. LOAD METADATA --------------------------------------------

metadat <- read_rds(path = "data_output/rapture_metadata_rabo_quant.rds")

# view summary of data:
metadat %>% group_by(HUC_8) %>% tally %>% print(n=Inf)

# need to make a new field to match the bam names (this is lame but whatever)
metadat <- metadat %>% separate(seqID, into = c("barcode", "wellcode"), drop=T)

# now make the merge field
metadat <- metadat %>% 
  mutate(Seq = paste0("SOMM163_", barcode, "_RA_GG", wellcode, "TGCAGG")) %>% 
  distinct(Seq, .keep_all = T) # without this there are 10 dups, should have 1103 tot

# check for duplicates:
metadat[duplicated(metadat$Seq),] %>% arrange(SampleID) %>% View()


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

# 04b. NORTHCOAST ---------------------------------------------------------

# set site name (will be appended into filename)
site <- "ncoast_rabo"

summary(as.factor(metadat$HU_8_NAME))

# By HUC8 (Wide Range)
# dat<- filter(metadat, HU_8_NAME %in% c("Mad-Redwood", "Mattole", "Lower Eel", "Russion", "Tomales-Drake Bays", "Trinity", "South Fork Trinity", "Smith", "Lower Klamath", "Gualala-Salmon"))

# By EcoRgions

summary(as.factor(metadat$EcoRegion))
dat <- filter(metadat, grepl("^North Coast|North Coast$", EcoRegion))
summary(as.factor(metadat$HU_8_NAME))


# 04c. SIERRAS/BASIN RANGE ------------------------------------------------

# Sierras/BASIN
# dat <- filter(metadat, SPP_ID=="RABO" | SPP_pc1=="RABO") %>%
#   filter(EcoRegion=="Sierra Nevada" | EcoRegion=="Sierra/Basin Range")

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

# check tally's of groups
dfout %>% group_by(River, SPP_ID) %>% tally
dfout %>% group_by(Locality, SPP_ID) %>% tally

# check tallys of NA
#dfout %>% filter(is.na(River)) %>% tally
#dfout %>% filter(is.na(Locality)) %>% tally

# 06. WRITE TO BAMLIST SINGLE ----------------------------------------

# Write to bamlist for angsd call (for IBS, no bamNo)
write_delim(as.data.frame(paste0("/home/rapeek/projects/rangewide/alignments/", dfout$Seq, ".sortflt.mrg.bam")), path = paste0("data_output/bamlists/",site,"_",bamNo,"k_thresh.bamlist"), col_names = F)

# 07. PUT BAMLIST ON CLUSTER -----------------------------------------

# go to local dir with bamlists
# cd data_output/bamlists/

# farmer2 (sftp)
# cd projects/rangewide/pop_gen/bamlists
paste0("put ",site,"_*",bamNo,"k*.bamlist") # (this goes from local to cluster)

# 08. ANGSD PCA (IBS) ------------------------------------------------

# create command:
lsite<- tolower(site)

# NEW IBS METHOD
paste0("sbatch -p high -t 2880 --mail-type ALL --mail-user rapeek@ucdavis.edu 03_pca_ibs.sh ",site,"_",bamNo,"k_thresh.bamlist", " ", lsite, "_",bamNo,"k")




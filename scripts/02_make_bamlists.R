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

# set site name (will be appended into filename)
site <- "ncoast_rabo"

summary(as.factor(metadat$HU_8_NAME))

# By HUC8 (Wide Range)
# dat<- filter(metadat, HU_8_NAME %in% c("Mad-Redwood", "Mattole", "Lower Eel", "Russion", "Tomales-Drake Bays", "Trinity", "South Fork Trinity", "Smith", "Lower Klamath", "Gualala-Salmon"))

# By EcoRgions

summary(as.factor(metadat$EcoRegion))
dat <- filter(metadat, grepl("^North Coast|North Coast$", EcoRegion))
summary(as.factor(metadat$HU_8_NAME))

# Sierras/BASIN
# dat <- filter(metadat, SPP_ID=="RABO" | SPP_pc1=="RABO") %>%
#   filter(EcoRegion=="Sierra Nevada" | EcoRegion=="Sierra/Basin Range")

# be aware outliers: 
## BEAR-MISS Canyon RAP-278
## SFY-SpringCk RAP 122


# 4c. FILTER BY SITE NAME -------------------------------------------------

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



# 04. FILTER METADAT TO SAMPLES OF INTEREST ------------------------

# filter out hybrids
#datHYB <- metadat %>% filter(grepl("RAP1745|RAP1587", SampleID))

# metadat %>% filter(!SPP_ID2=="RANA") %>% 
#   group_by(watershed, SPP_ID2) %>% 
#   tally

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

# # double check:
# dat %>% group_by(watershed, SPP_ID2) %>% tally

# 05. JOIN WITH RAW BAMLIST -----------------------------------------

# check col names for join...should be BAMFILE name with plate/well ID
dfout <- inner_join(dat, bams, by=c("Seq"="X1")) %>% arrange(Seq)

# dfout %>% group_by(watershed, SPP_ID2) %>% tally
# dfout %>% group_by(Locality, SPP_ID2) %>% tally

#dfout %>% filter(is.na(watershed)) %>% tally
# dfout %>% filter(is.na(Locality)) %>% tally

# 06. WRITE TO BAMLIST SINGLE ----------------------------------------

# Write to bamlist for angsd call (for IBS, no bamNo)
write_delim(as.data.frame(paste0("/home/rapeek/projects/rasi_hybrid/alignments/", dfout2$Seq, ".sortflt.mrg.bam")), path = paste0("data_output/bamlists/",site,"_",bamNo,"k_thresh.bamlist"), col_names = F)

# 07. PUT BAMLIST ON CLUSTER -----------------------------------------

# farmer (sftp)
# cd projects/....
paste0("put ",site,"_*",bamNo,"k*.bamlist") # (this goes from local to cluster)

# 08. ANGSD PCA (IBS) ------------------------------------------------

# create command:
lsite<- tolower(site)

# NEW IBS METHOD
paste0("sbatch -p high -t 48:00:00 --mail-type ALL --mail-user rapeek@ucdavis.edu scripts/03a_pca_new.sh ",site,"_",bamNo,"k_thresh.bamlist", " ", lsite, "_",bamNo,"k", " bamlists")




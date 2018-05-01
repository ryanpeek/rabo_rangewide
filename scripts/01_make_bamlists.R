## Create Bamlists for Analysis
## R. Peek 2018

# Updated Mon Apr 30 22:53:29 2018 ------------------------------

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
site <- "fea_rabo"

# 02. GET METADATA --------------------------------------------

# Read in metadata file here, with Seq column (i.e., the ID)
metadat<- read_csv("data_output/rapture06_metadata_revised.csv") %>% arrange(Seq)

# drop unneeded cols:
metadat <- select(metadat, -c(Regulated:HUC_12))

# filter out WESTUS Samples and all non-Sierra Samples:
# metadat <- filter(metadat, EcoRegion=="Sierra Nevada" | EcoRegion=="Sierra/Basin Range") %>%
#   filter(!HU_6_NAME=="San Joaquin")

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

# add watershed based on river/site name
metadat <- metadat %>%
  mutate(watershed= case_when(
    River == "BEAR" ~ "Bear",
    grepl("FEA|NFF|MFF", River)  ~ "Feather",
    grepl("NFY|MFY|SFY|DEER|FORD", River)  ~ "Yuba",
    grepl("NFA|MFA|RUB|NFMFA|SFA", River) ~ "American")
  )

metadat %>% group_by(watershed, SPP_ID2) %>% tally

#save(metadat, file = "data_output/rapture06_metadata_sierra.rda")

# 03. READ FULL BAMLIST --------------------------------------------------

# this can be a full bamlist (all bams) or subsampled bamlist

bams <- read_tsv(paste0("data/bamlists/bamlist_flt_mrg_",bamNo,"k"), col_names = F)

# remove the *000 component for join, requires fixing scipen for digits
subsamp<-bamNo*1000

# Fix the name in bamlist (based on subsample output)
bams$X1<-gsub(pattern = paste0(".sortflt.mrg_",subsamp,".bam"), replacement = "", bams$X1)

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




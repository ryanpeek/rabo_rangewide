## Clean Metadata File For Analysis of Rangewide Data
## R. Peek 2018

# Fri May  4 21:30:29 2018 ------------------------------

## Remove the Multispp panel from WSU
## Re-label missing EcoRegions/Localities
## Add "RIVER" Column
## Filter for ONLY RABO (or RANA that are RABO from WSU)

# 01. LOAD LIBRARIES & Global Settings ---------------------------

suppressPackageStartupMessages(library(tidyverse))
options(tibble.print_min = Inf) # print all rows 

# Global Setting: set digits to 12 to avoid paste function/joining issues
options(scipen = 12)
set.seed(111) # for repeatable random sampling

# LOAD METADATA --------------------------------------------

# Read in metadata file here, with Seq column (i.e., the ID)
metadat<- read_csv("data/rapture06_metadata_revised.csv") %>% arrange(Seq) %>%
  # drop unneeded cols:
  select(-c(Regulated:HUC_12))

# view a table interactively:
# DT::datatable(metadat)
table(metadat$EcoRegion)
table(metadat$SPP_ID)

# TIDY & CLEAN ------------------------------------------------------------

# FILTER OUT WESTUS (MULTISPP)  -------------------------------------------

# metadat_multispp <- filter(metadat, Locality=="WESTUS")
# write_csv(metadat_multispp, path = "data/rapture_metadata_multispp.csv")
metadat <- filter(metadat, !Locality=="WESTUS" , !SPP_ID=="ANBO" , !SPP_ID=="LICA")
table(metadat$SPP_ID)

# FIX SPP -----------------------------------------------------------------

# assign unknown spp
metadat <- metadat %>% 
  mutate(SPP_ID2=case_when(
    grepl("RAP1563|RAP1620|RAP1622|RAP1634|RAP1640|RAP1743|RAP1558|RAP1552|RAP1631|RAP1722|RAP1667", SampleID) ~ "RABO",
    grepl("RAP1745|RAP1587", SampleID) ~ "Hybrid",
    SPP_ID=="RANA" ~ "RASI",
    TRUE ~ SPP_ID))

table(metadat$SPP_ID2)

# FILTER TO RABO ONLY -----------------------------------------------------

metadat <- metadat %>% 
  filter(SPP_ID2=="RABO")

# FILTER OUT SAMPLES WITH NO LOC DATA -------------------------------------

# predominantly from EAST BRANCH FEATHER & PIT River
metadat %>% filter(is.na(lon)) %>% group_by(Locality) %>% tally

metadat <- filter(metadat, !is.na(lon))

table(metadat$HU_8_NAME)
table(metadat$EcoRegion)

# FIX HUC8's for PIT
# metadat <- metadat %>% 
#   mutate(HUC_8 = ifelse(Locality=="PIT", 18020003, HUC_8))
# summary(as.factor(metadat$HUC_8))

# ADD RIVER COLUMN (SEPARATE LOCALITY) ------------------------------------

metadat <- metadat %>%
  separate(Locality, c("River", "Site"), "-", remove = FALSE) #%>% 

metadat %>% group_by(River) %>% tally %>% print() 
# 42 different rivers 

metadat %>% group_by(HU_8_NAME) %>% tally %>% print()
# 36 different HUC 8s

metadat %>% group_by(HUC_8) %>% tally %>% print()
# 36 different HUC 8s

# SUMMARIZE ---------------------------------------------------------------

# how many samples per River
metadat %>% group_by(River) %>% distinct() %>% tally %>% print()

# how many samples per EcoRegion
metadat %>% group_by(EcoRegion) %>% distinct() %>% tally()


# SELECT COLS -------------------------------------------------------------

metadat_out <- metadat %>% select(-contains("barcode"), -RAD_ID, -PlateID, -Well, -SPP_ID, -SPP_pc1) %>% rename(SPP_ID=SPP_ID2)

# NEED TO JOIN WITH TISSUE TYPE AND DNA QUANT -----------------------------
quant_dna <- read_rds("data/merged_RANA_dna_metadata.rds") %>% select(SampleID, DNA.type)

# fix 3 values with quotation spaces:
quant_dna$SampleID[quant_dna$SampleID =="\"MH024\n\""] <- "MH024"
quant_dna$SampleID[quant_dna$SampleID =="\"RR012\n\""] <- "RR012"
quant_dna$SampleID[quant_dna$SampleID =="\"SM043\n\""] <- "SM043"

# now add the DNA.type
metadat_out <- left_join(metadat_out, quant_dna, by="SampleID")
quant_meta <- read_csv("data/rapture_quant_metadata.csv")
quant_wsu <- read_csv("data/rapture_quant_WSU_plates.csv")
# WRITE OUT ---------------------------------------------------------------

write_csv(metadat, path="data/rapture_metadata_rabo.csv")


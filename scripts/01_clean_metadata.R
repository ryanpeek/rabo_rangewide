## Clean Metadata File For Analysis of Rangewide Data
## R. Peek 2018

# Fri May  4 21:30:29 2018 ------------------------------

## Remove the Multispp panel from WSU
## Re-label missing EcoRegions/Localities
## Add "RIVER" Column
## Filter for ONLY RABO (or RANA that are RABO from WSU)

# 01. LOAD LIBRARIES & Global Settings ---------------------------

suppressPackageStartupMessages(library(tidyverse))

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
# summary(metadat$EcoRegion)

# TIDY & CLEAN ------------------------------------------------------------

## FILTER OUT WESTUS (MULTISPP) 

# metadat_multispp <- filter(metadat, Locality=="WESTUS")
# write_csv(metadat_multispp, path = "data/rapture_metadata_multispp.csv")
metadat <- filter(metadat, !Locality=="WESTUS")

## FIX PIT ECOREGION SAMPLES 

# fix the PIT samples so they are in Cascades EcoRegion
metadat <- metadat %>% 
  mutate(EcoRegion = case_when(
    EcoRegion=="UNKNOWN" ~ "Cascades",
    TRUE ~ EcoRegion),
    HU_8_NAME = case_when(
      Locality == "PIT" ~ "Lower Pit",
      TRUE ~ HU_8_NAME)
  )

## ADD RIVER as a column (SEPARATE LOCALITY)
metadat <- metadat %>%
  separate(Locality, c("River", "Site"), "-", remove = FALSE)

summary(as.factor(metadat$EcoRegion))
summary(as.factor(metadat$Locality))
summary(as.factor(metadat$River))

# FIX HUC8's for PIT
metadat <- metadat %>% 
  mutate(HUC_8 = ifelse(Locality=="PIT", 18020003, HUC_8))

summary(as.factor(metadat$HUC_8))

## FILTER OUT NON-RABO SAMPLES

metadat <- metadat %>% 
  filter(SPP_ID=="RABO" | SPP_pc1 == "RABO")

# SUMMARIZE ---------------------------------------------------------------

# how many samples per River
metadat %>% group_by(River) %>% distinct() %>% tally

# how many samples per EcoRegion
metadat %>% group_by(EcoRegion) %>% distinct() %>% tally

# WRITE OUT ---------------------------------------------------------------

write_csv(metadat, path="data/rapture_metadata_rabo.csv")


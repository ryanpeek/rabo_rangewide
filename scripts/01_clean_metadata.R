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
metadat04<- read_csv("data/rapture04_metadata_hucs_nhd_quant.csv") %>% 
  mutate(seqID = paste0(plate_barcode, "_", well_barcode)) %>% 
  select(seqID, PlateID:elev_m, Locality, -SPP_pc1, -PlateID, -Well) %>% 
  filter(!Locality=="WESTUS" , !SPP_ID=="ANBO" , !SPP_ID=="LICA")

# identify duplicates:
metadat04 %>% filter(duplicated(.[["LabID"]], incomparables = NA)) %>% tally
metadat04 %>% filter(duplicated(.[["LabID"]], incomparables = NA)) %>% View

# Read in metadata file here, with Seq column (i.e., the ID)
metadat<- read_csv("data/rapture06_metadata_revised.csv") %>% 
  mutate(seqID = paste0(plate_barcode, "_", well_barcode)) %>%
  # drop unneeded cols:
  select(-c(Seq, plate_barcode:well_barcode, SPP_pc1, Regulated:HUC_12)) %>% 
  select(seqID, everything()) %>% 
  # filter out WESTUS (multispp)
  filter(!Locality=="WESTUS" , !SPP_ID=="ANBO" , !SPP_ID=="LICA")

# to make specific WESTUS/diff spp:
# metadat_multispp <- filter(metadat, Locality=="WESTUS")
# write_csv(metadat_multispp, path = "data/rapture_metadata_multispp.csv")

# view a table interactively:
# DT::datatable(metadat)
table(metadat$EcoRegion)
table(metadat$SPP_ID)
table(metadat04$SPP_ID)

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

# MERGE WITH DNA CONC DATA ------------------------------------------------
#trim out cols we don't need:
metadat04 <- metadat04 %>% select(-c(SampleID:SPP_ID, YYYY, lat:elev_m, Locality))

# first see diff records:
anti_join(metadat, metadat04, by="seqID") %>% View

# now add the DNA_conc data
df1 <- left_join(metadat, metadat04, by="seqID")

table(df1$SPP_ID2)

# FILTER OUT SAMPLES WITH NO LOC DATA -------------------------------------

# predominantly from EAST BRANCH FEATHER & PIT River
df1 %>% filter(is.na(lon)) %>% group_by(Locality) %>% tally

datfilt <- filter(df1, !is.na(lon))

table(datfilt$HU_8_NAME)
table(datfilt$EcoRegion)

# check for duplicates:
datfilt %>% filter(duplicated(.[["LabID"]], incomparables = NA)) %>% tally # only Bean21 is RABO
#datfilt %>% filter(duplicated(.[["LabID"]], incomparables = NA)) %>% View

# ADD RIVER COLUMN (SEPARATE LOCALITY) ------------------------------------

datfilt <- datfilt %>%
  separate(Locality, c("River", "Site"), "-", remove = FALSE) #%>% 

datfilt %>% group_by(River) %>% tally %>% print() 
# 42 different rivers 

datfilt %>% group_by(HU_8_NAME) %>% tally %>% print()
# 36 different HUC 8s

datfilt %>% group_by(HUC_8) %>% tally %>% print()
# 36 different HUC 8s

# SUMMARIZE ---------------------------------------------------------------

# how many samples per River (42)
datfilt %>% group_by(River) %>% distinct() %>% tally %>% print()

# how many samples per EcoRegion (9)
datfilt %>% group_by(EcoRegion) %>% distinct() %>% tally()


# WRITE OUT RAW RABO DATA -------------------------------------------------

write_csv(datfilt, path = "data_output/rapture_metadata_rabo_raw.csv")
save(datfilt, file = "data_output/rapture_metadata_rabo_raw.rda")

# SELECT COLS -------------------------------------------------------------

datfilt_out <- datfilt %>% select(-RAD_ID, -PlateID, -Well, -SPP_ID) %>% rename(SPP_ID=SPP_ID2)

# NEED TO JOIN WITH TISSUE TYPE AND DNA QUANT -----------------------------

quant_dna <- read_rds("data/merged_RANA_dna_metadata.rds") %>% select(SampleID, DNA.type) %>% 
  rename(DNA_type=DNA.type)

# fix 3 values with quotation spaces:
quant_dna$SampleID[quant_dna$SampleID =="\"MH024\n\""] <- "MH024"
quant_dna$SampleID[quant_dna$SampleID =="\"RR012\n\""] <- "RR012"
quant_dna$SampleID[quant_dna$SampleID =="\"SM043\n\""] <- "SM043"

# now add the DNA.type and combine rm NAs
datfilt_out2 <- left_join(datfilt_out, quant_dna, by="SampleID") #%>% 
datfilt_out2 <- datfilt_out2 %>% 
  mutate(SampleTYPE=case_when(
    is.na(DNA_type) ~ SampleType, 
    is.na(SampleType) ~ DNA_type,
    TRUE ~ DNA_type)
  )

table(datfilt_out2$SampleType)
table(datfilt_out2$SampleTYPE)
table(datfilt_out2$DNA_type)

datfilt_out2 <- datfilt_out2 %>% select(-SampleType, -DNA_type) %>% 
  rename(DNA_Type=SampleTYPE)


# RECODE SAMPLE TYPES -----------------------------------------------------

datfilt_out2 <- datfilt_out2 %>% 
  mutate(DNA_Category=case_when(
    grepl("Tail|Tissue|Toe|Web|Tadpole", DNA_Type) ~ "Tissue",
    grepl("Buccal", DNA_Type) ~ "Buccal",
    grepl("Egg", DNA_Type) ~ "Egg",
    grepl("Skin", DNA_Type) ~ "Skin")
  )


table(datfilt_out2$DNA_Type)
table(datfilt_out2$DNA_Category)

datfilt_out2[duplicated(datfilt_out2$SampleID, incomparables = NA),] %>% select(seqID, SampleID, LabID, River, Site, SPP_ID:ng_ul) %>% View

datfilt_out2[duplicated(datfilt_out2$LabID, incomparables = NA),] %>% select(seqID, SampleID, LabID, River, Site, SPP_ID:ng_ul) %>% arrange(LabID) %>% View


# ADD DNA QUANT DATA ------------------------------------------------------

# quant_meta <- read_csv("data/rapture_quant_metadata.csv") %>% 
#   select(SampleID, SampleType:total_extract_ul)
# 
# quant_wsu <- read_csv("data/rapture_quant_WSU_plates.csv") %>% 
#   select(LabID, SampleType, Elevation_m, conc_ng_ul)


# WRITE OUT ---------------------------------------------------------------

write_csv(datfilt_out2, path="data_output/rapture_metadata_rabo_quant.csv")
save(datfilt_out2, file="data_output/rapture_metadata_rabo_quant.rda")

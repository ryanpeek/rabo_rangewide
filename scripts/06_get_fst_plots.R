
# GET FST Data ------------------------------------------------------------

# first we copy over from either farmer or using scp:
# scp -P 2022 rapeek@agri.cse.ucdavis.edu:/home/rapeek/projects/rangewide/pop_gen/results_fst/all_global_fst.txt .

# Load Libraries ----------------------------------------------------------

library(tidyverse)
library(stringr)

# GET DATA ----------------------------------------------------------------

dat <- read_csv("data_output/fst/all_global_fst.txt", col_names = "alldat")

# skip 1 line of every 3 rows           
dat <- dat[-seq(from=1, to=10599, by=3),]

# now pull out first line and add as a new col
datFiles <- dat[seq(from=1, nrow(dat), by=2),]
datFst <- dat[seq(from=2, nrow(dat), by=2),]

fsts <- bind_cols(datFiles, datFst) %>% 
  rename(fst=alldat1, filenames=alldat)

rm(datFiles, datFst)


# FIX WEIRD TXT -----------------------------------------------------------

fsts$filenames <- gsub(x = fsts$filenames, pattern = "_25k.folded.fst.gz", replacement = "")

# SEPARATE COLS -----------------------------------------------------------

fsts2 <- fsts %>% 
  separate(col = filenames, into = c("siteA", "siteB"), sep = "\\.") %>% 
  separate(col = fst, into = c("fst_unweight", "fst_weight"), sep=" ") %>% 
  separate(col=fst_unweight, into=c("label", "no_obs_loci", "fst_unweighted"), ":", remove = T) %>% 
  select(-label)

# sub/str fix text
fsts2$no_obs_loci <- gsub(x=fsts2$no_obs_loci, pattern="]", "")  
fsts2$fst_unweighted <- gsub(x=fsts2$fst_unweighted, pattern="FST.Unweight\\[", "")  
fsts2$fst_weight <- gsub(x=fsts2$fst_weight, pattern="Fst.Weight:", "")  


fsts2 <- fsts2 %>% 
  mutate_at(c("no_obs_loci", "fst_unweighted", "fst_weight"), .funs = as.numeric)

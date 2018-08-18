# RAPTURE bamlists for calculating Fst and Thetas
# 2017-Aug

# This script can be used to generate bamlists and bamlist_clst files for different subsampled groups, for Fst runs

# 1. LOAD LIBRARIES ----------------------------------------------------------
suppressMessages({
  library(tidyverse);
  library(lubridate);
  library(magrittr)
})

options(scipen = 12) # to avoid issues with paste functions, & joining on subsample number

# 2. GET DATA ----------------------------------------------------------------

# set subsample number (so need corresponding subsampled bamlist locally, see Steps 1A-1C)
bamNo<-100

# set site name (will be appended into filename)
site<-"SFY"

# METADATA
metadat<- read_csv("data_output/rapture06_metadata_revised.csv") %>% arrange(Seq)

# need to change Seq column to match the new merged one:
metadat$Seq <- gsub(pattern = "SOMM165", replacement = "SOMM163", x = metadat$Seq)

# metadat %>% 
#   filter(SPP_ID=="RABO" | SPP_ID=="RASI") %>% 
#   group_by(Locality, SPP_ID) %>% tally() %>% arrange(Locality) #%>% View()

# see sites:
metadat_sites <- sf::st_read("data/shps/fig1/metadat_sites.shp")
mapview::mapview(metadat_sites)

# 3. SUBSAMPLED BAMLIST --------------------------------------------------

bams <- read_tsv(paste0("data_output/bamlists/bamlist_flt_mrg_",bamNo,"k"), col_names = F)

# remove the *000 component for join, requires fixing scipen for digits
subsamp<-bamNo*1000
bams$X1<-gsub(pattern = paste0(".sortflt.mrg_",subsamp,".bam"), replacement = "", bams$X1)

# 4a. NFA Sites ------------------------------------------------------------

# this grabs all sites that start with NFA
dat <- filter(metadat, grepl('^NFA', Locality))

# Splits them by the Locality column into a list of unique dataframes
dat <- split(x = dat, dat$Locality)
names(dat)<-gsub(pattern = "-", x = tolower(names(dat)), replacement = "_")

# split apart into indiv data frames in GlobalEnv
list2env(x = dat, .GlobalEnv)

# rm all the data.frames starting with site
#rm(list = ls(pattern = "nfa"))

# 4b. MFA Sites ----------------------------------------------------------

# this grabs all sites that start with...
dat <- filter(metadat, SPP_ID=="RABO", grepl('^MFA|^RUB|^NFMFA', Locality)) 
 
# Splits them by the Locality column into a list of unique dataframes
dat <- split(x = dat, dat$Locality)
names(dat)<-gsub(pattern = "-", x = tolower(names(dat)), replacement = "_")
names(dat)

# split apart into indiv data frames in GlobalEnv
list2env(x = dat, .GlobalEnv)

# RENAME here if necessary
mfa_usox <- mfa_us_r
rm(mfa_us_r)
rub_lcus <- rub_lc_us
rm(rub_lc_us)

# rm all the data.frames starting with site
#rm(list = ls(pattern = site))

# 4c. Mossy Sites---------------------------------------------------------

# this grabs all sites that start with MFA 
dat <- filter(metadat, SPP_ID=="RASI", grepl('^FORD', Locality)) 

# Splits them by the Locality column into a list of unique dataframes
dat <- split(x = dat, dat$Locality)
names(dat)<-gsub(pattern = "-", x = tolower(names(dat)), replacement = "_")
names(dat)

# split apart into indiv data frames in GlobalEnv
list2env(x = dat, .GlobalEnv)


# 4d. FEA Sites---------------------------------------------------------

# this grabs all sites that start with... 
dat <- filter(metadat, SPP_ID=="RABO", grepl('^FEA|^NFF', Locality)) 

# Splits them by the Locality column into a list of unique dataframes
dat <- split(x = dat, dat$Locality)
names(dat)<-gsub(pattern = "-", x = tolower(names(dat)), replacement = "_")
names(dat)<-gsub(pattern = "fea_", x = tolower(names(dat)), replacement = "fea_Rb_")
names(dat)

# split apart into indiv data frames in GlobalEnv
list2env(x = dat, .GlobalEnv)

# rename/rm if need be:
rm(fea_Rb_spanish_wapaunsie)
rm(fea_Rb_spanish_rockck)
rm(fea_Rb_ebnff)
#rm(fea_Rs_rocklkbucksck)

fea_Rb_nff_poe <- nff_poe
rm(nff_poe)

# split out the Poe reach into two groups:
fea_Rb_nff_poe_pulga <- fea_Rb_nff_poe %>% filter(SampleID %in% c("RAP1580","RAP1628"))
fea_Rb_nff_poe_ds <- fea_Rb_nff_poe %>% filter(SampleID %in% c("RAP1568","RAP1604", "RAP1616"))

# rm all the data.frames starting with site
rm(list = ls(pattern = "fea"))
  
# 4e. BEAR Sites---------------------------------------------------------

# this grabs all sites that start with... 
dat <- filter(metadat, SPP_ID=="RABO", grepl('^BEAR', Locality)) 

# Splits them by the Locality column into a list of unique dataframes
dat <- split(x = dat, dat$Locality)
names(dat)<-gsub(pattern = "-", x = tolower(names(dat)), replacement = "_")
names(dat)

# split apart into indiv data frames in GlobalEnv
list2env(x = dat, .GlobalEnv)

# split stho into 
#sth2: 1:5, 9, 
#sthc: 6:8, 10:18


# sites are: sthc, sth2 (upper), stha (steep hollow @ hawkins)
# rename/rm if need be:
bear_sthc2 <- bear_sth2
rm(bear_sth2)
bear_sthc_haw <- bear_stha
rm(bear_stha)

# rm all the data.frames starting with site
rm(list = ls(pattern = "bear"))

# 4f. YUB Sites ---------------------------------------------------------

# this grabs all sites that start with... 
#dat <- filter(metadat, SPP_ID=="RABO", grepl('^SFY|^MFY|^NFY', Locality)) 
dat <- filter(metadat, SPP_ID=="RABO", grepl('^SFY', Locality)) 

# Splits them by the Locality column into a list of unique dataframes
dat <- split(x = dat, dat$Locality)
names(dat)<-gsub(pattern = "-", x = tolower(names(dat)), replacement = "_")
names(dat)

# split apart into indiv data frames in GlobalEnv
list2env(x = dat, .GlobalEnv)

# rename/rm if fewer than 3 samples:
rm(sfy_humbug); rm(sfy_thimc); # 100k/75k/50k

# combine nfy_slate with nfy_slate_cgrav
nfy_slate_cgrav <- bind_rows(nfy_slate, nfy_slate_cgrav)
rm(nfy_slate) # 100k/75k/50k  (So KEEP nfy_slate_cgrav and nfy_slate_onion)
#rm(mfy_remmington); # 100k/75k/50k
#rm(mfy_grizzly_ck); # only 2 obs

#rm(sfy_misc);
#rm(nfy_slate); rm(mfy_remmington); rm(sfy_scotchman);
#rm(nfy_slate_cgrav); rm(nfy_slate_onion)

# rm all the data.frames starting with site
#rm(list = ls(pattern = "sfy|nfy|mfy"))

# 4g. SIERRA Sites -------------------------------------------------------

# this grabs all sites that start with... 
dat <- filter(metadat, SPP_ID=="RABO", grepl('^STAN|^TUO|^SJOA|^SFA|^CALAV', Locality)) 

# Splits them by the Locality column into a list of unique dataframes
dat <- split(x = dat, dat$Locality)
names(dat)<-gsub(pattern = "-", x = tolower(names(dat)), replacement = "_")
names(dat)<-gsub(pattern = "^", x = tolower(names(dat)), replacement = "sierra_")
names(dat)

# split apart into indiv data frames in GlobalEnv
list2env(x = dat, .GlobalEnv)

# rename/rm if need be:
rm(sierra_stan_roseck); rm(sierra_stan_coyoteck); rm(sierra_sjoa_joseck);
rm(sierra_calav_esperck); rm(sierra_calav_jesusmariack)

# rm all the data.frames starting with site
rm(list = ls(pattern = "sierra"))


# 4h. TRIN Sites ---------------------------------------------------------

# this grabs all sites that start with... 
dat <- filter(metadat, SPP_ID=="RABO", grepl('^TRIN|KLAM', Locality)) 

# Splits them by the Locality column into a list of unique dataframes
dat <- split(x = dat, dat$Locality)
names(dat)<-gsub(pattern = "-", x = tolower(names(dat)), replacement = "_")
names(dat)<-gsub(pattern = " ", x = tolower(names(dat)), replacement = "_")
names(dat)

# split apart into indiv data frames in GlobalEnv
list2env(x = dat, .GlobalEnv)

# rename/rm if need be:
rm(trin_connerck); rm(trin_nftrinity); rm(trin_willowck); rm(trin_sftrinity)

# rm all the data.frames starting with site
rm(list = ls(pattern = "trin"))

# 4i. EEL Sites ---------------------------------------------------------

# this grabs all sites that start with... 
dat <- filter(metadat, SPP_ID=="RABO", grepl('^EEL|SFEEL', Locality)) 

# Splits them by the Locality column into a list of unique dataframes
dat <- split(x = dat, dat$Locality)
names(dat)<-gsub(pattern = "-", x = tolower(names(dat)), replacement = "_")
names(dat)<-gsub(pattern = " ", x = tolower(names(dat)), replacement = "_")
names(dat)

# split apart into indiv data frames in GlobalEnv
list2env(x = dat, .GlobalEnv)

# rename/rm if need be:
rm(sfeel_fish)

# rm all the data.frames starting with site
#rm(list = ls(pattern = "sfeel"))

# 4j. VANDUZ Sites ---------------------------------------------------------

# this grabs all sites that start with... 
dat <- filter(metadat, SPP_ID=="RABO", grepl('^VANDZ', Locality)) 

# Splits them by the Locality column into a list of unique dataframes
dat <- split(x = dat, dat$Locality)
names(dat)<-gsub(pattern = "-", x = tolower(names(dat)), replacement = "_")
names(dat)<-gsub(pattern = " ", x = tolower(names(dat)), replacement = "_")
names(dat)

# split apart into indiv data frames in GlobalEnv
list2env(x = dat, .GlobalEnv)

# rename/rm if need be:
rm(vandz_mill)

# rm all the data.frames starting with site
#rm(list = ls(pattern = "vandz"))

# 5. JOIN and MAKE BAMLIST-CLST FILES -------------------------------------

# need to make a list of site names:
sites<-ls(pattern = "nfy") 

# SIERRA: tuo / stan / calav / sfa / fea / vanduz / eel|sfeel / mfa|rub|nfmfa / trin|klam / sfy|mfy|nfy
# sjoa/ # no samples for sjoa worked for subsampled 50/75/100k

# COASTAL: paj / sal / dry / sancarp / soquel / ala / lagun / put / russ / mat / mad / smith / sumpqua

# function to join and write bamlist and clst
make_bamclst <- function(sitename) {
    dfout <- inner_join(get(sitename), bams, by=c("Seq"="X1")) %>% arrange(Seq)
    
    # Write to bamlist for angsd call:
    write_delim(as.data.frame(paste0(dfout$Seq, ".sortflt.mrg_",subsamp,".bam")),
                path = paste0("data_output/bamlists/bamlist_mrg_",
                              sitename,"_",bamNo,"k"), col_names = F)
    clst_out<-dfout %>%
      dplyr::select(HUC_10, SampleID, Locality) %>%
      dplyr::rename(FID=HUC_10, IID=SampleID, CLUSTER=Locality)
    head(clst_out)
    write_delim(clst_out, 
                path=paste0("data_output/bamlists/bamlist_mrg_",
                            sitename,"_",bamNo,"k_clst"))
}

map(sites, make_bamclst) # my god it works!!


# COARSE SCALE
# clst_out<-dfout %>%
#   dplyr::select(HU_8_NAME, SampleID, HUC_6) %>%
#   dplyr::rename(FID=HU_8_NAME, IID=SampleID, CLUSTER=HUC_6)
# length(unique(clst_out$CLUSTER))
# length(unique(clst_out$FID))

# FINE SCALE
# clst_out<-dfout %>%
#   dplyr::select(HUC_10, SampleID, Locality) %>%
#   dplyr::rename(FID=HUC_10, IID=SampleID, CLUSTER=Locality)
# length(unique(clst_out$CLUSTER))
# length(unique(clst_out$FID))
# clst_out %>% filter(is.na(FID)) %>% tally
# clst_out %>% filter(is.na(CLUSTER)) %>% tally


# 6. TERMINAL SFTP --------------------------------------------------------

# cd Documents/github/rabo_regulation/data_output/bamlists
# farmer #(sftp)
# cd projects/rana_rapture/MERGED/bamlists/subpops
# 'put bamlist_RABO_nfa_100k*' # (this goes from local to cluster)
sitename <- tolower(site)
paste0("put bamlist_mrg_",sitename,"_*",bamNo,"k*")

# 7. BASH: SAF.sh ---------------------------------------------------------

# calc sfs/SAF for FST
# create thetalists: 
# ls bamlist_mrg_nfa*k | grep "k$" | sed 's/bamlist_mrg_//g' > thetalist_yuba

# double check script for POPS or subpops (4 lines to change)

# sbatch -p high -t 12:00:00 06_theta_sfs.sh thetalist_yuba

# 8. BASH: Calc Pairwise FST ---------------------------------------------

# sbatch -t 24:00:00 -p high 07a_get_fst.sh pop1 pop2
# can check datestamp & sort by time: ls -lt FILE*

sites<-ls(pattern = "nfy") 

# get all unique pairs:
sitepairs <- combn(x = sites, m = 2)

bamNo<-100
for(i in 1:ncol(sitepairs)){
  cat(paste0("sbatch -p high -t 12:00:00 07_get_fst.sh ",sitepairs[,i][1],"_",bamNo,"k", " ",sitepairs[,i][2],"_",bamNo,"k"))
  cat("\n")
}

# sbatch -p high -t 12:00:00 07a_get_fst.sh mfy_oregck_50k mfy_us_oh_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh nfy_50k nfy_slate_cgrav_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh nfy_50k nfy_slate_onion_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh nfy_slate_cgrav_50k nfy_slate_onion_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh sfy_fallck_50k sfy_loga_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh sfy_fallck_50k sfy_mcki_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh sfy_fallck_50k sfy_misc_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh sfy_fallck_50k sfy_rockck_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh sfy_fallck_50k sfy_scotchman_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh sfy_fallck_50k sfy_shadyck_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh sfy_fallck_50k sfy_springck_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh sfy_loga_50k sfy_mcki_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh sfy_loga_50k sfy_misc_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh sfy_loga_50k sfy_rockck_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh sfy_loga_50k sfy_scotchman_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh sfy_loga_50k sfy_shadyck_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh sfy_loga_50k sfy_springck_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh sfy_mcki_50k sfy_misc_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh sfy_mcki_50k sfy_rockck_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh sfy_mcki_50k sfy_scotchman_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh sfy_mcki_50k sfy_shadyck_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh sfy_mcki_50k sfy_springck_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh sfy_misc_50k sfy_rockck_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh sfy_misc_50k sfy_scotchman_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh sfy_misc_50k sfy_shadyck_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh sfy_misc_50k sfy_springck_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh sfy_rockck_50k sfy_scotchman_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh sfy_rockck_50k sfy_shadyck_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh sfy_rockck_50k sfy_springck_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh sfy_scotchman_50k sfy_shadyck_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh sfy_scotchman_50k sfy_springck_50k
# sbatch -p high -t 12:00:00 07a_get_fst.sh sfy_shadyck_50k sfy_springck_50k


# 9. THETAS ----------------------------------------------------------------

# get list of all pops for thetas:
#ls bamlist_mrg*k | grep "k$" | sed 's/bamlist_mrg_//g' > thetalist_yuba


# 10. Get Files from SFTP -------------------------------------------------

# use sftp and copy over all the *finalfstout files
# can sort by most recent and list only those from specific date (here the 22nd)
# ls -lt *finalfstout | awk '$7=="22"{print $9}'

# then cat the files together for use below:
# cat *finalfstout > FST_all
# ls *finalfstout > FST_names
# paste FST_all FST_names > ALL_FST


# GET FILES FROM RSYNC ----------------------------------------------------

# rsync -avh -e "ssh -p 2022" rapeek@farm.cse.ucdavis.edu:~rapeek/projects/rangewide/pop_gen/results_fst/all_global_fst.txt


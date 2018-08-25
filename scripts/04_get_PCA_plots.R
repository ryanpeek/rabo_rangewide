

# GET covMat DATA ---------------------------------------------------------

# rsync -avh -e "ssh -p 2022" rapeek@farm.cse.ucdavis.edu:home/rapeek/projects/rangewide/pop_gen/all_rabo_filt10_1_100k_thresh.covMat

# scp -P 2022 rapeek@agri.cse.ucdavis.edu:/home/rapeek/projects/rangewide/pop_gen/all_rabo_filt10_1_100k_thresh.covMat .

# sftp: farmer
# cd projects/rangewide/pop_gen/results_pca
# get *covMat

# libraries ---------------------------------------------------------------

library(here)
library(tidyverse)

# GET PCA -----------------------------------------------------------------

# source the function to plot
source(paste0(here(),"/scripts/functions/f_read_covar.R"))
source(paste0(here(),"/scripts/functions/f_read_covar_rangewide.R"))

# load the metadata
metadat <- read_rds(paste0(here(), "/data_output/rapture_metadata_rabo_quant.rds"))

# need to make a new field to match the bam names (this is lame but whatever)
metadat <- metadat %>% 
  separate(seqID, into = c("barcode", "wellcode"), drop=T) %>% 
  mutate(Seq = paste0("SOMM163_", barcode, "_RA_GG", wellcode, "TGCAGG"))

# add groups based on Shaffer and PCA splits:
metadat<- metadat %>% 
  mutate(admix_groups = case_when(
    grepl("STAN|TUO|SFA", River) ~ "East", # southern siera
    grepl("ANTV|BEAR|DEER|MFA|MFY|NFA|NFMFA|NFY|SFY|RUB", River) ~ "North-East", # northern sierra
    grepl("CHETCO|SFEEL|VANDZ|TRIN|MAT|KLAM|SSANTIAM|PUT|MAD|LAGUN|SUMPQUA|RUSS|SMITH|EEL", River) ~ "North-West", # North Coast
    grepl("NFF|FEA", River) ~ "Feather-North", # feather
    grepl("PAJ|ALA|DRY|SOQUEL", River) ~ "West", # Central Coast
    grepl("SANCARP", River) ~ "South-West") # South Coast
  )



# set site/reads for bamlist/covar filepaths:
reads <- "100k_thresh"
site <- "rabo_nofeath_filt10_1"
(covarpath<- paste0(here(), "/data_output/pca/", site, "_", reads, ".covMat"))
(bampath <- paste0(here(), "/data_output/bamlists/", site, "_", reads, ".bamlist"))

# run function
(read_covar_range(covarpath, bampath, metadat, pcs = c(1,3), colvar = "ecoreg", plotlyplot = T))

#ggsave(filename = paste0("figs/pca_", site, "_", reads, "_pc1-3.png"), width = 8, height = 5, 
#       units = "in", dpi = 300)
        

# THE INNER BITS ----------------------------------------------------------

covar <- read.table(covarpath, stringsAsFactors = F)

# get bamlist for spp
bams <- read.table(bampath, stringsAsFactors = F, header = F) # read in the bamlist
bams$V2 <- sub('\\..*$', '', basename(bams$V1)) # remove the path and file extension
annot <- left_join(bams, metadat, by=c("V2"="Seq")) %>% 
  distinct(V2, .keep_all = T) %>% # drop any duplicates
  select(-V1) # join with the metadata

# Eigenvalues
eig <- eigen(covar, symm=TRUE)
eig$val <- eig$val/sum(eig$val)
PC <- as.data.frame(eig$vectors)
colnames(PC) <- gsub("V", "PC", colnames(PC))

# these are all local variables that can be added...ultimately up to personal pref for plot
PC$Pop <- factor(annot$River)
PC$Locality <- factor(annot$Locality) 
PC$ID <- factor(annot$SampleID)
PC$HU_12_NAME <- factor(annot$HU_12_NAME)
PC$HU_8_NAME <- factor(annot$HU_8_NAME)
PC$HUC6 <- factor(annot$HUC_6)
PC$LabID <- factor(annot$LabID)
PC$ecoreg <- factor(annot$EcoRegion)

# set up plot
pcs <- c(1,2)
colvar <- "ecoreg"

# PC's:
pc1 <- pcs[1]
pc2 <- pcs[2]


# Title: (% explained)
title <- paste("PC",pc1," (",signif(eig$val[pc1], digits=3)*100,"%)"," / PC", pc2," (",signif(eig$val[pc2], digits=3)*100,"%)", sep="",collapse="")

plotpca <- ggplotly(p = 
                      ggplot(data=PC, 
                             aes_string(x=paste0("PC",pc1), 
                                        y=paste0("PC",pc2), 
                                        color=colvar, 
                                        #shape=shapevar, # can be Pop, HU_8_NAME, HUC6, ecoreg etc 
                                        text=quote((paste0("ID: ", ID, 
                                                           " <br> River: ",
                                                           Locality))))) +
                      geom_point(size=4, alpha=0.8) + 
                      theme_bw(base_size = 9) +
                      scale_color_viridis_d() + 
                      ggtitle(paste0(title)))

plotpca

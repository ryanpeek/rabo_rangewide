
# libraries ---------------------------------------------------------------

library(here)
library(tidyverse)
library(plotly)

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


# set site/reads for bamlist/covar filepaths:
reads <- "25k"
site <- "all_rabo_filt01" # "all_rabo_filt01", "all_rabo_n4", "all_rabo_25k"
covarpath<- paste0(here(), "/data_output/pca/", site, "_", reads, ".covMat")
bampath <- paste0(here(), "/data_output/bamlists/", site, "_", reads, "_thresh.bamlist")

# run function
#(read_covar_range(covarpath, bampath, metadat, c(2,3),plotlyplot = TRUE))


# THE INNER BITS ----------------------------------------------------------

pcs <- c(1,2)
covar <- read.table(covarpath, stringsAsFactors = F)

# get bamlist for spp
bams <- read.table(bampath, stringsAsFactors = F, header = F) # read in the bamlist
bams$V2 <- sub('\\..*$', '', basename(bams$V1)) # remove the path and file extension
annot <- left_join(bams, metadat, by=c("V2"="Seq")) %>% 
  #row.names(annot[duplicated(annot$V2),])
  slice(c(-251,-419)) %>% # drop duplicated rows
  #distinct(V2, .keep_all = T) %>% # drop any duplicates
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

# PC's:
pc1 <- pcs[1]
pc2 <- pcs[2]

# Title: (% explained)
title <- paste("PC",pc1," (",signif(eig$val[pc1], digits=3)*100,"%)"," / PC", pc2," (",signif(eig$val[pc2], digits=3)*100,"%)", sep="",collapse="")

# PLOT
(plotpca <- ggplotly(p = 
                       ggplot(data=PC, 
                              aes_string(x=paste0("PC",pc1), 
                                         y=paste0("PC",pc2), 
                                         color="HUC6", 
                                         #shape="HUC6", 
                                         # can be Pop, HU_8_NAME, HUC6, ecoreg etc 
                                         text=quote((paste0("ID: ", ID, 
                                                            " <br> River: ",
                                                            Locality))))) +
                       geom_point(size=4, alpha=0.8) + 
                       theme_bw(base_size = 9) +
                       scale_color_viridis_d() + 
                       ggtitle(paste0(title))))
  

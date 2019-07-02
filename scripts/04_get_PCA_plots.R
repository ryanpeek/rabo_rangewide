

# GET covMat DATA ---------------------------------------------------------

# rsync -avh -e "ssh -p 2022" rapeek@farm.cse.ucdavis.edu:home/rapeek/projects/rangewide/pop_gen/all_rabo_filt10_1_100k_thresh.covMat

# scp -P 2022 rapeek@agri.cse.ucdavis.edu:/home/rapeek/projects/rangewide/pop_gen/all_rabo_filt10_1_100k_thresh.covMat .

# sftp: farmer
# cd projects/rangewide/pop_gen/results_pca
# get *covMat

# libraries ---------------------------------------------------------------

library(here)
library(tidyverse)

#devtools::install_github('thomasp85/ggforce')
library(ggforce)

# GET DATA AND FORMAT -----------------------------------------------------

#source(paste0(here(),"/scripts/functions/f_read_covar.R"))
source(paste0(here(),"/scripts/functions/f_read_covar_rangewide.R"))

# load the metadata
metadat <- read_rds(paste0(here(), "/data_output/rapture_metadata_rabo_quant.rds"))

metadat$Locality<-gsub(pattern = "[[:space:]]", replacement = "-", x = metadat$Locality)

# need to make a new field to match the bam names (this is lame but whatever)
metadat <- metadat %>% 
  separate(seqID, into = c("barcode", "wellcode"), drop=T) %>% 
  mutate(Seq = paste0("SOMM163_", barcode, "_RA_GG", wellcode, "TGCAGG"))

# unique localities?
metadat %>% distinct(Locality) %>% tally
#metadat %>% distinct(Locality) %>% arrange(Locality) %>% View

# add groups based on Shaffer and PCA splits:
metadat<- metadat %>% 
  mutate(
    admix_orig = case_when(
      grepl("STAN|TUO|CALAV", River) ~ "East", # southern siera
      grepl("SFA|ANTV", River) ~ "Unknown",
      grepl("NFF|FEA|BEAR|DEER|MFA|MFY|NFA|NFMFA|NFY|SFY|RUB", River) ~ "North-East", # northern sierra
      grepl("CHETCO|SFEEL|COW|VANDZ|TRIN|MAT|KLAM|SSANTIAM|PUT|^MAD$|LAGUN|SUMPQUA|RUSS|SMITH|EEL", River) ~ "North-West", # North Coast
      #grepl("NFF|FEA", River) ~ "North-Feather", # feather
      grepl("PAJ|ALA|DRY|SOQUEL", River) ~ "West", # Central Coast
      grepl("SANCARP|SALIN", River) ~ "South-West"), # South Coast
    admix_groups = case_when(
      grepl("STAN|TUO|CALAV|SFA", River) ~ "East", # southern siera
      grepl("BEAR|DEER|MFA|MFY|NFA|NFMFA|NFY|SFY|RUB", River) ~ "North-East", # northern sierra
      grepl("CHETCO|SFEEL|COW|VANDZ|TRIN|MAT|KLAM|SSANTIAM|PUT|^MAD$|LAGUN|SUMPQUA|RUSS|SMITH|EEL", River) ~ "North-West", # North Coast
      grepl("NFF|FEA", River) ~ "North-Feather", # feather
      grepl("PAJ|ALA|DRY|SOQUEL", River) ~ "West", # Central Coast
      grepl("SANCARP|SALIN", River) ~ "South-West") # South Coast
  )



# SET BAMLIST AND COVMAT --------------------------------------------------

# set site/reads for bamlist/covar filepaths:
reads <- "100k_thresh"
site <-  "all_rabo_filt10_1" 
# all_rabo_filt : all rabo samples minus outliers
# all_rabo_filt_1 : all rabo minus outliers and all localties
# all_rabo_filt10_1 : all rabo minus outliers and all subsampled to 10 or less samples
# all_rabo_filt10 : all rabo minus outliers and samples > 2 per location

(covarpath<- paste0(here(), "/data_output/pca/", site, "_", reads, ".covMat"))
(bampath <- paste0(here(), "/data_output/bamlists/", site, "_", reads, ".bamlist"))

# run function
(read_covar_range(covarpath, bampath, metadat, pcs = c(1,2), colvar = "admix_orig", plotlyplot = T))

# PARSE COVAR FILE FOR PCA ------------------------------------------------

covar <- read.table(covarpath, stringsAsFactors = F)

# get bamlist for spp
bams <- read.table(bampath, stringsAsFactors = F, header = F) # read in the bamlist
bams$V2 <- sub('\\..*$', '', basename(bams$V1)) # remove the path and file extension
annot <- left_join(bams, metadat, by=c("V2"="Seq")) %>% select(-V1) # join with the metadata

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
PC$admix_groups<-factor(annot$admix_groups)
PC$admix_orig<-factor(annot$admix_orig)

# weird issue with a single sample from North-West RAP-111, (MAD, clusters wiht Sierras near SFY). Likely mis-labeled or ID'd.
PC <- PC %>% dplyr::filter(!ID=="RAP-111")

# SETUP PLOT --------------------------------------------------------------

# define colors
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#CFCFCF")

# define PCs
pcs <- c(1,2)
colvar <- "admix_orig" #"admix_groups"

# PC's:
pc1 <- pcs[1]
pc2 <- pcs[2]

# Title: (% explained)
title <- paste("PC",pc1," (",signif(eig$val[pc1], digits=3)*100,"%)"," / PC", pc2," (",signif(eig$val[pc2], digits=3)*100,"%)", sep="",collapse="")

# PLOT AND SAVE -----------------------------------------------------------

#plotly::ggplotly(
(ggpca <- ggplot(data=PC, aes_string(x=paste0("PC",pc1), y=paste0("PC",pc2),
                                     color=colvar, #shape=shapevar)
                                     text=quote((paste0("ID: ", ID, " <br> River: ", Locality))))) +
  
  geom_point(size=4, alpha=0.8) +
   # need to do this to get the lassen point to show up
   geom_point(data=PC[PC$ID=="RAP-327",], aes_string(x=paste0("PC",pc1), y=paste0("PC",pc2),
                                  color=colvar), size=4, alpha=0.8) +
  #ggforce::geom_mark_ellipsis(data=PC, aes_string(x=paste0("PC",pc1), y=paste0("PC",pc2), color=colvar, group=colvar))+
  theme_bw(base_family = "Roboto Condensed") +
  theme(legend.position="bottom") +
  # fix colors
  scale_color_manual("Admix Groups", 
                     values = c("East"=cbbPalette[1], # E
                                "North-East"=cbbPalette[2], # NE 
                                "North-West"=cbbPalette[3], # NW
                                #"North-Feather"=cbbPalette[4], #N Feather
                                "West"=cbbPalette[5],  # W
                                "South-West"=cbbPalette[6],# SW
                                "Unknown"=cbbPalette[9]), # SFAmerican
                     labels = c("S. Sierra (E)", # E
                                "N. Sierra (NE)", # NE 
                                "N. Coast (NW)", # NW
                                #"N. Sierra-Feather", #N Feather
                                "S. Coast (W)",  # W
                                "Unknown",
                                "C. Coast (SW)")) +
  guides(col=guide_legend(nrow = 2, direction = "horizontal")) +#byrow = TRUE)) +
  ggtitle(paste0(title))
)
#)


# SAVE OUT ----------------------------------------------------------------

ggpca_12 <-ggpca
#ggpca_34 <-ggpca
#ggpca_56 <-ggpca

plotly::ggplotly(ggpca_12)

ggpca_12
ggpca_34

ggsave(filename = paste0("figs/pca_", site, "_", reads, "_pc",pcs[1],"-",pcs[2], "_orig.png"), width = 7, height = 5, 
       units = "in", dpi = 300)
#ggsave(filename = paste0("figs/pca_", site, "_", reads, "_pc",pcs[1],"-",pcs[2], "_orig_ellipsis.png"), width = 7, height = 5, 
#       units = "in", dpi = 300)


# STACK PLOTS -------------------------------------------------------------

# https://cran.r-project.org/web/packages/cowplot/vignettes/shared_legends.html        

library(cowplot)

# arrange in one row with no legend
(pcastack <- plot_grid(ggpca_12 + theme(legend.position="none"), ggpca_34 + theme(legend.position = "none"), 
                       align='vh', hjust=-1, nrow=1, labels = "AUTO"))

#extract legend from a single plot
plotLeg <- get_legend(ggpca_12)

# now add to the plot
(pcastack_leg <- plot_grid(pcastack, plotLeg, nrow = 2, rel_heights = c(1, .2))) # fuck yes!!

# filt_10_1_100k try annotate? wow this is cool
(pca_final <- pcastack_leg + 
    annotate("text", x = 0.3, y=.85, label="Feather", col="gray30", size=3.5, fontface="italic") +
    annotate("text", x = 0.4, y=0.65, col="gray30", label="SF American",fontface="italic", size=3.5) +
    annotate("text", x = 0.85, y=0.45, label="Feather", col="gray30", size=3.5, fontface="italic") +
    annotate("text", x = 0.75, y=0.85, col="gray30", label="SF American",fontface="italic", size=3.5) +
    annotate("text", x = 0.71, y=0.35, col="gray30", label="Lassen",fontface="italic", size=3.5) +
    annotate("text", x = 0.37, y=0.78, col="gray30", label="Lassen",fontface="italic", size=3.5))


# filt_100k try annotate? wow this is cool
(pca_final <- pcastack_leg + 
    annotate("text", x = 0.27, y=.85, label="Feather", col="gray30", size=3.5, fontface="italic") +
    annotate("text", x = 0.23, y=0.6, col="gray30", label="SF American",fontface="italic", size=3.5) +
    annotate("text", x = 0.87, y=0.6, label="Feather", col="gray30", size=3.5, fontface="italic") +
    annotate("text", x = 0.74, y=0.7, col="gray30", label="SF American",fontface="italic", size=3.5))

save_plot(plot = pca_final, filename = paste0("figs/pca_12_34_", site, "_", reads, ".png"), base_width = 6, base_height = 4, base_aspect_ratio = 1.3, dpi=300)

#save_plot(plot = pcaquad, filename = paste0("figs/pca_12_34_", site, "_", reads, ".png"), base_width = 4.5, base_height = 4.5, dpi=300)

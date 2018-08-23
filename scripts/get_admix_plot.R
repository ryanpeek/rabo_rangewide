
# GET ADMIX AND PLOT ------------------------------------------------------


# Load Libraries ----------------------------------------------------------

library(viridis)
library(tidyverse)

# GET BAMLIST AND METADATA ------------------------------------------------

bamfile <- "all_rabo_100k"

# Get ID and pop info for each individual from bamlists
bams <- read.table(paste0("data_output/bamlists/",bamfile, "_thresh.bamlist"),stringsAsFactors = F, header = F)
bams$V2 <- sub('\\..*$', '', basename(bams$V1)) # remove the path and file extension

# get metadat
metadat <- read_rds(path = "data_output/rapture_metadata_rabo_quant.rds")

# join together
annot <- left_join(bams, metadat, by=c("V2"="Seq")) %>% select(-V1) # join with the metadata

# check for duplicates?
annot[duplicated(annot$V2),]
metadat[duplicated(metadat$Seq),]

# change ecoregion to factor
annot$EcoRegion <- as.factor(annot$EcoRegion)

# reorder factor:
annot$EcoRegion <- factor(annot$EcoRegion, levels=c("Cascades", "North Coast", "Klamath/North Coast", "Northern CA Coastal Foothills", "Sierra/Basin Range" ,"Sierra Nevada", "Central CA Coastal Foothills",
                                                    "Central Coast"))

levels(annot$EcoRegion)
#annot %>% arrange(EcoRegion) %>% select(SampleID, EcoRegion, Locality, HU_6_NAME, HU_8_NAME) %>% View
# recode levels based on Shaffer paper (E=Southern Sierra Nevada, NE=Northern Sierra Nevada, NW=North Coast, W=Central Coast, SW=South Coast, NEW::::Sierra/Basin Range)

#annot<- annot %>% 
#  mutate(admix_groups = case_when()


# GET QOPT ADMIX DATA -----------------------------------------------------

k <- 7

#admixfile <- paste0("data_output/admix/", qopt_file, "k", k,"_admix", admixk, ".qopt")
# Read inferred admixture proportions
q<-read.table(paste0("data_output/admix/", bamfile, "_k", k,"_admix.qopt"))

# join with annot data
ord <- order(annot$EcoRegion, q$V1)
#ord <- order(q$V1,q$V2, q$V4)

# PLOT Ks ----------------------------------------------------------------

# Plot them (ordered by population)
png(filename = paste0("figs/admix/",bamfile, "_k",k, "_admix.png"), width = 9, height = 6, units = "in", res = 300)

#setup
par(mar=c(10,4.5,2,1))

# test to see alignment of groups:
#tst <- barplot(t(q)[,ord],col=viridis(k),names=annot$EcoRegion[ord], space=0, las=2,ylab="Admixture proportions",cex.names=0.6, border=NA, main= paste0("Admixture k=", k, " (", bamfile,")"))

# now replot for saving
barplot(t(q)[,ord],col=viridis(k),names=annot$EcoRegion[ord], space=0, las=2,ylab="Admixture proportions",cex.names=0.6, border=NA, main= paste0("Admixture k=", k, " (", bamfile,")"), xaxt='n')

#unique(annot$EcoRegion[ord])
text(x = c(2, 40, 120, 185, 250, 440, 600, 630), par("usr")[3] - 0.2, labels = unique(annot$EcoRegion[ord]) , srt = 90, pos = 1, xpd = TRUE, cex = 0.8) #col = "blue")

# ecoregion breaks
abline(v = 7, col="black", lty=2)
abline(v = 71, col="black", lty=2)
abline(v = 171, col="black", lty=2)
abline(v = 204, col="black", lty=2)
abline(v = 296, col="black", lty=2)
abline(v = 582, col="black", lty=2)
abline(v = 582, col="black", lty=2)
abline(v = 624, col="black", lty=2)

dev.off()


# POPHELPER ---------------------------------------------------------------

library(pophelper)

#k <- 5
#admixk <- 5
#qopt_file <- "all_rabo_n4_25k_"
admixfiles <- list.files(path = "data_output/admix/", "*.qopt", full.names = T)
#admixfile <- paste0("data_output/admix/", qopt_file, "k", k,"_admix", admixk, ".qopt")

# read in multiple files
poplist <- readQ(admixfiles)

# tabulate
tabulateQ(qlist=poplist, sorttable = F)

# summarize
summariseQ(tabulateQ(qlist=poplist))

# GET LABELS:
admix_labels <- annot[,c("EcoRegion", "HU_6_NAME")]

# replace NA's
admix_labels <- admix_labels %>% replace(., is.na(.), "unknown")

# make factors
# admix_labels <- admix_labels %>% 
#   mutate_at(.vars = c("HU_6_NAME", "River", "EcoRegion"), as.factor)
# summary(admix_labels)

par(mar=c(6,1,1,1))
p1 <- plotQ(poplist[9], 
            #imgoutput="join", 
            returnplot = F, exportplot = T, quiet=T,
            width = 6, height=5, units = "in", 
            font = "Roboto Condensed", 
            #sharedindlab = T, 
            sortind = 'all', 
            splab = paste0("K=",sapply(poplist[c(8)],ncol)),splabsize = 5, 
            #subsetgrp=c("North Coast","Sierra Nevada", "Central Coast"), 
            selgrp="EcoRegion",  
            showlegend=T, showtitle=F,showsubtitle=F,titlelab="Structure",
            grplabangle = 90, grplabsize = 2, grplab=admix_labels, ordergrp=T, 
            grplabpos = 0.4, grplabheight = 1, grplabjust = 0.2, grplabspacer = .1,
            linepos = 0.9,
            #showindlab=T, indlabsize=2, indlabheight=0.1, indlabspacer=-1,indlabangle=70,
            barbordercolour="white", barbordersize=0, 
            outputfilename="figs/admix/all_n4_25k_k3_k5", imgtype="png")

dev.off()
print(p1$plot[[1]])

# ADMIXTURE 25k thresh ------------------------------------------------------

colnames(q) <- c("k", "proportion")
q2 <- bind_cols(annot, q)

# join with annot data
ord <- order(q2$proportion)

#write_csv(q2, path = "data_output/admixture_25k_thresh_meta_qopt.csv")

# Plot them (ordered by population)
#par(mar=c(7,3,1,1))
barplot(t(q2[,c("k","proportion")])[,ord],col=viridis(2),names=annot$SPP_ID2[ord],las=2,ylab="Admixture proportions",cex.names=0.6, border=NA)

# ggplot version
#ggplot() + geom_col(data=q2[ord,], aes(x=row.names(q2), y=proportion, 
#                                 fill=SPP_ID2), position = "stack"#) + 
#  scale_fill_viridis_d() + theme_classic()



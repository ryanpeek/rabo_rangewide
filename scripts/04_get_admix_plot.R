
# GET ADMIX AND PLOT ------------------------------------------------------


# Load Libraries ----------------------------------------------------------

library(viridis)
library(tidyverse)

# ADMIXTURE 25k ------------------------------------------------------

# Get ID and pop info for each individual from clst file
annot <- read.table("data_output/bamlists/bamlist_mrg_RASI_all_25k_clst",stringsAsFactors = F, header = TRUE)

metadat25k <- read_rds("data_output/bamlist_dat_RASI_all_25k_thresh.rds")

annot <- inner_join(metadat25k, annot, by=c("SampleID"="IID")) %>% arrange(Seq) %>% distinct(SampleID, .keep_all = T)

# Read inferred admixture proportions
q<-read.table("data_output/admix/admix_rasi_all_25k_out.qopt")

# join with annot data
ord <- order(annot$SPP_ID)

# Plot them (ordered by population)
#par(mar=c(7,3,1,1))
barplot(t(q)[,ord],col=viridis(2),names=annot$watershed[ord],las=2,ylab="Admixture proportions",cex.names=0.6, border=NA)

# ADMIXTURE 25k thresh ------------------------------------------------------

# Get ID and pop info for each individual from clst file
annot <- read.table("data_output/bamlists/bamlist_mrg_RASI_all_25k_clst",stringsAsFactors = F, header = TRUE)

metadat25k <- read_rds("data_output/bamlist_dat_RASI_all_25k_thresh.rds")

annot <- inner_join(metadat25k, annot, by=c("SampleID"="IID")) %>% arrange(Seq) %>% distinct(SampleID, .keep_all = T)

# refactor
annot$SampleID <- factor(annot$SampleID)
str(q2$SampleID)

# Read inferred admixture proportions
q<-read.table("data_output/admix/admix_rasi_all_25k_thresh_out.qopt")

colnames(q) <- c("k", "proportion")
q2 <- bind_cols(annot, q)

# join with annot data
ord <- order(q2$proportion)

write_csv(q2, path = "data_output/admixture_25k_thresh_meta_qopt.csv")

# Plot them (ordered by population)
#par(mar=c(7,3,1,1))
barplot(t(q2[,c("k","proportion")])[,ord],col=viridis(2),names=annot$SPP_ID2[ord],las=2,ylab="Admixture proportions",cex.names=0.6, border=NA)

# ggplot version
#ggplot() + geom_col(data=q2[ord,], aes(x=row.names(q2), y=proportion, 
#                                 fill=SPP_ID2), position = "stack"#) + 
#  scale_fill_viridis_d() + theme_classic()

# ADMIXTURE 50k ------------------------------------------------------

# Get ID and pop info for each individual from clst file
annot <- read.table("data_output/bamlists/bamlist_mrg_RASI_all_50k_clst",stringsAsFactors = F, header = TRUE)

# combine with metadata
#suppressMessages(metadat<- read_csv("data_output/rapture06_metadata_revised.csv"))

#annot <- inner_join(metadat, annot, by=c("SampleID"="IID")) %>% arrange(Seq) %>% distinct(SampleID, .keep_all = T)

# Read inferred admixture proportions
q<-read.table("data_output/admix/admix_rasi_all_50k_out.qopt")

# join with annot data
ord <- order(annot$FID)

# Plot them (ordered by population)
#par(mar=c(7,3,1,1))
barplot(t(q)[,ord],col=viridis(2),names=annot$CLUSTER[ord],las=2,ylab="Admixture proportions",cex.names=0.6, border=NA)




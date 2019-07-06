# basic PCA

# set up packages
if(!require(viridis)) { install.packages("viridis"); require(viridis, warn.conflicts = F)}
if(!require(ggplot2)) { install.packages("ggplot2"); require(ggplot2, warn.conflicts = F)}
if(!require(dplyr)) { install.packages("dplyr"); require(dplyr, warn.conflicts = F)}
if(!require(plotly)) { install.packages("plotly"); require(plotly, warn.conflicts = F)}
library(readr) # read in csv
library(janitor) # fix col names

# get covMat and bamlist
#covmat <- "data_output/pca/rothstein.covMat"
covmat <- "data_output/pca/rothstein_pca.cov.npy"
bamlist <- "data_output/bamlists/rothstein.bamlist"
barcodes <- "data_output/bamlists/96BestRADBarcodesCorrected.txt"

# read in metadata/platemap
metadata <- read_csv("data/rapture_pilot_arothstein_samples.csv")
metadata <- janitor::clean_names(metadata)
barcodeLS <- read_delim(barcodes,delim = "\t", col_names = c("well","barcode"))

# add spp info
metadata$col <- as.integer(substr(metadata$well, 2,3))
metadata$spp <- ifelse(metadata$col > 9, "panama", "rasi")

# read in data

# if using pcAngsd
library(RcppCNPy)
covar <- npyLoad(covmat) # Reads in estimated covariance matrix

# if using IBS
# covar <- read.table("data_output/pca/all_rabo_100k_thresh.covMat", stringsAsFactors = F)

# simple test
#e <- eigen(covar)
#plot(e2$vectors[,1:2],lwd=2,ylab="PC 2",xlab="PC 2",main="Principal components",col="red",pch=21)

# read in bamlist and trim to barcode
bams <- read.table(bamlist, stringsAsFactors = F, header = F) # read in the bamlist
bams$V2 <- sub('\\..*$', '', basename(bams$V1)) # remove the path and file extension
bams$barcode <- stringr::str_remove_all(string = bams$V2, pattern = "SOMM446_TGCGCT_GG")
bams$barcode <- stringr::str_remove_all(string = bams$barcode, pattern = "TGCAGG")

# now join with metadat by well
metadata <- left_join(metadata, barcodeLS, by="well")
annot <- left_join(bams, metadata, by=c("barcode")) %>% select(-V1) # join

# get pca vals
eig <- eigen(covar, symm=TRUE)
eig$val <- eig$val/sum(eig$val)
PC <- as.data.frame(eig$vectors)
colnames(PC) <- gsub("V", "PC", colnames(PC))

# these are all local variables that can be added...ultimately up to personal pref for plot
PC$ID <- factor(annot$sample_id)
PC$ng_ul <- factor(annot$qubit_ng_u_l)
PC$spp <- factor(annot$spp)
PC$type <- factor(annot$sample_type)

# set up plot
pcs <- c(1,2)

# PC's:
pc1 <- pcs[1]
pc2 <- pcs[2]

# Title: (% explained)
title <- paste("PC",pc1," (",signif(eig$val[pc1], digits=3)*100,"%)"," / PC", pc2," (",signif(eig$val[pc2], digits=3)*100,"%)", sep="",collapse="")

# static pca: type and ng_ul
(gg12 <- ggplot(data=PC, 
       aes_string(x=paste0("PC",pc1), y=paste0("PC",pc2),
                  color=quote(as.numeric(ng_ul)), shape=quote(type),
                  text=quote(paste0("ID=",ID)))) +
  geom_point(size=4, alpha=0.8) +
  #ggforce::geom_mark_circle(aes(group=watershed))+
  theme_bw() +
  #theme(legend.position="bottom") +
  scale_color_viridis_c("DNA Conc.\n (ng/ul)") + 
  #scale_color_viridis_d() + 
  ggtitle(paste0(title)))
ggsave(filename = "figs/rothstein_pca_1_2_dnaconc.pdf", width = 11, height = 8, units = "in")

ggplotly(gg1)

# static pca 2: spp and type
(gg2 <- ggplot(data=PC, 
       aes_string(x=paste0("PC",pc1), y=paste0("PC",pc2),
                  color=quote(spp), shape=quote(type), text=quote(paste0("ID=",ID)))) +
  geom_point(size=4, alpha=0.8) +
  #ggforce::geom_mark_circle(aes(group=watershed))+
  theme_bw() +
  #theme(legend.position="bottom") +
  #scale_color_viridis_c("DNA Conc.\n (ng/ul)") + 
  scale_color_viridis_d("Spp") + 
  ggtitle(paste0(title)))

ggsave(filename = "figs/rothstein_pca_1_2_spp.pdf", width = 11, height = 8, units = "in")

# plotly
ggplotly(gg2)


# 2 v 3
# set up plot
pcs <- c(1,3)

pc1 <- pcs[1]
pc2 <- pcs[2]

# Title: (% explained)
title <- paste("PC",pc1," (",signif(eig$val[pc1], digits=3)*100,"%)"," / PC", pc2," (",signif(eig$val[pc2], digits=3)*100,"%)", sep="",collapse="")


(gg23a <- ggplot(data=PC, 
                 aes_string(x=paste0("PC",pc1), y=paste0("PC",pc2),
                            color=quote(spp), shape=quote(type), text=quote(paste0("ID=",ID)))) +
    geom_point(size=4, alpha=0.8) +
    theme_bw() +
    #theme(legend.position="bottom") +
    scale_color_viridis_d("Spp") + 
    ggtitle(paste0(title)))

(gg23b <- ggplot(data=PC, 
                aes_string(x=paste0("PC",pc1), y=paste0("PC",pc2),
                           color=quote(as.numeric(ng_ul)), shape=quote(type),
                           text=quote(paste0("ID=", ID)))) +
    geom_point(size=4, alpha=0.8) +
    #ggforce::geom_mark_circle(aes(group=watershed))+
    theme_bw() +
    #theme(legend.position="bottom") +
    scale_color_viridis_c("DNA Conc.\n (ng/ul)") + 
    #scale_color_viridis_d() + 
    ggtitle(paste0(title)))


ggsave(filename = "figs/rothstein_pca_1_2_dnaconc.pdf", width = 11, height = 8, units = "in")

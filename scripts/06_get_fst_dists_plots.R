
# GET FST Data ------------------------------------------------------------

# first we copy over from either farmer or using scp:
# scp -P 2022 rapeek@agri.cse.ucdavis.edu:/home/rapeek/projects/rangewide/pop_gen/results_fst/all_global_fst.txt .

# Load Libraries ----------------------------------------------------------

library(tidyverse)
library(stringr)
library(here)
library(viridis)

# GET DATA ----------------------------------------------------------------

# make spaces to tab delimited
# :%s/ /\t/g
# :%s/_100k\.folded\.finalFSTout//g
# need pop 1, pop 2, fst for perlscript

dat <- read_tsv("data_output/fst/all_rabo_100k_folded_adj_fst.txt", col_names = c("fst_unad", "fst_adj", "filenames"))

# FIX WEIRD TXT -----------------------------------------------------------

dat$filenames <- gsub(x = dat$filenames, pattern = "_100k.folded.finalFSTout", replacement = "")

# SEPARATE COLS -----------------------------------------------------------

# separate and cleanup
fsts <- dat %>% 
  separate(col = filenames, into = c("siteA", "siteB"), sep = "\\.") %>% 
  mutate(sitepair = paste0(siteA, "_", siteB))

fsts %>% group_by(siteA) %>% tally
fsts %>% group_by(siteA) %>% tally

#fsts_only <- fsts %>% select(siteA, siteB, fst_adj)

# Get Metadata and Clean --------------------------------------------------

# load the metadata
metadat <- read_rds(paste0(here(), "/data_output/rapture_metadata_rabo_quant.rds"))

# fix trinity spaces
metadat$Locality<-tolower(gsub(pattern = "[[:space:]]", replacement = "-", x = metadat$Locality))

# fix deer-clearck/ deer-clec
metadat$Locality <- gsub(pattern="deer-clearck", replacement = "deer-clec", x=metadat$Locality)

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
    grepl("SANCARP|SALIN", River) ~ "South-West") # South Coast
  ) %>% 
  select(Seq, admix_groups, Locality, lat, lon: NHD_Tot_DA_sqkm, River, Site, EcoRegion)

# set order in way that you want
ords_admix_grps <- c("East", "North-East", "Feather-North", "North-West", "South-West", "West")

metadat$admix_groups <- factor(metadat$admix_groups, levels = ords_admix_grps)


# GET BAMILES -------------------------------------------------------------

bamfile <- "all_rabo_filt_100k"

# Get ID and pop info for each individual from bamlists
bams <- read.table(paste0("data_output/bamlists/",bamfile, "_thresh.bamlist"),stringsAsFactors = F, header = F)
bams$V2 <- sub('\\..*$', '', basename(bams$V1)) # remove the path and file extension

annot <- left_join(bams, metadat, by=c("V2"="Seq")) %>% select(-V1) # join with the metadata

annot %>% group_by(Locality) %>% tally

# get only the sites for pairings
annot <- annot %>% distinct(Locality, .keep_all = T)
annot$Locality <- tolower(annot$Locality)

# JOIN DATA ---------------------------------------------------------------

# now join with fst data:
annot_fst <- annot %>% select(admix_groups, Locality, lat, lon, HUC_6, EcoRegion)

fstsA <- left_join(fsts, annot_fst, by=c("siteA"="Locality")) %>% 
  rename_at(c("admix_groups","lat","lon", "HUC_6","EcoRegion"), funs( paste0(., "_A")))

fsts_out <- left_join(fstsA, annot_fst, by=c("siteB"="Locality")) %>% 
  rename_at(c("admix_groups","lat","lon", "HUC_6","EcoRegion"), funs( paste0(., "_B"))) %>% arrange(admix_groups_A) %>% 
  as.data.frame()

#fsts_out$admix_groups_A <- factor(fsts_out$admix_groups_A)

# append to site names based on admix groups:
head(fsts_out)
fsts_out <- fsts_out %>% 
  mutate(sitea = paste0(admix_groups_A,"_",siteA),
         siteb = paste0(admix_groups_B, "_", siteB))


# fsts_out %>% group_by(admix_groups_A) %>% select(siteA) %>% as.data.frame() %>% select(siteA) %>% distinct(siteA) %>% 
#   as.factor -> siteA_ord
# fsts_out %>% group_by(admix_groups_B) %>% select(siteB) %>% as.data.frame() %>% select(siteB) %>% distinct(siteB) %>% 
#   as.factor -> siteB_ord
# siteA_ord

siteA_ord <- c("sfa-cami", "tuo-clavey", "bear", "bear-grhn", "bear-sth2", "bear-stha", "bear-sthc", "deer-clec", "mfa-amec", "mfa-gasc", "mfa-todc", "mfy-oregck", "mfy-us-oh", "nfa-bunc", "nfa-euchds", "nfa-euchus", "nfa-indc", "nfa", "nfa-pond", "nfa-robr", "nfa-saic", "nfa-shic", "nfa-slar", "nfmfa-sc", "nfy", "nfy-slate-cgrav", "rub-lc-us", "rub-usph", "sfy-fallck", "sfy-loga", "sfy-misc", "sfy-rockck", "sfy-shadyck", "sfy-springck", "fea-beanck", "fea-spanish-bgulch", "fea-spanishck", "fea-spanish-silverck", "nff-poe", "eel", "mad-mapleck", "mad-redwoodck2", "mad-redwoodck", "mat-bearriver", "mat-lowermattole", "put-cold", "put-wildhorseck", "russ-mwsprngck", "sfeel-cedar", "smith-hurdygurdyck", "ssantiam", "trin-sftrinity-sandybar", "trin-sftrinity", "trin-tishtang", "vandz-dinsmore", "vandz-rootcreek", "salin-losburros", "sancarp-dutrack", "ala-arroyomocho", "paj-clearck", "paj-sanbenito", "soquel-ebsc")

siteB_ord <- c("sfeel-cedar", "sfy-fallck", "sfy-loga", "sfy-misc", "sfy-rockck", "sfy-shadyck", "sfy-springck", "smith-hurdygurdyck", "soquel-ebsc", "ssantiam", "trin-sftrinity", "trin-sftrinity-sandybar", "trin-tishtang", "tuo-clavey", "vandz-dinsmore", "vandz-rootcreek", "vandz-shanty", "bear-grhn", "bear-sth2", "bear-stha", "bear-sthc", "deer-clec", "eel", "fea-beanck", "fea-spanish-bgulch", "fea-spanishck", "fea-spanish-silverck", "mad-mapleck", "mad-redwoodck", "mad-redwoodck2", "mat-bearriver", "mat-lowermattole", "mfa-amec", "mfa-gasc", "mfa-todc", "mfy-oregck", "mfy-us-oh", "nfa", "nfa-bunc", "nfa-euchds", "nfa-euchus", "nfa-indc", "nfa-pond", "nfa-robr", "nfa-saic", "nfa-shic", "nfa-slar", "nff-poe", "nfmfa-sc", "nfy", "nfy-slate-cgrav", "paj-clearck", "paj-sanbenito", "put-cold", "put-wildhorseck", "rub-lc-us", "rub-usph", "russ-mwsprngck", "salin-losburros", "sancarp-dutrack", "sfa-cami", "bear") 

#anti_join(siteA_ord, siteB_ord, by=c("siteA"="siteB"))

fsts_out$siteA <- factor(fsts_out$siteA, levels=siteA_ord)
fsts_out$siteB <- factor(fsts_out$siteB, levels=siteB_ord)
#fsts_out$siteA <- factor(fsts_out$siteA, levels=fsts_out$siteA, labels=fsts_out$admix_groups_A)
#fsts_out$siteB <- factor(fsts_out$siteB, levels=fsts_out$siteB, labels=fsts_out$admix_groups_B)


# HEATMAP PLOT ------------------------------------------------------------

# plot
ggplot() + 
  geom_tile(data = fsts_out, aes(x=sitea, y=siteb, fill=fst_adj)) + ylab("") + xlab("")+
  theme_minimal(base_size = 8, base_family = "Roboto Condensed") +
  scale_fill_viridis("Fst Weighted", limit = c(0,.4)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+
  coord_fixed() + 
  #geom_text(data=melted_fst, aes(x=Var1, y=Var2, label = round(value, digits = 3)), color = "black", size = 1.2) +
  # add text
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    #legend.justification = c(1, 0),
    #legend.position = c(0.9, 0.1),
    legend.direction = "vertical")#+
    #legend.direction = "horizontal")+
  #guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
   #                            title.position = "top", title.hjust = 0.5))

#ggsave(filename = "figs/fst_matrix_heatmap_all_rabo_filt_100k.svg", width = 8, height=7, units = "in", dpi=600)
ggsave(filename = "figs/fst_matrix_heatmap_all_rabo_filt_100k.png", width = 8, height=7, units = "in", dpi=600)

# EVERYTHING AFTER THIS IS JUST WHO THE FUCK KNOWS... ---------------------

# tried to reorder things through here but no luck...or weird patterns. should be this hard but it is. Will need to deal with this later if it's worth it.


# MAKE MATRIX -------------------------------------------------------------

fst_mat <- fsts_out %>% select(siteA, siteB, fst_adj, admix_groups_A)
fst_mat$siteA <- as.factor(siteA)
fst_mat$admix_groups_A <- as.factor(fst_mat$admix_groups_A)
#fst_mat <- fst_mat[order(fsts_out$admix_groups_A),]

fst_mat <- fst_mat[,c(1, 2, 3)] %>% as.data.frame()

#library(spaa)
#fst_mat <- list2dist(dat = fst_mat)

# convert to a matrix for other stuff:
fst_mat <- fsts_out %>% select(siteA, siteB, fst_adj)
fst_mat <- with(fst_mat,tapply(fst_adj,list(siteA,siteB),"[[",1))

melted_fst <- melt(fst_mat, na.rm = T)

library(reshape2)
# # Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
cormat[upper.tri(cormat)] <- NA
return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
cormat[lower.tri(cormat)]<- NA
return(cormat)
}
fst_mat2 <- get_upper_tri(fst_mat)


# PLOT --------------------------------------------------------------------


# plot
ggplot() + 
  geom_tile(data = melted_fst, aes(x=Var1, y=Var2, fill=value)) + ylab("") + xlab("")+
  theme_minimal(base_size = 8, base_family = "Roboto Condensed") +
  scale_fill_viridis("Fst Weighted", limit = c(0,.4)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+
  coord_fixed() + 
  #geom_text(data=melted_fst, aes(x=Var1, y=Var2, label = round(value, digits = 3)), color = "black", size = 1.2) +
# add text
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.9, 0.1),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

#ggsave(filename = "figs/fst_matrix_heatmap_all_rabo_filt_100k.png", width = 8, height=7, units = "in", dpi=300)


# Calculate distances -----------------------------------------------------

# load the metadata and join to sites:
library(sf)
library(Imap)

sites <- st_read("data_output/sites_25k_n3.shp") %>% arrange(Localty)
sites$River <- tolower(sites$River)
sites$Site <- tolower(sites$Site)
sites$Localty <- tolower(sites$Localty)

# fix space:
unique(sites$Localty)
sites$Localty<-gsub(pattern = "[[:space:]]", replacement = "-", x = sites$Localty)

# Functions ---------------------------------------------------------------

# calculate the matrix of distances between each pairwise comparison (each set of sites).
ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
  # If triangle.to.replace="lower", replaces the lower triangle of a square matrix with its upper triangle.
  # If triangle.to.replace="upper", replaces the upper triangle of a square matrix with its lower triangle.
  
  if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
  if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
  else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
  else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
  m[tri] <- t(m)[tri]
  return(m)
}

# calculate the matrix of distances between each pairwise comparison (each set of sites).
GeoDistanceInMetresMatrix <- function(df.geopoints){
  # Returns a matrix (M) of distances between geographic points.
  # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
  # (df.geopoints$lat[j], df.geopoints$lon[j]).
  # The row and column names are given by df.geopoints$name.
  
  GeoDistanceInMetres <- function(g1, g2){
    # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
    # The 1st value in the returned vector is the distance between g1[[1]] and g2[[1]].
    # The 2nd value in the returned vector is the distance between g1[[2]] and g2[[2]]. Etc.
    # Each g1[[x]] or g2[[x]] must be a list with named elements "index", "lat" and "lon".
    # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
    DistM <- function(g1, g2){
      require("Imap")
      return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$lat, lon.1=g1$lon, lat.2=g2$lat, lon.2=g2$lon, units="m")))
    }
    return(mapply(DistM, g1, g2))
  }
  
  n.geopoints <- nrow(df.geopoints)
  
  # The index column is used to ensure we only do calculations for the upper triangle of points
  df.geopoints$index <- 1:n.geopoints
  
  # Create a list of lists
  list.geopoints <- by(df.geopoints[,c("index", "lat", "lon")], 1:n.geopoints, function(x){return(list(x))})
  
  # Get a matrix of distances (in metres)
  mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")
  
  # Set the row and column names
  rownames(mat.distances) <- df.geopoints$Locality
  colnames(mat.distances) <- df.geopoints$Locality
  
  return(mat.distances)
}


# calc dist in km
site_distmatrix<-round(GeoDistanceInMetresMatrix(sites[,c(26,6:7)]) / 1000, digits = 2)
head(site_distmatrix)

# add names
colnames(site_distmatrix)<-sites$Localty
rownames(site_distmatrix)<-sites$Localty

# Plot
image(1:nrow(site_distmatrix), 1:ncol(site_distmatrix), site_distmatrix, axes = FALSE, 
      xlab="", ylab="", col = terrain.colors(100))

axis(1, 1:nrow(site_distmatrix), rownames(site_distmatrix), cex.axis = 0.7, las=3, family="Roboto Condensed")
axis(2, 1:nrow(site_distmatrix), colnames(site_distmatrix), cex.axis = 0.5, las=1, family="Roboto Condensed")
#text(expand.grid(1:nrow(site_distmatrix), 1:ncol(site_distmatrix)), sprintf("%0.1f", site_distmatrix), cex=0.6, family="Roboto Condensed")
title("Euclidean Distance matrix (km) \n for RABO Sites", cex.main=0.8, family="Roboto Condensed")


# transform to longwise:
site_dist_df <- data.frame(site= t(combn(colnames(site_distmatrix), 2)), dist=t(site_distmatrix)[lower.tri(site_distmatrix)])

# rename/format
site_dist_df <- site_dist_df %>% 
  rename(siteA = site.1, siteB=site.2, dist_km = dist) %>% 
  mutate_at(.vars = c("siteA", "siteB"), .funs = as.character) %>% 
  mutate(sitepair = paste0(siteA, "_", siteB))


# JOIN DATA ---------------------------------------------------------------

final_fst <- left_join(fsts_out, site_dist_df, by="sitepair") %>% select(siteA.x:sitepair, dist_km) %>% 
  rename(siteA=siteA.x, siteB=siteB.x)

# SAVE OUT ----------------------------------------------------------------

save(final_fst, file = "data_output/fst/final_fst_25k_n3.rda") 


# PLOT FST WITH IMAGE -----------------------------------------------------

fsts_only <- fsts_out %>% select(siteA, siteB, fst_adj) %>% as.data.frame()
str(fsts_only)

tst <- with(fsts_only, sort(unique(c(as.character(siteA_ord),
                                     as.character(siteB_ord)))))
#tst <- with(fsts_only, sort(unique(c(as.character(siteA),
#                              as.character(siteB)))))
M <- array(0, c(length(tst), length(tst)), list(tst, tst))
i <- match(fsts_only$siteA, tst)
j <- match(fsts_only$siteB, tst)
M[cbind(i,j)] <- M[cbind(j,i)] <- fsts_only$fst_adj



values = data.frame(Eel_M = values.o$Eel_M, Tri_M = values.o$Tri_M, Tri_P = values.o$Tri_P, Sal_M = values.o$Sal_M, Sal_P = values.o$Sal_P, Rog_M = values.o$Rog_M, Rog_P = values.o$Rog_P, Sou_M = values.o$Sou_M, Nor_P = values.o$Nor_P, Sil_M = values.o$Sil_M, Sil_P = values.o$Sil_P, Puy_M = values.o$Puy_M, Puy_P = values.o$Puy_P, Noo_M = values.o$Noo_M, Noo_P = values.o$Noo_P, Sou_P = values.o$Sou_P)
rownames(fsts_values) <- siteA_ord
  fst_values$Sou_P <- NULL
  fst_values <- head(fst_values,-1)
  siteA_ord <- head(siteA_ord, -1)

n1<-nrow(values)-1
n2<-ncol(values)-1
key <- c(seq(from = 0.015, to = 0.22, by = (0.205/n1)))
values$key <- key
values <- rbind(values, key)
values <- as.matrix(values)
key <- round(key,2)
colv <- colorRampPalette(c("black", "white"))


# Plot
output.file = paste("plots/",spp, "_fst.pdf", sep = "")
pdf(file = output.file, height = 6.5/2.54, width = 6.5/2.54)
par(mar=c(5,5,5,5), cex= 0.3)
image(x=seq(0,n1+1), y=seq(0,n2+1), col=colv(1000), z=(-values), axes = F, xlab = "", ylab = "")
axis(1, at = c(0:(n1)), labels = order, las =2)
axis(2, at = c(0:(n1)), labels = order, las =2)
axis(1, at = c(n1+1), labels = c("key"), las =2)
axis(2, at = c(n1+1), labels = c("key"), las =2)
axis(3, at = c(0:(n1)), labels = key, las =2, cex.lab= 0.0001)
axis(4, at = c(0:(n1)), labels = key, las =2, cex.lab= 0.0001)
dev.off();




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

metadat<- metadat %>% 
  mutate(admix_groups = case_when(
    grepl("STAN|TUO|SFA", River) ~ "E", # southern siera
    grepl("ANTV|BEAR|DEER|MFA|MFY|NFA|NFMFA|NFY|SFY|RUB", River) ~ "NE", # northern sierra
    grepl("CHETCO|SFEEL|VANDZ|TRIN|MAT|KLAM|SSANTIAM|PUT|MAD|LAGUN|SUMPQUA|RUSS|SMITH|EEL", River) ~ "NW", # North Coast
    grepl("NFF|FEA", River) ~ "N-Fea", # feather
    grepl("PAJ|ALA|DRY|SOQUEL", River) ~ "W", # Central Coast
    grepl("SANCARP|SALIN", River) ~ "SW") # South Coast
  ) %>% 
  select(Seq, admix_groups, Locality, lat, lon: NHD_Tot_DA_sqkm, River, Site, EcoRegion)

# add McCartney Clades
metadat <- metadat %>%  
  mutate(mccartClades = case_when(
    grepl("sfa-cami", Locality) ~ "NE",
    grepl("stan-roseck|tuo-clavey", Locality) ~ "E",
    grepl("North-Feather", admix_groups) ~ "NE",
    grepl("North-East", admix_groups) ~ "NE",
    grepl("North-West", admix_groups) ~ "NW",
    grepl("South-West", admix_groups) ~"SW",
    grepl("West", admix_groups) ~ "W")
  ) %>%   select(Seq, admix_groups, mccartClades, Locality, lat, lon: NHD_Tot_DA_sqkm, River, Site, EcoRegion)

# set order in way that you want
ords_admix_grps <- c("E", "NE", "N-Fea","NW", "SW", "W")

metadat$admix_groups <- factor(metadat$admix_groups, levels = ords_admix_grps)
levels(metadat$admix_groups)


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
annot_fst <- annot %>% select(admix_groups, Locality) #lat, lon, HUC_6, EcoRegion)

fstsA <- left_join(fsts, annot_fst, by=c("siteA"="Locality")) %>% 
  rename_at(c("admix_groups"), funs( paste0(., "_A")))
  #rename_at(c("admix_groups","lat","lon", "HUC_6","EcoRegion"), funs( paste0(., "_A")))

fsts_out <- left_join(fstsA, annot_fst, by=c("siteB"="Locality")) %>% 
  rename_at(c("admix_groups"), funs( paste0(., "_B"))) %>% arrange(admix_groups_A) %>% 
  as.data.frame()

# append to site names based on admix groups:
fsts_out <- fsts_out %>% 
  mutate(sitea = paste0(admix_groups_A,"_",siteA),
         siteb = paste0(admix_groups_B, "_", siteB))

# could add locality NUMBER here instead of site name

head(fsts_out)

rm(fstsA)

# siteA_ord <- c("sfa-cami", "tuo-clavey", "bear", "bear-grhn", "bear-sth2", "bear-stha", "bear-sthc", "deer-clec", "mfa-amec", "mfa-gasc", "mfa-todc", "mfy-oregck", "mfy-us-oh", "nfa-bunc", "nfa-euchds", "nfa-euchus", "nfa-indc", "nfa", "nfa-pond", "nfa-robr", "nfa-saic", "nfa-shic", "nfa-slar", "nfmfa-sc", "nfy", "nfy-slate-cgrav", "rub-lc-us", "rub-usph", "sfy-fallck", "sfy-loga", "sfy-misc", "sfy-rockck", "sfy-shadyck", "sfy-springck", "fea-beanck", "fea-spanish-bgulch", "fea-spanishck", "fea-spanish-silverck", "nff-poe", "eel", "mad-mapleck", "mad-redwoodck2", "mad-redwoodck", "mat-bearriver", "mat-lowermattole", "put-cold", "put-wildhorseck", "russ-mwsprngck", "sfeel-cedar", "smith-hurdygurdyck", "ssantiam", "trin-sftrinity-sandybar", "trin-sftrinity", "trin-tishtang", "vandz-dinsmore", "vandz-rootcreek", "salin-losburros", "sancarp-dutrack", "ala-arroyomocho", "paj-clearck", "paj-sanbenito", "soquel-ebsc")
# 
# siteB_ord <- c("sfeel-cedar", "sfy-fallck", "sfy-loga", "sfy-misc", "sfy-rockck", "sfy-shadyck", "sfy-springck", "smith-hurdygurdyck", "soquel-ebsc", "ssantiam", "trin-sftrinity", "trin-sftrinity-sandybar", "trin-tishtang", "tuo-clavey", "vandz-dinsmore", "vandz-rootcreek", "vandz-shanty", "bear-grhn", "bear-sth2", "bear-stha", "bear-sthc", "deer-clec", "eel", "fea-beanck", "fea-spanish-bgulch", "fea-spanishck", "fea-spanish-silverck", "mad-mapleck", "mad-redwoodck", "mad-redwoodck2", "mat-bearriver", "mat-lowermattole", "mfa-amec", "mfa-gasc", "mfa-todc", "mfy-oregck", "mfy-us-oh", "nfa", "nfa-bunc", "nfa-euchds", "nfa-euchus", "nfa-indc", "nfa-pond", "nfa-robr", "nfa-saic", "nfa-shic", "nfa-slar", "nff-poe", "nfmfa-sc", "nfy", "nfy-slate-cgrav", "paj-clearck", "paj-sanbenito", "put-cold", "put-wildhorseck", "rub-lc-us", "rub-usph", "russ-mwsprngck", "salin-losburros", "sancarp-dutrack", "sfa-cami", "bear") 

#anti_join(siteA_ord, siteB_ord, by=c("siteA"="siteB"))



# CONVERT TO MATRIX -------------------------------------------------------

fsts_only <- fsts_out %>% select(sitea, siteb, fst_adj) %>% as.data.frame()

fst <- with(fsts_only, sort(unique(c(as.character(sitea),
                                     as.character(siteb)))))

fst_M <- array(0, c(length(fst), length(fst)), list(fst, fst))
i <- match(fsts_only$sitea, fst)
j <- match(fsts_only$siteb, fst)
fst_M[cbind(i,j)] <- fst_M[cbind(j,i)] <- fsts_only$fst_adj
fst_M_df <- fst_M %>% as.data.frame(fst_M)
#fst_M_df$sites <- row.names(fst_M)
#fst_M_df <- fst_M_df %>% select(sites, everything())

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

# get upper triangle
fst_M_upper <- get_upper_tri(fst_M)


# MELT FOR LONG FORMAT ----------------------------------------------------

library(reshape2)
melted_fst <- melt(fst_M, na.rm = T)

melted_fst_upper <- melt(fst_M_upper, na.rm=T)

# GGPLOT --------------------------------------------------------------------

# plot
ggplot() + 
  geom_tile(data = melted_fst_upper, aes(x=Var1, y=Var2, fill=value)) + ylab("") + xlab("")+
  theme_minimal(base_size = 8, base_family = "Roboto Condensed") +
  scale_fill_viridis("Fst Weighted", limit = c(0,.4)) +
  scale_x_discrete(position = "top") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 0))+
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

ggsave(filename = "figs/fst_matrix_heatmap_all_rabo_filt_100k.png", width = 8, height=7, units = "in", dpi=300)


# IMAGE HEAT --------------------------------------------------------------

image(1:nrow(M), 1:ncol(M), M, axes = FALSE, 
      xlab="", ylab="", col = viridis(20))

axis(1, 1:nrow(M), rownames(M), cex.axis = 0.7, las=3, family="Roboto Condensed")
axis(2, 1:nrow(M), colnames(M), cex.axis = 0.5, las=1, family="Roboto Condensed")
axis(3, 1:nrow(M), colnames(M), cex.axis = 0.5, las=2, family="Roboto Condensed", cex.lab=0.001)


library(fields)

pdf(file = paste0("figs/fst_matrix_heatmap_all_rabo_filt_100k.png"), height = 6.5/2.54, width = 6.5/2.54)
par(mar=c(3,8,8,3))
image.plot(1:nrow(M), 1:ncol(M), M, axes = FALSE, 
           xlab="", ylab="", col = viridis(200))

#axis(1, 1:nrow(M), rownames(M), cex.axis = 0.7, las=2, family="Roboto Condensed", cex.lab=0.001)
axis(2, 1:nrow(M), colnames(M), cex.axis = 0.7, las=2, family="Roboto Condensed", cex.lab=0.001)
axis(3, 1:nrow(M), colnames(M), cex.axis = 0.7, las=2, family="Roboto Condensed", cex.lab=0.001)
#axis(4, 1:nrow(M), colnames(M), cex.axis = 0.5, las=1, family="Roboto Condensed")
dev.off()



# Calculate distances -----------------------------------------------------

# load the metadata and join to sites:
library(sf)
library(Imap)

sites <- st_read("data_output/sites_all_rabo_filt10_1_100k.shp") %>% arrange(Localty)
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


# calc dist in km (need the grouping ID and the lat lon cols)
site_distmatrix<-round(GeoDistanceInMetresMatrix(sites[,c(28,8:9)]) / 1000, digits = 2)
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

final_fst_dist <- left_join(fsts, site_dist_df, by="sitepair") %>% select(siteA.x:sitepair, dist_km, fst_adj) %>% 
  rename(siteA=siteA.x, siteB=siteB.x)

# join with admix groups:
final_fst_dist <- left_join(final_fst_dist, fsts_out[,c("sitepair", "sitea", "siteb")], by="sitepair")

# SAVE OUT ----------------------------------------------------------------

save(final_fst_dist, file = "data_output/fst/final_fst_all_rabo_filt_100k.rda") 

# LOOK AT PLOTTING BOTH? --------------------------------------------------



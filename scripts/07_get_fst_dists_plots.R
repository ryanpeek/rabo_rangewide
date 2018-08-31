
# GET FST Data ------------------------------------------------------------

# first we copy over from either farmer or using scp:
# scp -P 2022 rapeek@agri.cse.ucdavis.edu:/home/rapeek/projects/rangewide/pop_gen/results_fst/all_global_fst.txt .

# Load Libraries ----------------------------------------------------------

library(tidyverse)
library(stringr)
library(here)
library(viridis)

# skip steps 00-07 and load data at step 08.

# 00. DATA PREP -----------------------------------------------------------

## GET DATA

# make spaces to tab delimited
# :%s/ /\t/g
# :%s/_100k\.folded\.finalFSTout//g
# need pop 1, pop 2, fst for perlscript

dat <- read_tsv("data_output/fst/all_rabo_100k_folded_adj_fst.txt", 
                col_names = c("fst_unad", "fst_adj", "filenames"))

## Fix Text

dat$filenames <- gsub(x = dat$filenames, pattern = "_100k.folded.finalFSTout", replacement = "")

## Separate and Cleanup
fsts <- dat %>% 
  separate(col = filenames, into = c("siteA", "siteB"), sep = "\\.") %>% 
  mutate(sitepair = paste0(siteA, "_", siteB))

# quick tally
fsts %>% group_by(siteA) %>% tally
fsts %>% group_by(siteA) %>% tally

## MERGE WITH METADATA

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
    grepl("N-Fea", admix_groups) ~ "NE",
    grepl("NE", admix_groups) ~ "NE",
    grepl("NW", admix_groups) ~ "NW",
    grepl("SW", admix_groups) ~"SW",
    grepl("W", admix_groups) ~ "W")
  ) %>%   select(Seq, admix_groups, mccartClades, Locality, lat, lon: NHD_Tot_DA_sqkm, River, Site, EcoRegion)

# set order in way that you want
ords_admix_grps <- c("E", "NE", "N-Fea","NW", "SW", "W")

metadat$admix_groups <- factor(metadat$admix_groups, levels = ords_admix_grps)
levels(metadat$admix_groups)

## GET AND MERGE WITH BAMFILES

bamfile <- "all_rabo_filt_100k"

# Get ID and pop info for each individual from bamlists
bams <- read.table(paste0("data_output/bamlists/",bamfile, "_thresh.bamlist"),stringsAsFactors = F, header = F)
bams$V2 <- sub('\\..*$', '', basename(bams$V1)) # remove the path and file extension

annot <- left_join(bams, metadat, by=c("V2"="Seq")) %>% select(-V1) # join with the metadata

annot %>% group_by(Locality) %>% tally
annot %>% distinct(Locality) %>% arrange

# get only the sites for pairings
annot <- annot %>% distinct(Locality, .keep_all = T)
annot$Locality <- tolower(annot$Locality)

# 01. JOIN DATA ----------------------------------------------------

# get LOCALITY ID and make label column
loc_ID  <- read_csv("data_output/table_site_localities_clades.csv") %>% 
  mutate(Locality = tolower(Locality),
         Locality = gsub(pattern="deer-clearck", replacement = "deer-clec", Locality))
  

# join with annot
annot2 <- left_join(annot, loc_ID[,c(1:2)], by="Locality")


annot_fst <- annot2 %>% select(admix_groups, Locality, siteID) #lat, lon, HUC_6, EcoRegion)

# now join with fst data:
fstsA <- left_join(fsts, annot_fst, by=c("siteA"="Locality")) %>% 
  rename_at(c("admix_groups", "siteID"), funs( paste0(., "_A")))
  #rename_at(c("admix_groups","lat","lon", "HUC_6","EcoRegion"), funs( paste0(., "_A")))

fsts_out <- left_join(fstsA, annot_fst, by=c("siteB"="Locality")) %>% 
  rename_at(c("admix_groups", "siteID"), funs( paste0(., "_B"))) %>% arrange(admix_groups_A) %>% 
  as.data.frame()

# append to site names based on admix groups:
fsts_out <- fsts_out %>% 
  mutate(sitea = paste0(admix_groups_A,"_",siteA),
         siteb = paste0(admix_groups_B, "_", siteB),
         siteida = paste0(admix_groups_A, "-",siteID_A),
         siteidb = paste0(admix_groups_B, "-",siteID_B))

# remove unused objects
rm(fstsA, bams, dat, loc_ID, annot, annot2, annot_fst, fsts)

# 02. CONVERT TO MATRIX --------------------------------------------------

# select only sites and fsts value
fsts_only <- fsts_out %>% select(siteida, siteidb, fst_adj) %>% as.data.frame()

# get list of site names for row/cols
fst <- with(fsts_only, sort(unique(c(as.character(siteida),
                                     as.character(siteidb)))))

# set up a matrix array based on length/width of data
fst_M <- array(0, c(length(fst), length(fst)), list(fst, fst))
i <- match(fsts_only$siteida, fst) # fill vectors with matching names
j <- match(fsts_only$siteidb, fst)
fst_M[cbind(i,j)] <- fst_M[cbind(j,i)] <- fsts_only$fst_adj # now add value for matrix
fst_M_df <- fst_M %>% as.data.frame(fst_M) # back to a wide dataframe if you want

#fst_M_df$sites <- row.names(fst_M)
#fst_M_df <- fst_M_df %>% select(sites, everything())


# 03. MELT TO LONG FORMAT -------------------------------------------------


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

# melt data from wide to long
library(reshape2)
melted_fst <- melt(fst_M, na.rm = T)

melted_fst_upper <- melt(fst_M_upper, na.rm=T)


# 04. GGPLOT OF PAIRWISE MATRIX -------------------------------------------

# plot
ggfstmat <- ggplot() + 
  geom_tile(data = melted_fst_upper, aes(x=Var1, y=Var2, fill=(value/(1-value)))) + ylab("") + xlab("")+
  theme_minimal(base_size = 9, base_family = "Roboto Condensed") +
  scale_fill_viridis("Fst Weighted") + #limit = c(0,.4)) +
  scale_x_discrete(position = "top") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 0))+
  coord_fixed() + 
  #geom_text(data=melted_fst, aes(x=Var1, y=Var2, label = round(value, digits = 3)), color = "black", size = 1.2) +
# add text
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 6),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 6),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(1, 0.65),
    legend.direction = "vertical")
  #guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
  #                             title.position = "top", title.hjust = 0.5))
ggfstmat

ggsave(filename = "figs/fst_matrix_heatmap_all_rabo_filt_100k.png", width = 8, height=7, units = "in", dpi=300)

# IMAGE HEAT --------------------------------------------------------------

# alternate option for plotting
image(1:nrow(fst_M), 1:ncol(fst_M), fst_M, axes = FALSE, 
      xlab="", ylab="", col = viridis(20))

axis(1, 1:nrow(fst_M), rownames(fst_M), cex.axis = 0.7, las=3, family="Roboto Condensed")
axis(2, 1:nrow(fst_M), colnames(fst_M), cex.axis = 0.5, las=1, family="Roboto Condensed")
axis(3, 1:nrow(fst_M), colnames(fst_M), cex.axis = 0.5, las=2, family="Roboto Condensed", cex.lab=0.001)


library(fields)

#pdf(file = paste0("figs/fst_matrix_heatmap_all_rabo_filt_100k.png"), height = 6.5/2.54, width = 6.5/2.54)
par(mar=c(3,8,8,3))
image.plot(1:nrow(fst_M), 1:ncol(fst_M), fst_M, axes = FALSE, 
           xlab="", ylab="", col = viridis(200))

#axis(1, 1:nrow(fst_M), rownames(fst_M), cex.axis = 0.7, las=2, family="Roboto Condensed", cex.lab=0.001)
axis(2, 1:nrow(fst_M), colnames(fst_M), cex.axis = 0.7, las=2, family="Roboto Condensed", cex.lab=0.001)
axis(3, 1:nrow(fst_M), colnames(fst_M), cex.axis = 0.7, las=2, family="Roboto Condensed", cex.lab=0.001)
#axis(4, 1:nrow(fst_M), colnames(fst_M), cex.axis = 0.5, las=1, family="Roboto Condensed")
#dev.off()

# 05. Calculate distances -----------------------------------------------------

# load the metadata and join to sites:
library(sf)
library(Imap)

sites <- st_read("data_output/sites_all_rabo_filt_100k.shp") %>% arrange(Localty)
sites$River <- tolower(sites$River)
sites$Site <- tolower(sites$Site)
sites$Localty <- tolower(sites$Localty)

# fix space:
unique(sites$Localty)
sites$Localty<-gsub(pattern = "[[:space:]]", replacement = "-", x = sites$Localty)
sites$Localty<-gsub("deer-clearck", replacement = "deer-clec", x = sites$Localty)

## Functions 

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

# Plot distances in matrix
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

# 06. JOIN DATA ---------------------------------------------------------------

final_fst_dist <- left_join(fsts_out, site_dist_df, by="sitepair") %>% 
  select(siteA.x:sitepair, siteida, siteidb, dist_km, fst_adj) %>% 
  rename(siteA=siteA.x, siteB=siteB.x)

glimpse(final_fst_dist)

# 07. SAVE OUT -----------------------------------------------------------

#save(melted_fst_upper, melted_fst, fst_M, ggfstmat, final_fst_dist, file = "data_output/fst/final_fst_all_rabo_filt_100k.rda") 

# 08. FINAL COMBINE PLOTS ------------------------------------------------

load("data_output/fst/final_fst_all_rabo_filt_100k.rda")

fsts_dist <- final_fst_dist %>% select(siteida, siteidb, dist_km) %>% as.data.frame()

fst <- with(fsts_dist, sort(unique(c(as.character(siteida),
                                     as.character(siteidb)))))
# fst_dist
fst_Mdist <- array(0, c(length(fst), length(fst)), list(fst, fst))
i <- match(fsts_dist$sitea, fst)
j <- match(fsts_dist$siteb, fst)
fst_Mdist[cbind(i,j)] <- fst_Mdist[cbind(j,i)] <- fsts_dist[,3]
fst_Mdist_df <- fst_Mdist %>% as.data.frame(fst_Mdist)

# fst_fst
fsts_gen <- final_fst_dist %>% select(siteida, siteidb, fst_adj) %>% as.data.frame()

fst <- with(fsts_gen, sort(unique(c(as.character(siteida),
                                     as.character(siteidb)))))

# fst_dist
fst_Mfst <- array(0, c(length(fst), length(fst)), list(fst, fst))
i <- match(fsts_gen$siteida, fst)
j <- match(fsts_gen$siteidb, fst)
fst_Mfst[cbind(i,j)] <- fst_Mfst[cbind(j,i)] <- fsts_gen[,3]
fst_Mfst_df <- fst_Mfst %>% as.data.frame(fst_Mdist)

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
fst_M_upper <- get_upper_tri(fst_Mfst) # genetic dist
fst_M_lower <- get_lower_tri(fst_Mdist) # geographic dist

# melt for plots
library(reshape2)
melted_fst_upper <- melt(fst_M_upper, na.rm=T) # genetic
melted_dist_lower <- melt(fst_M_lower, na.rm=T) # geographic


# 08a. Make point plot --------------------------------------------------------------

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# make admix cols with tidyr
fst_summary <- final_fst_dist %>% 
  separate(siteida, c("admixA"), sep = "-", extra = "drop", remove = F) %>% 
  separate(siteidb, c("admixB"), sep = "-", extra="drop", remove=F) %>% 
  select(siteida, siteidb, admixA, admixB, dist_km, fst_adj) %>% 
  mutate(within = if_else(admixA == admixB, "Y", "N"),
         admix_pair = paste0(admixA, "_", admixB))

# PLOT

#plotly::ggplotly(

(pt_fst <- ggplot(data=fst_summary, aes(x=dist_km, y=(fst_adj/(1-fst_adj)),
                             text=paste0("siteA: ", siteida, " <br> siteB: ", siteidb),
                             color=within)) + #color=as.factor(admix_pair))) + 
  geom_point(size=4, alpha=0.7) +
  #ggrepel::geom_label_repel(data=fst_summary, aes(x=dist_km, y=fst_adj, label=admix_pair)) +
  scale_color_manual("", values = c("Y"=cbbPalette[3], "N"=cbbPalette[2]), labels=c("Between Clade", "Within Clade")) +
  
  theme_bw(base_family = "Helvetica", base_size = 9) +
  labs(#title=expression(paste("Mean F" ["ST"], " vs Mean Distance (km)")),
    y=expression(paste("Mean F" ["ST"], " / (1 - Mean F" ["ST"],")")),
    x="Euclidean Distance (km)") +
  theme(legend.position = c(0.75, 0.18)))

#)

#ggsave(filename = "figs/fst_vs_dist_by_clade.png", width = 3.42, height= 2.9, units = "in", dpi = 600, scale=1.3)
ggsave(filename = "figs/fst_vs_dist_by_clade_facet_admixB.png", width = 6, height= 5, units = "in", dpi = 300)


# 09. COMBINE WITH COWPLOT ------------------------------------------------
library(cowplot)

pt_fst
ggfstmat

#extract legend from a single plot
fst_legend <- get_legend(ggfstmat)

# make into one plot
(combined_fst_cowplot <- ggdraw() +
  draw_plot(ggfstmat + theme(legend.justification = c(0, 0.6),
                             legend.position = c(0.85, 0.75)), 0, 0, 1, 1) + 
              #theme(legend.justification = "bottom"), 0, 0, 1, 1) +
  draw_plot(pt_fst, 0.54, 0.02, 0.43, 0.42, scale = 1.1) +
  draw_plot_label(c("A", "B"), c(0.07, 0.52), c(0.98, 0.46), size = 12))

# save it
save_plot(combined_fst_cowplot , filename = "figs/combined_fst_cowplot.png", base_width = 8, base_height = 7,
          #base_aspect_ratio = 1.5, 
          dpi = 300)


# OTHER STUFF -------------------------------------------------------------
# summarize
# fst_summaryA <- fst_summary %>% 
#   group_by(admixA, within) %>% 
#   summarize(fst_admixA = mean(fst_adj, na.rm=T),
#          dist_admixA = mean(dist_km, na.rm = T)) 
# 
# fst_summaryB <- fst_summary %>% 
#   group_by(admixB, within) %>% 
#   summarize(fst_admixB = mean(fst_adj, na.rm=T),
#          dist_admixB = mean(dist_km, na.rm = T))

# ggplot() + geom_point(data=fst_summaryA, aes(x=dist_admixA, y=fst_admixA, fill=fst_admixA), pch=21, size=4, alpha=0.7) + 
#   ggrepel::geom_label_repel(data=fst_summaryA, aes(x=dist_admixA, y=fst_admixA, label=admixA)) + 
#   scale_fill_viridis() +
#   geom_point(data=fst_summaryB, aes(x=dist_admixB, y=fst_admixB, fill=fst_admixB), pch=21, size=4, alpha=0.8) + 
#   ggrepel::geom_label_repel(data=fst_summaryB, aes(x=dist_admixB, y=fst_admixB, label=admixB), fill="cornsilk") + 
#   scale_fill_viridis()
# ggsave(filename = "figs/fst_vs_eucDist_by_admix_groupsAB.png", width = 7, height = 5, units = "in", dpi=300)










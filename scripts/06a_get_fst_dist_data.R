
# GET FST Data ------------------------------------------------------------

# first we copy over from either farmer or using scp:
# scp -P 2022 rapeek@agri.cse.ucdavis.edu:/home/rapeek/projects/rangewide/pop_gen/results_fst/all_global_fst.txt .

# Load Libraries ----------------------------------------------------------

library(tidyverse)
library(stringr)

# GET DATA ----------------------------------------------------------------

# create in bash in farm:
# tail -n+2 -q *fst | cat > all_global_fst_DATA.txt
# visual block delete space 
# sed out crap
# :%s/\.folded\.fst\.gz//g
# :%s/ST\.Unweight\[//g
# :%s/ssuming\ \.fst\.gz\ file\:\ results_fst\///g

# finally, fix file that's messed up: missing (calav-esperck.nfa-slar_25k.folded.fst.idx), manually deleted line.


dat <- read_csv("data_output/fst/all_global_fst_25k_n3.txt", col_names = "alldat")

# now pull out first line and add as a new col
datFiles <- dat[seq(from=1, nrow(dat), by=2),]
datFst <- dat[seq(from=2, nrow(dat), by=2),]

fsts <- bind_cols(datFiles, datFst) %>% 
  rename(fst=alldat1, filenames=alldat)

rm(datFiles, datFst)

# FIX WEIRD TXT -----------------------------------------------------------

fsts$filenames <- gsub(x = fsts$filenames, pattern = "_25k", replacement = "")

# SEPARATE COLS -----------------------------------------------------------

# separate and cleanup
fsts2 <- fsts %>% 
  separate(col = filenames, into = c("siteA", "siteB"), sep = "\\.") %>% 
  separate(col = fst, into = c("no_obs_loci", "fst_unweight", "fst_weight"), sep=" ") %>% 
  mutate(no_obs_loci = gsub(no_obs_loci, pattern="nObs:", replacement=""),
         fst_unweight = gsub(fst_unweight, pattern=":", replacement=""),
         fst_weight = gsub(fst_weight, pattern="Fst.Weight:", replacement="")) %>% 
  mutate_at(c("no_obs_loci", "fst_unweight", "fst_weight"), .funs = as.numeric) %>%
  mutate(sitepair = paste0(siteA, "_", siteB))



# FILTER OUT SITES WITH <3 samples ----------------------------------------

# 25k samples:
# filtout <- c("cache-wfsulphurck", "fea-spanish-rockck", "fea-spanish-wapaunsie",
#              "nfa-nfnfa", "nfy-slate", "sancarp-sancarpoforock", "sfeel", "stan-roseck")

#fsts_out <- fsts2 %>% filter(!siteA %in% filtout, !siteB %in% filtout)

#head(fsts_out)
fsts_out <- fsts2

# convert to a matrix for other stuff:
fst_mat <- fsts_out %>% select(siteA, siteB, fst_weight)

fst_mat <- with(fst_mat,tapply(fst_weight,list(siteA,siteB),"[[",1)) 

#write.csv(fst_mat, file = "data_output/fst/fst_matrix_25k_n3.csv")

# fst_long <- select(fsts_out, fst_weight, siteA, siteB) %>% 
#   gather(pair, site, -fst_weight) %>% group_by(site) %>% 
#   summarize(fst_mean=mean(fst_weight)) %>% rename(sites = site)


# TEST PLOTS: POINTS ------------------------------------------------------

library(viridis)

ggplot() + 
  geom_point(data=fsts_out, aes(y=fst_weight/(1-fst_weight), x=sitepair, color=fst_weight)) + 
  theme_bw(base_family = "Roboto Condensed") + 
  scale_color_viridis_c() + 
  theme(axis.text.x = element_blank())

library(plotly)
ggplotly(
  ggplot() + 
    geom_point(data=fsts_out, aes(y=fst_weight/(1-fst_weight), x=sitepair, color=fst_weight)) + 
    theme_bw(base_family = "Roboto Condensed") + 
    scale_color_viridis_c() + 
    theme(axis.text.x = element_blank())
)


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


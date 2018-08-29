
# GET ADMIX AND PLOT ------------------------------------------------------


# Load Libraries ----------------------------------------------------------

library(viridis)
library(tidyverse)
library(scales)

# GET BAMLIST AND METADATA ------------------------------------------------

bamfile <- "all_rabo_100k"

# Get ID and pop info for each individual from bamlists
bams <- read.table(paste0("data_output/bamlists/",bamfile, "_thresh.bamlist"),stringsAsFactors = F, header = F)
bams$V2 <- sub('\\..*$', '', basename(bams$V1)) # remove the path and file extension

# get metadat
metadat <- read_rds(path = "data_output/rapture_metadata_rabo_quant.rds")
# fix trinity spaces
metadat$Locality<-tolower(gsub(pattern = "[[:space:]]", replacement = "-", x = metadat$Locality))
# fix deer-clearck/ deer-clec
metadat$Locality <- gsub(pattern="deer-clearck", replacement = "deer-clec", x=metadat$Locality)


metadat<- metadat %>% 
  mutate(admix_groups = case_when(
    grepl("STAN|TUO|SFA", River) ~ "East", # southern siera
    grepl("ANTV|BEAR|DEER|MFA|MFY|NFA|NFMFA|NFY|SFY|RUB", River) ~ "North-East", # northern sierra
    grepl("CHETCO|SFEEL|VANDZ|TRIN|MAT|KLAM|SSANTIAM|PUT|MAD|LAGUN|SUMPQUA|RUSS|SMITH|EEL", River) ~ "North-West", # North Coast
    grepl("NFF|FEA", River) ~ "North-Feather", # feather
    grepl("PAJ|ALA|DRY|SOQUEL", River) ~ "West", # Central Coast
    grepl("SANCARP|SALIN", River) ~ "South-West") # South Coast
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
)

# set order in way that you want
#ords_admix_grps <- c("North-Feather", "North-East","East", "North-West", "South-West", "West")

ords_admix_grps <- c("East", "North-East", "North-Feather","North-West", "South-West", "West")

metadat$admix_groups <- factor(metadat$admix_groups, levels = ords_admix_grps)
levels(metadat$admix_groups)

# join together
annot <- left_join(bams, metadat, by=c("V2"="Seq")) %>% select(-V1)  %>% # join with the metadata
  mutate(annotID = row_number()) # original order of bamlist

# check for duplicates?
#annot[duplicated(annot$V2),]
#metadat[duplicated(metadat$Seq),]

# add id and trim
annot <- annot %>% 
  arrange(admix_groups, Locality) %>% # now arrange and number localities by admix groups
  mutate(siteID=row_number(),
         Locality=tolower(Locality)) %>% 
  select(siteID, annotID, admix_groups, mccartClades, Locality:HUC_10,county, EcoRegion) %>% 
  arrange(annotID)

# GGPLOT VERSION ----------------------------------------------------------

# add color palette: 
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#CFCFCF")

# the palette assignment to match the map
# scale_fill_manual(values = c("East"=cbbPalette[1], 
# "North-East"=cbbPalette[2], 
# "North-West"=cbbPalette[3],
# "North-Feather"=cbbPalette[4],
# "West"=cbbPalette[5], 
# "South-West"=cbbPalette[6])) +

# K3 ----------------------------------------------------------------------

k <- 3

# Read inferred admixture proportions
q<-read.table(paste0("data_output/admix/", bamfile, "_k", k,"_admix.qopt"))

admix <- bind_cols(q, annot)

admix_gg <- admix %>% gather(clustname, "clust", V1:V3) %>% 
  arrange(clustname, clust) %>% #group_by(clustname) %>% 
  mutate(
    clustname = fct_reorder2(factor(clustname, levels = c("V1", "V2", "V3")), clustname, clust),
    annotID = fct_reorder2(as.factor(annotID), clustname, clust),
    admix_groups = fct_reorder2(admix_groups, clustname, clust))

ords_admix_grps <- c("East", "North-East", "North-Feather","North-West", "South-West", "West")
admix_gg$admix_groups <- factor(admix_gg$admix_groups, levels = ords_admix_grps)
levels(admix_gg$admix_groups)

# replace V with Clust
admix_gg$clustname <- gsub("V", "Clust", admix_gg$clustname)

# the plot
(plotk3 <- ggplot(data = admix_gg,
                  aes(x = annotID , y = clust, fill = clustname)) + 
   geom_col(position = "fill", width=1, show.legend = F) + # switch legend off for print/save
   scale_x_discrete(expand=c(0,0)) +
   scale_y_continuous(expand=c(0,0), labels = percent, breaks = seq(0,1,0.1)) +
   #scale_fill_viridis_d("Cluster") + 
   scale_fill_manual(values = c("Clust1"=cbbPalette[1], # East
                                 #"Clust3"=cbbPalette[2], # North-East
                                 "Clust2"=cbbPalette[5], #North-West
                                 "Clust3"=cbbPalette[4]))+ # North-Feather
                                 #"Clust2"=cbbPalette[5])) + # West / South-West
   labs(title=paste0("NGSadmix k=",k)) +
   theme_classic(base_size = 8) +
   theme(axis.text.x = element_text(angle = 80, hjust = 1, size = 7)) +
   labs(x="", y="fraction ancestry") +
   facet_grid(.~ admix_groups, scales = "free_x") +
   theme(strip.background = element_blank(), 
         panel.spacing = unit(0, "lines"),
         panel.border = element_rect(fill = NA, color = "gray40"),
         axis.text.x=element_blank(),
         axis.ticks.x = element_blank()))

#ggsave(filename = paste0("figs/admix/",bamfile, "_k",k, "_admix.png"), width = 4.5, height = 2.7, scale = 1.3, units = "in", dpi = 300)


# K4 ----------------------------------------------------------------------

k <- 4

# Read inferred admixture proportions
q<-read.table(paste0("data_output/admix/", bamfile, "_k", k,"_admix.qopt"))

# bind with metadat
admix <- bind_cols(q, annot)

admix_gg <- admix %>% gather(clustname, "clust", V1:V4) %>% 
  arrange(clustname, clust) %>% #group_by(clustname) %>% 
  mutate(
    clustname = fct_reorder2(factor(clustname, levels = c("V1", "V2", "V3", "V4")), clustname, clust),
    annotID = fct_reorder2(as.factor(annotID), clustname, clust),
    admix_groups = fct_reorder2(admix_groups, clustname, clust))

ords_admix_grps <- c("East", "North-East", "North-Feather","North-West", "South-West", "West")
admix_gg$admix_groups <- factor(admix_gg$admix_groups, levels = ords_admix_grps)
levels(admix_gg$admix_groups)

# replace V with Clust
admix_gg$clustname <- gsub("V", "Clust", admix_gg$clustname)

(plotk4 <- ggplot(data = admix_gg,
                  aes(x = annotID , y = clust, fill = clustname)) + 
    geom_col(position = "fill", width=1, show.legend = F) + # turn off legend
    scale_x_discrete(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0), labels = percent, breaks = seq(0,1,0.1)) +
    #scale_fill_viridis_d("Cluster") + 
    scale_fill_manual(values = c("Clust2"=cbbPalette[1], # East
                                 #"Clust5"=cbbPalette[2], # North-East
                                 "Clust4"=cbbPalette[3], #North-West
                                 "Clust1"=cbbPalette[4], # North-Feather
                                 "Clust3"=cbbPalette[5])) + # West / South-West
    labs(title=paste0("NGSadmix k=",k)) +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 80, hjust = 1, size = 7)) +
    labs(x="", y="fraction ancestry") +
    facet_grid(.~ admix_groups, scales = "free_x") +
    theme(strip.background = element_blank(), 
          panel.spacing = unit(0, "lines"),
          panel.border = element_rect(fill = NA, color = "gray40"),
          axis.text.x=element_blank(),
          axis.ticks.x = element_blank()))

#ggsave(filename = paste0("figs/admix/",bamfile, "_k",k, "_admix.png"), width = 4.5, height = 2.7, scale = 1.3, units = "in", dpi = 300)
# width = 9, height = 6, units = "in", res = 300


# K5 ----------------------------------------------------------------------

k <- 5

# Read inferred admixture proportions
q<-read.table(paste0("data_output/admix/", bamfile, "_k", k,"_admix.qopt"))

# bind with metadat
admix <- bind_cols(q, annot)

admix_gg <- admix %>% gather(clustname, "clust", V1:V5) %>% 
  arrange(clustname, clust) %>% #group_by(clustname) %>% 
  mutate(
    clustname = fct_reorder2(factor(clustname, levels = c("V1", "V2", "V3", "V4", "V5")), clustname, clust),
    annotID = fct_reorder2(as.factor(annotID), clustname, clust),
    admix_groups = fct_reorder2(admix_groups, clustname, clust))

ords_admix_grps <- c("East", "North-East", "North-Feather","North-West", "South-West", "West")
admix_gg$admix_groups <- factor(admix_gg$admix_groups, levels = ords_admix_grps)
levels(admix_gg$admix_groups)

# replace V with Clust
admix_gg$clustname <- gsub("V", "Clust", admix_gg$clustname)

(plotk5 <- ggplot(data = admix_gg,
                 aes(x = annotID , y = clust, fill = clustname)) + 
   geom_col(position = "fill", width=1, show.legend = F) +
   scale_x_discrete(expand=c(0,0)) +
   scale_y_continuous(expand=c(0,0), labels = percent, breaks = seq(0,1,0.1)) +
   #scale_fill_viridis_d("Cluster") + 
   scale_fill_manual(values = c("Clust4"=cbbPalette[1], # East
                                "Clust5"=cbbPalette[2], # North-East
                                "Clust3"=cbbPalette[4], # North-Feather
                                "Clust1"=cbbPalette[3], #North-West
                                "Clust2"=cbbPalette[5])) + # West / South-West
   labs(title=paste0("NGSadmix k=",k)) +
   theme_classic(base_size = 8) +
   theme(axis.text.x = element_text(angle = 80, hjust = 1, size = 7)) +
   labs(x="", y="fraction ancestry") +
   facet_grid(.~ admix_groups, scales = "free_x") +
   theme(strip.background = element_blank(), 
         panel.spacing = unit(0, "lines"),
         panel.border = element_rect(fill = NA, color = "gray40"),
         axis.text.x=element_blank(),
         axis.ticks.x = element_blank()))

#ggsave(filename = paste0("figs/admix/",bamfile, "_k",k, "_admix.png"), width = 4.5, height = 2.7, scale = 1.3, units = "in", dpi = 300)
# width = 9, height = 6, units = "in", res = 300

# K6 ----------------------------------------------------------------------

k <- 6

# Read inferred admixture proportions
q<-read.table(paste0("data_output/admix/", bamfile, "_k", k,"_admix.qopt"))

# bind with metadat
admix <- bind_cols(q, annot)

admix_gg <- admix %>% gather(clustname, "clust", V1:V6) %>% 
  arrange(clustname, clust) %>% #group_by(clustname) %>% 
  mutate(
    clustname = fct_reorder2(factor(clustname, levels = c("V1", "V2", "V3", "V4", "V5", "V6")), clustname, clust),
    annotID = fct_reorder2(as.factor(annotID), clustname, clust),
    admix_groups = fct_reorder2(admix_groups, clustname, clust))

ords_admix_grps <- c("East", "North-East", "North-Feather","North-West", "South-West", "West")
admix_gg$admix_groups <- factor(admix_gg$admix_groups, levels = ords_admix_grps)
levels(admix_gg$admix_groups)

# replace V with Clust
admix_gg$clustname <- gsub("V", "Clust", admix_gg$clustname)

(plotk6 <- ggplot(data = admix_gg,
                  aes(x = annotID , y = clust, fill = clustname)) + 
    geom_col(position = "fill", width=1, show.legend = F) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0), labels = percent, breaks = seq(0,1,0.1)) +
    #scale_fill_viridis_d("Cluster") + 
    scale_fill_manual(values = c("Clust3"=cbbPalette[1], # East
                                 "Clust4"=cbbPalette[2], # North-East
                                 "Clust5"=cbbPalette[3], # North-West
                                 "Clust2"=cbbPalette[4], # North-Feather
                                 "Clust6"=cbbPalette[5], # West / South-West
                                 "Clust1"=cbbPalette[6]))+ # ??
    labs(title=paste0("NGSadmix k=",k)) +
    theme_classic(base_size = 8)+ #base_family = "Roboto") +
    theme(axis.text.x = element_text(angle = 80, hjust = 1, size = 7)) +
    labs(x="", y="fraction ancestry") +
    facet_grid(.~ admix_groups, scales = "free_x") +
    theme(strip.background = element_blank(), 
          panel.spacing = unit(0, "lines"),
          panel.border = element_rect(fill = NA, color = "gray40"),
          axis.text.x=element_blank(),
          axis.ticks.x = element_blank()))

#ggsave(filename = paste0("figs/admix/",bamfile, "_k",k, "_admix.png"), width = 4.5, height = 2.7, scale = 1.3, units = "in", dpi = 300)
# width = 9, height = 6, units = "in", res = 300

# K7 ----------------------------------------------------------------------

k <- 7

# Read inferred admixture proportions
q<-read.table(paste0("data_output/admix/", bamfile, "_k", k,"_admix.qopt"))

# bind with metadat
admix <- bind_cols(q, annot)

admix_gg <- admix %>% gather(clustname, "clust", V1:V7) %>% 
  arrange(clustname, clust) %>% #group_by(clustname) %>% 
  mutate(
    clustname = fct_reorder2(factor(clustname, levels = c("V1", "V2", "V3", "V4", "V5","V6", "V7")), clustname, clust),
    annotID = fct_reorder2(as.factor(annotID), clustname, clust),
    admix_groups = fct_reorder2(admix_groups, clustname, clust))

ords_admix_grps <- c("East", "North-East", "North-Feather","North-West", "South-West", "West")
admix_gg$admix_groups <- factor(admix_gg$admix_groups, levels = ords_admix_grps)
levels(admix_gg$admix_groups)

# replace V with Clust
admix_gg$clustname <- gsub("V", "Clust", admix_gg$clustname)

(plotk7 <- ggplot(data = admix_gg,
                  aes(x = annotID , y = clust, fill = clustname)) + 
    geom_col(position = "fill", width=1, show.legend = F) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0), labels = percent, breaks = seq(0,1,0.1)) +
     scale_fill_manual(values = c("Clust6"=cbbPalette[1], # East
                                  "Clust1"=cbbPalette[2], # North-East
                                  "Clust2"=cbbPalette[3], # North-West
                                  "Clust4"=cbbPalette[4], # North-Feather
                                  "Clust5"=cbbPalette[5], # West / South-West
                                  "Clust3"=cbbPalette[6], # branch of NorthWest
                                  "Clust7"=cbbPalette[7]))+ # Branch of Northwest
    labs(title=paste0("NGSadmix k=",k)) +
    theme_classic(base_size = 8) +
    #theme(axis.text.x = element_text(angle = 80, hjust = 1, size = 7)) +
    labs(x="", y="fraction ancestry") +
    facet_grid(.~ admix_groups, scales = "free_x") +
    theme(strip.background = element_blank(), 
          panel.spacing = unit(0, "lines"),
          panel.border = element_rect(fill = NA, color = "gray40"),
          axis.text.x=element_blank(),
          axis.ticks.x = element_blank()))

#ggsave(filename = paste0("figs/admix/",bamfile, "_k",k, "_admix.png"), width = 4.5, height = 2.7, scale = 1.3, units = "in", dpi = 300)
# width = 9, height = 6, units = "in", res = 300

# K8 ----------------------------------------------------------------------

k <- 8

# Read inferred admixture proportions
q<-read.table(paste0("data_output/admix/", bamfile, "_k", k,"_admix.qopt"))

# bind with metadat
admix <- bind_cols(q, annot)

admix_gg <- admix %>% gather(clustname, "clust", V1:V8) %>% 
  arrange(clustname, clust) %>% #group_by(clustname) %>% 
  mutate(
    clustname = fct_reorder2(factor(clustname, levels = c("V1", "V2", "V3", "V4", "V5","V6","V7","V8")), clustname, clust),
    annotID = fct_reorder2(as.factor(annotID), clustname, clust),
    admix_groups = fct_reorder2(admix_groups, clustname, clust))

ords_admix_grps <- c("East", "North-East", "North-Feather","North-West", "South-West", "West")
admix_gg$admix_groups <- factor(admix_gg$admix_groups, levels = ords_admix_grps)
levels(admix_gg$admix_groups)

# replace V with Clust
admix_gg$clustname <- gsub("V", "Clust", admix_gg$clustname)

# plot
(plotk8 <- ggplot(data = admix_gg,
                  aes(x = annotID , y = clust, fill = clustname)) + 
    geom_col(position = "fill", width=1, show.legend = F) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0), labels = percent, breaks = seq(0,1,0.1)) +
    scale_fill_manual(values = c("Clust7"=cbbPalette[1], # East
                                 "Clust2"=cbbPalette[2], # North-East
                                 "Clust3"=cbbPalette[3], # North-West
                                 "Clust1"=cbbPalette[4], # North-Feather
                                 "Clust5"=cbbPalette[5], # West / South-West
                                 "Clust8"=cbbPalette[6], # branch of NorthWest
                                 "Clust4"=cbbPalette[7],
                                 "Clust6"=cbbPalette[8]))+ # Branch of Northwest
    
    labs(title=paste0("NGSadmix k=",k)) +
    theme_classic(base_size = 8) +
    #theme(axis.text.x = element_text(angle = 80, hjust = 1, size = 7)) +
    labs(x="", y="fraction ancestry") +
    facet_grid(.~ admix_groups, scales = "free_x") +
    theme(strip.background = element_blank(), 
          panel.spacing = unit(0, "lines"),
          panel.border = element_rect(fill = NA, color = "gray40"),
          axis.text.x=element_blank(),
          axis.ticks.x = element_blank()))

#ggsave(filename = paste0("figs/admix/",bamfile, "_k",k, "_admix.png"), width = 4.5, height = 2.7, scale = 1.3, units = "in", dpi = 300)
# width = 9, height = 6, units = "in", res = 300

# K9 ----------------------------------------------------------------------

k <- 9

# Read inferred admixture proportions
q<-read.table(paste0("data_output/admix/", bamfile, "_k", k,"_admix.qopt"))

# bind with metadat
admix <- bind_cols(q, annot)

admix_gg <- admix %>% gather(clustname, "clust", V1:V9) %>% 
  arrange(clustname, clust) %>% #group_by(clustname) %>% 
  mutate(
    clustname = fct_reorder2(factor(clustname, levels = c("V1", "V2", "V3", "V4", "V5","V6","V7","V8","V9")), clustname, clust),
    annotID = fct_reorder2(as.factor(annotID), clustname, clust),
    admix_groups = fct_reorder2(admix_groups, clustname, clust))

ords_admix_grps <- c("East", "North-East", "North-Feather","North-West", "South-West", "West")
admix_gg$admix_groups <- factor(admix_gg$admix_groups, levels = ords_admix_grps)
levels(admix_gg$admix_groups)

# replace V with Clust
admix_gg$clustname <- gsub("V", "Clust", admix_gg$clustname)

# plot
(plotk9 <- ggplot(data = admix_gg,
                  aes(x = annotID , y = clust, fill = clustname)) + 
    geom_col(position = "fill", width=1, show.legend = F) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0), labels = percent, breaks = seq(0,1,0.1)) +
    scale_fill_manual(values = c("Clust8"=cbbPalette[1], # East
                                 "Clust3"=cbbPalette[2], # North-East
                                 "Clust9"=cbbPalette[3], # North-West
                                 "Clust7"=cbbPalette[4], # North-Feather
                                 "Clust2"=cbbPalette[5], # West / South-West
                                 "Clust5"=cbbPalette[6], # branch of NorthWest
                                 "Clust6"=cbbPalette[7],
                                 "Clust1"=cbbPalette[8],
                                 "Clust4"=cbbPalette[9]))+ # Branch of Northwest
    
    labs(title=paste0("NGSadmix k=",k)) +
    theme_classic(base_size = 8) +
    #theme(axis.text.x = element_text(angle = 80, hjust = 1, size = 7)) +
    labs(x="", y="fraction ancestry") +
    facet_grid(.~ admix_groups, scales = "free_x") +
    theme(strip.background = element_blank(), 
          panel.spacing = unit(0, "lines"),
          panel.border = element_rect(fill = NA, color = "gray40"),
          axis.text.x=element_blank(),
          axis.ticks.x = element_blank()))

#ggsave(filename = paste0("figs/admix/",bamfile, "_k",k, "_admix.png"), width = 4.5, height = 2.7, scale = 1.3, units = "in", dpi = 300)
# width = 9, height = 6, units = "in", res = 300


# SAVE PLOTS cowplot -----------------------------------------------------------

library(cowplot)

(kplots <- plot_grid(plotk3, plotk4, plotk5, plotk6, align = "hv", nrow = 2, labels = "AUTO"))

save_plot(plot = kplots , filename = paste0("figs/admix/admix_", bamfile, "_k3-6",".png"), base_width = 8, base_height = 4, base_aspect_ratio = 1.5, dpi=300)

(kplots <- plot_grid(plotk7, plotk8, plotk9, nrow=3, labels = "AUTO"))

save_plot(plot = kplots , filename = paste0("figs/admix/admix_", bamfile, "_k7-9",".png"), base_width = 8, base_height = 4, base_aspect_ratio = 1.5, dpi=300)

# ggsave(filename = "figs/figure_03_admix_rana_by_watershed_25k_k2.pdf", width = 7, height = 4, scale = 1.1, units = "in", dpi = 300)
# 
# ggsave(filename = "figs/figure_03_admix_rana_by_watershed_25k_k2_scaled.pdf", width = 4.5, height = 2.7, scale=1.3, units = "in", dpi = 300)
# 
# ggsave(filename = "figs/figure_03_admix_rana_by_watershed_25k_k2_scaled.png", width = 4.5, height = 2.7, scale = 1.3, units = "in", dpi = 300)


# PLOT LOGLIKS ------------------------------------------------------------

k <- c(2:12)
logs <- c(-33148006.008206, -31793778.761774, -30793224.392646, 
          -30546187.662216, -30316365.727090, -30121297.356706, -29978921.546022,
          -29737643.919022, -29825944.804667, -29498362.645870, -29431334.806164)
#nsites 69425 # nind 635
admixLogs <- tibble(k=k, bestLik=logs)

ggplot() + geom_point(data=admixLogs, aes(x=as.factor(k), y=bestLik), col="maroon", size=3) +
  theme_bw() + 
  labs(x="Number of groups (k)", y="Likelihood Scores")

ggsave(filename = "figs/admix_likelihood_scores.png", width = 4, height=2.7, units="in", dpi=300)


# OLDER EXPERIMENTS THAT BASICALLY FAILED ---------------------------------

# everything below here is old code I didn't like (or failed writing).

# PLOT Ks ----------------------------------------------------------------

# Plot them (ordered by population)
#png(filename = paste0("figs/admix/",bamfile, "_k",k, "_admix.png"), width = 9, height = 6, units = "in", res = 300)

#setup
# par(mar=c(10,4.5,2,1))
# 
# # test to see alignment of groups:
# #tst <- barplot(t(q)[,ord],col=viridis(k),names=annot$EcoRegion[ord], space=0, las=2,ylab="Admixture proportions",cex.names=0.6, border=NA, main= paste0("Admixture k=", k, " (", bamfile,")"))
# 
# # now replot for saving
# barplot(t(q)[,ord],col=viridis(k),names=annot$admix_groups[ord], space=0, las=2,ylab="Admixture proportions",cex.names=0.6, border=NA, main= paste0("Admixture k=", k, " (", bamfile,")"), xaxt='n')
# 
# #unique(annot$EcoRegion[ord])
# text(x = c(2, 40, 120, 185, 250, 440, 600, 630), par("usr")[3] - 0.2, labels = unique(annot$admix_groups[ord]) , srt = 90, pos = 1, xpd = TRUE, cex = 0.8) #col = "blue")
# 
# # ecoregion breaks
# abline(v = 7, col="black", lty=2)
# abline(v = 71, col="black", lty=2)
# abline(v = 171, col="black", lty=2)
# abline(v = 204, col="black", lty=2)
# abline(v = 296, col="black", lty=2)
# abline(v = 582, col="black", lty=2)
# abline(v = 582, col="black", lty=2)
# abline(v = 624, col="black", lty=2)
# 
# dev.off()




# POPHELPER ---------------------------------------------------------------

library(pophelper)

admixfiles <- list.files(path = "data_output/admix/", bamfile, full.names = T)
#admixfile <- paste0("data_output/admix/", qopt_file, "k", k,"_admix", admixk, ".qopt")

# read in multiple files
poplist <- readQ(admixfiles)

# tabulate
tabulateQ(qlist=poplist, sorttable = F)

# summarize
summariseQ(tabulateQ(qlist=poplist))

# GET LABELS:
admix_labels <- annot %>% select("admix_groups") %>% 
  mutate_at(c("admix_groups"), as.character)

# replace NA's if they exist:
if(length(admix_labels[is.na(admix_labels)])){
  admix_labels <- admix_labels %>% replace(., is.na(.), "unknown")
} else {
  print("no NAs")
}


#par(mar=c(11,3,2,2))
p1 <- plotQ(poplist[c(6:7)], 
            imgoutput="join", 
            returnplot = T, exportplot = T, quiet=T,
            #width = 7, height=4, units = "in", 
            font = "Roboto Condensed", 
            sharedindlab = F, 
            sortind = 'all', 
            splab = paste0("K=",sapply(poplist[c(6:7)],ncol)),splabsize = 5, 
            #subsetgrp=c("North Coast","Sierra Nevada", "Central Coast"), 
            selgrp="admix_groups",  
            showlegend=T, showtitle=T,showsubtitle=F, titlelab="All 100k Admixture",
            grplabangle = 0, grplabsize = 1, grplab=admix_labels, ordergrp=T, 
            grplabpos = 1, grplabheight = 3, grplabjust = 1, grplabspacer = .2,
            linepos = 1.1,
            #showindlab=T, indlabsize=2, indlabheight=0.1, indlabspacer=-1,indlabangle=70,
            barbordercolour="white", barbordersize=0, 
            outputfilename=paste0("figs/admix/",bamfile,"_k8-k10"), imgtype="png")

dev.off()
#print(p1$plot[[1]])




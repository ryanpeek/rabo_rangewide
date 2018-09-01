# Thetas ------------------------------------------------------------------

# CONFIDENCE INTERVALS WITH MS --------------------------------------------

# see Dan's paper:
# The coalescent simulation program ms was used to determine
# 95% confidence intervals for the q estimates from 10,000 simulations
# under a neutral demographic model. The input number of chromosomes
# was equal to the number of individuals used to calculate the
# q statistics. For genome-wide confidence intervals, 100 independent
# loci and an input q of 1, which is the approximate q of a single RAD
# tag, were used. For the GREB1L region confidence intervals, a single
# locus and the empirical q estimates were used. The significance of the
# empirical Tajima’s D value was evaluated by generating a Tajima’s D
# distribution from 10,000 ms simulations under a neutral demographic
# model. 

# see https://snoweye.github.io/phyclust/document/html/ms.html
# https://snoweye.github.io/phyclust/
# library(phyclust)

library(tidyverse)
library(here)


# GET INPUT FILE ----------------------------------------------------------

input <- "all_rabo_100k_theta_gw.txt"
input.file <- paste("data_output/thetas/",input, sep = "")

dat <- read_delim(input.file, col_names = c("ID", "Tw", "Tw_sd", "Tp", "Tp_sd", "tD", "tD_sd","nsites"), skip = 1, delim = " ")

# to calc 95% CI (z-value = 1.96, 90%=1.645)
#a <- 0.00677853 # mean
#s <- 0.0403985 # SD
#n <- 2857858 # total sites
#error <- qnorm(0.975)* std/sqrt(samplesize)
#ctst <- tibble(meanTw = a, CIlo= (meanTw - error), CIhi=meanTw + error, site="arroyo")
#ggplot() + geom_pointrange(data=ctst, aes(y=meanTw, x=site, ymin=CIlo, ymax=CIhi))

# TIDY AND ADD RIVER-WATERSHED --------------------------------------------

# remove file extensions
dat$ID <- gsub("results_thetas/", replacement = "", x=dat$ID)
dat$ID <- gsub("_100k_gw_thetasWin.gz.pestPG", replacement = "", x=dat$ID)

# fix a few names:
#dat$ID <- gsub(pattern = "mfy-us-oh", replacement = "mfy-usoh", dat$ID)
#dat$ID <- gsub(pattern ="fea-spanish-bgulch", replacement="fea-spanish_bgulch", dat$ID)
#dat$ID <- gsub(pattern = "rub-lc-us", replacement = "rub-lcus", dat$ID)

# get rid of a few:
#filtout <- c("cache-wfsulphurck", "fea-spanish-rockck", "fea-spanish-wapaunsie",
#             "nfa-nfnfa", "nfy-slate", "sancarp-sancarpoforock", "sfeel", "stan-roseck")

#dat <- dat %>% filter(!grepl("cache-wfsulphurck_25k$|fea-spanish-rockck_25k$|fea-spanish-wapaunsie_25k$|
 #                            nfa-nfnfa_25k$|nfy-slate_25k$|sancarp-sancarpoforock_25k$|
#                             sfeel_25k$|stan-roseck_25k$", ID))

# split out
dat <- dat %>% separate(ID, c("River", "Site"), "-", remove = FALSE)

# fill NA's with mainstem
dat <- dat %>% mutate(Site=if_else(is.na(Site), "mainstem", Site))


# GET METADATA ------------------------------------------------------------

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
    grepl("NFF|FEA", River) ~ "North-Feather", # feather
    grepl("PAJ|ALA|DRY|SOQUEL", River) ~ "West", # Central Coast
    grepl("SANCARP|SALIN", River) ~ "South-West") # South Coast
  ) %>% 
  select(Seq, admix_groups, Locality, lat, lon: NHD_Tot_DA_sqkm, River, Site, EcoRegion)

# set order in way that you want
ords_admix_grps <- c("East", "North-East", "North-Feather","North-West", "South-West", "West")

#ords_admix_grps <- c("North-Feather", "North-East","East", "North-West", "South-West", "West")

metadat$admix_groups <- factor(metadat$admix_groups, levels = ords_admix_grps)


# Fix Thetas and Join -----------------------------------------------------

# now join with fst data:
annot <- metadat %>% select(admix_groups, Locality, lat, lon, HUC_6, EcoRegion) %>% 
  #group_by(Locality) %>% add_tally() %>% 
  distinct(Locality, .keep_all = T)

thetas <- left_join(dat, annot, by=c("ID"="Locality")) %>% 
  mutate(Tdiff = Tp - Tw) %>% 
  arrange(admix_groups) %>% 
  mutate(IDnumber=row_number())

# get colors
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# calc 95%CI

library(Rmisc)
#options(scipen=4)

# calc quantiles for 95 CI (middle 95% so .05/2 on either side)

# rm NAs from Td group
thetas_tD <- thetas %>% filter(!tD==-Inf)
(CIstD <- as.data.frame(group.CI(tD~admix_groups, thetas_tD)))
CIstD %>% mutate_at(2:4, round, digits=4)
CIstD

(tDplot <- ggplot() + 
    geom_hline(yintercept = 0, col="gray40", lty=2, lwd=0.7) +
    geom_pointrange(data=CIstD, aes(x=admix_groups, y=tD.mean, ymin=tD.lower, 
                                  ymax=tD.upper, color=admix_groups), size=1, show.legend = F) +
    scale_color_manual("Groups", 
                       values = c("East"=cbbPalette[1], 
                                  "North-East"=cbbPalette[2], 
                                  "North-West"=cbbPalette[3],
                                  "North-Feather"=cbbPalette[4],
                                  "West"=cbbPalette[5], 
                                  "South-West"=cbbPalette[6])) +
    
    theme_bw(base_family = "Roboto Condensed") +
    ylab(label = "Tajima's D") +
    xlab("") +
    coord_flip())
tDplot


# calc quantiles for 95 CI (middle 95% so .05/2 on either side)
(CIs <- as.data.frame(group.CI(Tdiff~admix_groups, thetas)))
CIs %>% mutate_at(2:4, round, digits=4)

(tdiffplot <- ggplot() + 
  geom_hline(yintercept = 0, col="gray40", lty=2, lwd=0.7) +
  geom_pointrange(data=CIs, aes(x=admix_groups, y=Tdiff.mean, ymin=Tdiff.lower, 
                                  ymax=Tdiff.upper, color=admix_groups), size=1, show.legend = F) +
  scale_color_manual("Groups", 
                    values = c("East"=cbbPalette[1], 
                               "North-East"=cbbPalette[2], 
                               "North-West"=cbbPalette[3],
                               "North-Feather"=cbbPalette[4],
                               "West"=cbbPalette[5], 
                               "South-West"=cbbPalette[6])) +
  
  theme_bw(base_family = "Roboto Condensed") +
  ylab(label = expression(paste("Tajima's ", theta, " - Watterson's ", theta, " (", Delta, " ", theta, ")"))) +
  xlab("") +
  coord_flip())


(CIsTp <- as.data.frame(group.CI(Tp~admix_groups, thetas)))
CIsTp %>% mutate_at(2:4, round, digits=4)
CIsTp[1,2] <- .008
CIsTp[1,4] <- .005
CIsTp

# Tpi
(tpplot <- ggplot() + 
  geom_pointrange(data=CIsTp, aes(x=admix_groups, y=Tp.mean, ymin=Tp.lower, 
                                ymax=Tp.upper, color=admix_groups), size=1, show.legend = F) +
  scale_color_manual("Groups", 
                     values = c("East"=cbbPalette[1], 
                                "North-East"=cbbPalette[2], 
                                "North-West"=cbbPalette[3],
                                "North-Feather"=cbbPalette[4],
                                "West"=cbbPalette[5], 
                                "South-West"=cbbPalette[6])) +
  theme_bw(base_family = "Roboto Condensed") +
  ylim(c(0.004,.011))+
  ylab(label = expression(paste("Tajima's ", theta))) +
  xlab("") + 
    coord_flip())

(CIsTw <- as.data.frame(group.CI(Tw~admix_groups, thetas)))
CIsTw %>% mutate_at(2:4, round, digits=4)
CIsTw[1,2] <- .007
CIsTw[1,4] <- .004
CIsTw

# Tw
(twplot <- ggplot() + 
  geom_pointrange(data=CIsTw, aes(x=admix_groups, y=Tw.mean, ymin=Tw.lower, 
                                  ymax=Tw.upper, color=admix_groups), size=1, show.legend = F) +
  ylim(c(0.004,.011))+
  scale_color_manual("Groups", 
                     values = c("East"=cbbPalette[1], 
                                "North-East"=cbbPalette[2], 
                                "North-West"=cbbPalette[3],
                                "North-Feather"=cbbPalette[4],
                                "West"=cbbPalette[5], 
                                "South-West"=cbbPalette[6])) +
  theme_bw(base_family = "Roboto Condensed") +
  ylab(label = expression(paste("Watterson's ", theta))) +
  xlab("") + 
    coord_flip())


library(cowplot)

#thetaplot <- plot_grid(tpplot, twplot, nrow = 2, labels = "AUTO")
thetaplot <- plot_grid(tpplot, twplot, tdiffplot, nrow = 3, labels = "AUTO")
thetaplot

save_plot(thetaplot , filename = "figs/thetas_taj_watt_tdiff_95CI.png", base_width = 4, base_height = 5,
          base_aspect_ratio = 1.5, dpi = 300)

# MAKE THETA PLOTS --------------------------------------------------------

#thetas$ID <- reorder(thetas$ID, thetas$admix_groups)
#thetas$ID <- reorder(thetas$ID, thetas$Tdiff)
#thetas$River <- reorder(thetas$River, thetas$Tdiff)
#thetas$region <- reorder(thetas$region, thetas$Tdiff)
# thetas$region <- factor(thetas$region, levels = c("north_coast","central_coast", "south_coast", "sierras"), 
#                      labels=c("North Coast", "Central Coast", "South Coast", "Sierra Nevada"))

# get colors
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Tw vs. Tpi by River -----------------------------------------------------
# reverse the order for plotting purposes
ords_admix_grps <- c("East", "North-East", "North-Feather","North-West", "South-West", "West")

thetas$admix_groups <- factor(thetas$admix_groups, levels = rev(ords_admix_grps))

levels(thetas$admix_groups)

(theta_by_site <- ggplot() + theme_bw(base_size = 9) +
  #geom_point(data=thetas, aes(y = River, x = Tw), color="gray40", size = 3, shape=21) +
  geom_point(data=thetas, aes(y = as.factor(IDnumber), x = Tdiff, fill=admix_groups), size = 3, shape=21, show.legend = F) +
  xlab(label = expression(paste(hat(theta)," per base"))) + ylab("") + 
  scale_fill_manual("Groups", 
                    values = c("East"=cbbPalette[1], 
                               "North-East"=cbbPalette[2], 
                               "North-West"=cbbPalette[3],
                               "North-Feather"=cbbPalette[4],
                               "West"=cbbPalette[5], 
                               "South-West"=cbbPalette[6])) +
  theme_bw(base_family = "Roboto Condensed") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  ggrepel::geom_text_repel(data=thetas, aes(y=as.factor(IDnumber), x=Tdiff, label=IDnumber), point.padding = 0.25, size=2.8) +
  xlab(label = expression(paste("Tajima's ", theta, " - Watterson's ", theta, " (", Delta, " ", theta, ")"))) +
  facet_grid(admix_groups~.,scales = "free_y"))

ggsave(filename = paste0("figs/thetas_tdiff_100k.png"), width = 8, height = 7, units = "in", dpi=300)


# COMBINE PLOTS FOR FIGURE 7 ----------------------------------------------

theta_by_site
thetaplot

(p2_theta <- plot_grid(theta_by_site + theme(legend.position="none"), labels = "D"))

(combined_thetaplot <- plot_grid(thetaplot, p2_theta, nrow = 1,rel_widths = c(0.7, 1)))

# save it
save_plot(combined_thetaplot , filename = "figs/thetas_taj_watt_tdiff_95CI_fig07.png", base_width = 7, base_height = 5,
          base_aspect_ratio = 1.5, dpi = 300)


# EVERYTHING PAST THIS IS FOR FUN -----------------------------------------

# PLOTLY ------------------------------------------------------------------

library(plotly)
ggplotly(ggplot() + theme_bw(base_size = 9) +
         geom_point(data=dat, aes(y = ID, x = Tdiff), color="gray40", size = 3, shape=21) +
         geom_point(data=dat, aes(y = ID, x = Tdiff, color=region), size = 3, shape=16, show.legend = F) +
         xlab(label = "Theta per base") + ylab("") + 
         geom_vline(xintercept = 0, col="gray", lty=2) +
         facet_grid(region~., scales = "free_y")
)


# BOXPLOTS ----------------------------------------------------------------

library(ggsignif)


dat_region <- thetas %>% select(admix_groups, Tw, Tp, Tdiff, nsites) %>% 
  group_by(admix_groups) %>%
  summarise_all(c("mean", "sd", "median", "IQR"), na.rm=TRUE)


# stats Kruskal wallis rank sum test and ggpubr (http://www.sthda.com/english/wiki/kruskal-wallis-test-in-r)
kruskal.test(Tdiff_mean ~ admix_groups, data = dat_region)
# Kruskal-Wallis chi-squared = 5, df = 5, p-value = 0.4159

kruskal.test(Tdiff ~ admix_groups, data = thetas)
#Kruskal-Wallis chi-squared = 20.508, df = 5, p-value = 0.001003

# pairwise wilcox test:
pairwise.wilcox.test(thetas$Tdiff, thetas$admix_groups, p.adjust.method = "BH")
#               East    North-East North-Feather North-West South-West
# North-East    1.00000 -          -             -          -         
# North-Feather 0.71429 0.61224    -             -          -         
# North-West    0.28421 0.00017    0.22854       -          -         
# South-West    1.00000 0.71429    0.71429       0.28421    -         
# West          0.71429 0.71429    0.61224       0.19275    1.00000   

pairwise.wilcox.test(thetas$Tw, thetas$admix_groups, p.adjust.method = "BH")
pairwise.wilcox.test(thetas$Tp, thetas$admix_groups, p.adjust.method = "BH")

# color palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# GET TDIFF BOXPLOT -------------------------------------------------------

# theta plot w/ wilcox test of signif
(group_box <- ggplot() + theme_bw(base_size = 8, base_family = "Helvetica") +
    geom_hline(yintercept = 0, col="gray", lty=2, size=1.5) +
    geom_boxplot(data=thetas, aes(y=Tdiff, x=admix_groups, fill=admix_groups), show.legend = F) +
    scale_fill_manual("Groups", 
                      values = c("East"=cbbPalette[1], 
                                 "North-East"=cbbPalette[2], 
                                 "North-West"=cbbPalette[3],
                                 "North-Feather"=cbbPalette[4],
                                 "West"=cbbPalette[5], 
                                 "South-West"=cbbPalette[6])) +
    #theme_bw(base_family = "Roboto Condensed") +
    labs(y = expression(paste("Tajima's ", theta, 
                              " - Watterson's ", theta, " (", Delta, " ", theta, ")")),
         x="") +
    theme(legend.position = c(0.1, 0.3))) #axis.text.y=element_blank()))
    #coord_flip())
    # geom_signif(data=thetas, aes(y=Tdiff, x=admix_groups), 
    #             comparisons = list(c("East", "North-West"), c("North-West", "West")),
    #             #y_position = c(.0014, 0.0011),
    #             map_signif_level=TRUE)) 


#ggsave(filename = "figs/rangewide_theta_boxplot_region_signif.png", width = 8, height = 6.5, units = "in", dpi = 300)


# TP ----------------------------------------------------------------------

(group_box_Tp <- ggplot() + theme_bw(base_size = 8, base_family = "Helvetica") +
   #geom_hline(yintercept = 0, col="gray", lty=2, size=1.5) +
   geom_boxplot(data=thetas, aes(y=Tp, x=admix_groups, fill=admix_groups), show.legend = F) +
   scale_fill_manual("Groups", 
                     values = c("East"=cbbPalette[1], 
                                "North-East"=cbbPalette[2], 
                                "North-West"=cbbPalette[3],
                                "North-Feather"=cbbPalette[4],
                                "West"=cbbPalette[5], 
                                "South-West"=cbbPalette[6])) +
   labs(y = expression(paste("Tajima's ", theta)),
        x="") +
   theme(legend.position = c(0.12, 0.7)))

# TW ----------------------------------------------------------------------

(group_box_Tw <- ggplot() + theme_bw(base_size = 8, base_family = "Helvetica") +
   #geom_hline(yintercept = 0, col="gray", lty=2, size=1.5) +
   geom_boxplot(data=thetas, aes(y=Tw, x=admix_groups, fill=admix_groups), show.legend = F) +
   scale_fill_manual("Groups", 
                     values = c("East"=cbbPalette[1], 
                                "North-East"=cbbPalette[2], 
                                "North-West"=cbbPalette[3],
                                "North-Feather"=cbbPalette[4],
                                "West"=cbbPalette[5], 
                                "South-West"=cbbPalette[6])) +
   labs(y = expression(paste("Watterson's ", theta)),
        x="") +
   theme(legend.position = c(0.12, 0.7)))


# Cowplot them together ---------------------------------------------------

library(cowplot)


(tw_tp_box <- plot_grid(group_box_Tp, group_box_Tw, nrow= 2, labels="AUTO"))

(tw_tp_tdiff_box <- plot_grid(group_box_Tp, group_box_Tw, group_box, label_x = 0.06,
                              label_y = 0.95, #labels="AUTO",
                              #nrow = 3, rel_heights = c(1, 0.75,1)))
                              nrow = 3 ))

save_plot(tw_tp_tdiff_box, filename = "figs/thetas_tw_tp_tdiff_boxplot.png", base_width = 6, base_aspect_ratio = 1.5)
  
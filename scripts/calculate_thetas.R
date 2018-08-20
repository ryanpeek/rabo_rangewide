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


# GET INPUT FILE ----------------------------------------------------------

input <- "theta_counts_gw_25k_n3.txt"
input.file <- paste("data_output/thetas/",input, sep = "")

dat <- read_delim(input.file, col_names = c("ID", "Tw", "Tp", "tD","nsites"), skip = 1, delim = " ")

# TIDY AND ADD RIVER-WATERSHED --------------------------------------------

# remove file extensions
dat$ID <- gsub("results_thetas/", replacement = "", x=dat$ID)
dat$ID <- gsub("_gw_thetasWin.gz.pestPG", replacement = "", x=dat$ID)

# fix a few names:
dat$ID <- gsub(pattern = "mfy-us-oh", replacement = "mfy-usoh", dat$ID)
dat$ID <- gsub(pattern ="fea-spanish-bgulch", replacement="fea-spanish_bgulch", dat$ID)
dat$ID <- gsub(pattern = "rub-lc-us", replacement = "rub-lcus", dat$ID)

# get rid of a few:
#filtout <- c("cache-wfsulphurck", "fea-spanish-rockck", "fea-spanish-wapaunsie",
#             "nfa-nfnfa", "nfy-slate", "sancarp-sancarpoforock", "sfeel", "stan-roseck")

dat <- dat %>% filter(!grepl("cache-wfsulphurck_25k$|fea-spanish-rockck_25k$|fea-spanish-wapaunsie_25k$|
                             nfa-nfnfa_25k$|nfy-slate_25k$|sancarp-sancarpoforock_25k$|
                             sfeel_25k$|stan-roseck_25k$", ID))

# split out
dat <- dat %>% separate(ID, c("River", "Site"), "-", remove = FALSE)

# fill NA's with mainstem
dat <- dat %>% mutate(Site=if_else(is.na(Site), "mainstem", Site))

# add region
dat <- dat %>% 
  mutate(region = as.factor(case_when(
    grepl("eel|cache|lagun|klam|mad|mat|russ|sfeel|smith|ssantiam|sumpqua|trin|vandz", ID) ~ "north_coast",
    grepl("sfa|tuo|deer|nff|sfy|rub|nfy|mfy|nfmfa|nfa|mfa|fea|calav|bear", ID) ~ "sierras",
    grepl("ala|put|soquel", ID) ~ "central_coast",
    grepl("paj|salin|sancarp|dry", ID) ~ "south_coast"
    ))) %>% 
  select(ID:Site, region, Tw:nsites) %>% 
  mutate(Tdiff = Tp - Tw)


# MAKE THETA PLOTS --------------------------------------------------------

dat$ID <- reorder(dat$ID, dat$Tdiff)
dat$River <- reorder(dat$River, dat$Tdiff)
#dat$region <- reorder(dat$region, dat$Tdiff)
dat$region <- factor(dat$region, levels = c("north_coast","central_coast", "south_coast", "sierras"), 
                     labels=c("North Coast", "Central Coast", "South Coast", "Sierra Nevada"))

# Tw vs. Tpi by River -----------------------------------------------------

ggplot() + theme_bw(base_size = 9) +
  geom_point(data=dat, aes(y = River, x = Tw), color="gray40", size = 3, shape=21) +
  geom_point(data=dat, aes(y = River, x = Tp, color=region), size = 3, shape=16, show.legend = F) +
  xlab(label = expression(paste(hat(theta)," per base"))) + ylab("") + 
  labs(title=expression(paste("GW: ", theta[w]," = open circles,", theta[pi]," = filled"))) + 
  facet_grid(region~., scales = "free_y")

#ggsave(filename = paste0("figs/thetas_Tw_Tp_",subsample,".png"), width = 8, height = 7, units = "in", dpi=200)


# Tdiff by River ----------------------------------------------------------

ggplot() + theme_bw(base_size = 9) +
  geom_point(data=dat, aes(y = River, x = Tdiff), color="gray40", size = 3, shape=21) +
  geom_point(data=dat, aes(y = River, x = Tdiff, color=region), size = 3, shape=16, show.legend = F) +
  xlab(label = expression(paste(hat(theta)," per base"))) + ylab("") + 
  labs(title=expression(paste("GW: ", theta[w]," = open circles,", theta[pi]," = filled"))) + 
  geom_vline(xintercept = 0, col="gray", lty=2) +
  facet_grid(region~., scales = "free_y")


# Tw vs. Tpi by Locality -----------------------------------------------------

ggplot() + theme_bw(base_size = 9) +
  geom_point(data=dat, aes(y = ID, x = Tw), color="gray40", size = 3, shape=21) +
  geom_point(data=dat, aes(y = ID, x = Tp, color=region), size = 3, shape=16, show.legend = F) +
  xlab(label = expression(paste(hat(theta)," per base"))) + ylab("") + 
  labs(title=expression(paste("GW: ", theta[w]," = open circles,", theta[pi]," = filled"))) + 
  facet_grid(region~., scales = "free_y")

ggsave(filename = paste0("figs/thetas_Tp-Tw_rangewide_n5-9.png"), width = 11, height = 8, units = "in", dpi=300)

# Tdiff by Locality ----------------------------------------------------------

ggplot() + theme_bw(base_size = 9) +
  geom_point(data=dat, aes(y = ID, x = Tdiff), color="gray40", size = 3, shape=21) +
  geom_point(data=dat, aes(y = ID, x = Tdiff, color=region), size = 3, shape=16, show.legend = F) +
  xlab(label = expression(paste(hat(theta)," per base"))) + ylab("") + 
  labs(title=expression(paste("GW: ", theta[w]," = open circles,", theta[pi]," = filled"))) + 
  geom_vline(xintercept = 0, col="gray", lty=2) +
  facet_grid(region~., scales = "free_y")

ggsave(filename = paste0("figs/thetas_Tdiff_p-w_rangewide_n5-9.png"), width = 11, height = 8, units = "in", dpi=300)

# WOW

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


dat_region <- dat %>% select(region, Tw, Tp, Tdiff, nsites) %>% 
  group_by(region) %>%
  summarise_all(c("mean", "sd", "median", "IQR"), na.rm=TRUE)


# stats Kruskal wallis rank sum test and ggpubr (http://www.sthda.com/english/wiki/kruskal-wallis-test-in-r)
kruskal.test(Tdiff_mean ~ region, data = dat_region)
#Kruskal-Wallis chi-squared = 10.698, df = 2, p-value = 0.004754

# pairwise wilcox test:
pairwise.wilcox.test(dat$Tdiff, dat$region, p.adjust.method = "BH")
pairwise.wilcox.test(dat$Tdiff, dat$region, p.adjust.method = "bonferroni")

# color palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# theta plot w/ wilcox test of signif
(group_box <- ggplot() + theme_bw(base_size = 8, base_family = "Helvetica") +
    
    geom_hline(yintercept = 0, col="gray", lty=2, size=1.5) +
    #ylim(-.0009,0.00055)+
    geom_boxplot(data=dat, aes(y=Tdiff, x=region, fill=region), show.legend = T) +
    scale_fill_viridis_d()+
    # scale_fill_manual("Regulation Type", values = c("Unregulated"=cbbPalette[3], 
    #                                                 "Bypass"=cbbPalette[2], 
    #                                                 "Hydropeaking"=cbbPalette[7]),
    #                   guide=guide_legend(reverse=TRUE)) +
    labs(y = expression(paste("Tajima's ", theta, 
                              " - Watterson's ", theta, " (", Delta, " ", theta, ")")),
         x="") +
    #theme(legend.position = c(0.1, 0.75), axis.text.y=element_blank()) +
    geom_signif(data=dat, aes(y=Tdiff, x=region), comparisons = list(c("Sierra Nevada", "North Coast"), 
                                                                     c("Sierra Nevada", "Central Coast")), 
                y_position = c(.0014, 0.0011),
                map_signif_level=TRUE)) #+ coord_flip())


ggsave(filename = "figs/rangewide_theta_boxplot_region_signif.png", width = 8, height = 6.5, units = "in", dpi = 300)

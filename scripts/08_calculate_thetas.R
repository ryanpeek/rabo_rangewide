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

input <- "theta_counts_gw.txt"
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

# split out
dat <- dat %>% separate(ID, c("River", "Site"), "-", remove = FALSE)

# fill NA's with mainstem
dat <- dat %>% mutate(Site=if_else(is.na(Site), "mainstem", Site))

# add region
dat <- dat %>% 
  mutate(region = as.factor(case_when(
    grepl("eel|lagun|klam|mad|mat|russ|sfeel|smith|ssantiam|sumpqua|trin|vandz", ID) ~ "north_coast",
    grepl("sfy|rub|nfy|mfy|nfmfa|nfa|mfa|fea|calav|bear", ID) ~ "sierras",
    grepl("ala|put|soquel", ID) ~ "central_coast",
    grepl("paj|salin|sancarp|dry", ID) ~ "south_coast"
    ))) %>% 
  select(ID:Site, region, Tw:nsites) %>% 
  mutate(Tdiff = Tp - Tw)


# MAKE THETA PLOTS --------------------------------------------------------

dat$ID <- reorder(dat$ID, dat$Tdiff)
dat$River <- reorder(dat$River, dat$Tdiff)
#dat$region <- reorder(dat$region, dat$Tdiff)
dat$region <- factor(dat$region, levels = c("north_coast","central_coast", "south_coast", "sierras"))

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


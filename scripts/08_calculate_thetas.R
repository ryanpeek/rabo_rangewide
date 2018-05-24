# Thetas ------------------------------------------------------------------



# see Dan's paper:
# The coalescent simulation program ms (92) was used to determine
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
library(phyclust)


input <- "theta_counts_gw.txt"
input.file <- paste("data_output/thetas/",input, sep = "")

#dat <- read_delim(input.file, col_names = c("ID", "Tw", "Tp", "tD","nsites"),skip = 1, delim = " ")
dat <- read_delim(input.file, delim = " ") %>% rename(ID=`#id`, Tw=tw, Tp=tp, nsites=ns)

dat$site<- gsub("(?>=results_thetas/).(?<=_gw_thetasWin.gz.pestPG)", "", x = dat$ID, perl=TRUE)
dat$site<- stringi::stri_extract(dat$ID, regex='[^_]*')

dat$site <- gsub("results_thetas/", replacement = "", x=dat$ID)
dat$site <- gsub("_gw_thetasWin.gz.pestPG", replacement = "", x=dat$site)

# MAKE THETA PLOTS --------------------------------------------------------

# Tw vs. Tpi (all sites)
ggplot() + theme_bw(base_size = 9) +
  geom_point(data=dat, aes(y = site, x = Tw), color="gray40", size = 3, shape=21) +
  geom_point(data=dat, aes(y = site, x = Tp, color=site), size = 3, shape=16, show.legend = F) +
  #scale_x_continuous(trans = 'log2', breaks=c(0.00025, 0.0005, 0.001,0.002,0.004,0.008, 0.016), labels=c("0.00025", "0.0005","0.001","0.002","0.004","0.008", "0.016")) +
  xlab(label = expression(paste(hat(theta)," per base"))) + ylab("") + 
  labs(title=expression(paste("GW: ", theta[w]," = open circles,", theta[pi]," = filled")))
#annotate("text", label="Tw = open circles \n Tp = filled", x = 0.010, y=38)

#ggsave(filename = paste0("figs/thetas_Tw_Tp_",subsample,".png"), width = 8, height = 7, units = "in", dpi=200)


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


input <- "theta_gw"
input.file <- paste("data_output/thetas/",input, sep = "")

#dat <- read_delim(input.file, col_names = c("ID", "Tw", "Tp", "tD","nsites"),skip = 1, delim = " ")
dat <- read_delim(input.file, delim = " ") %>% rename(ID=`#id`, Tw=tw, Tp=tp, nsites=ns)
dat$site<- stringi::stri_extract(dat$ID, regex='[^_]*')

# filter to single subsample (100k)
subsample <- "100k"
datF <- filter(dat, grepl(pattern = subsample, ID))
datF$ID <- gsub("_thetasWin.gz.pestPG", replacement = "", x = datF$ID) # for gw

# for gw
datF$ID <- gsub(paste0("_",subsample,"_gw"), replacement="", x=datF$ID)

# make a siteID 100k for joining
# df100.id <- tibble("ID"=datF$ID, "region"= NA, "reg"=NA)
# head(df100.id)
# df100.id$region[c(1,7,14:16,38:39,42:45,47:48,57:59,62:64)] <- "coastal"
# df100.id$region <- ifelse(is.na(df100.id$region),"sierra",df100.id$region)
# df100.id$reg[c(2:6,9:10,17:21,40:41,46,49:56,60,61)] <- "R" # in a reg mainstem watershed (but may be unreg trib)
# df100.id$reg <- ifelse(is.na(df100.id$reg),"U",df100.id$reg)
# save(df100.id, file = "data_output/theta_df100_ids.rda")
load("data_output/theta_df100_ids.rda")

# make a siteID 75k for joining
# df75.id <- tibble("ID"=datF$ID)
# head(df75.id)
# df75.id <- left_join(df75.id, df100.id, by="ID")
# df75.id$region[c(7,11)] <- "sierra"
# df75.id$region <- ifelse(is.na(df75.id$region),"coastal",df75.id$region)
# df75.id$reg[c(7,11)] <- "R" # in a reg mainstem watershed (but may be unreg trib)
# df75.id$reg <- ifelse(is.na(df75.id$reg),"U",df75.id$reg)
# save(df75.id, file = "data_output/theta_df75_ids.rda")
load("data_output/theta_df75_ids.rda")

# make a siteID 50k for joining
# df50.id <- tibble("ID"=datF$ID)
# head(df50.id)
# df50.id <- left_join(df50.id, df75.id, by="ID")
# df50.id$region <- ifelse(is.na(df50.id$region),"coastal",df50.id$region)
# df50.id$reg <- ifelse(is.na(df50.id$reg),"U",df50.id$reg)
# save(df50.id, file = "data_output/theta_df50_ids.rda")
load("data_output/theta_df50_ids.rda")

datF <- left_join(datF, df100.id, by="ID")
#datF <- left_join(datF, df75.id, by="ID")
#datF <- left_join(datF, df50.id, by="ID")

# MAKE THETA PLOTS --------------------------------------------------------

# drop lyon's peak (RASI)
datF <- filter(datF, ID!="nfa_lyonspk")

# Tw vs. Tpi (all sites)
ggplot() + theme_bw(base_size = 9) +
  geom_point(data=datF, aes(y = ID, x = Tw), color="gray40", size = 3, shape=21) +
  geom_point(data=datF, aes(y = ID, x = Tp, color=region), size = 3, shape=16, show.legend = T) +
  #scale_x_continuous(trans = 'log2', breaks=c(0.00025, 0.0005, 0.001,0.002,0.004,0.008, 0.016), labels=c("0.00025", "0.0005","0.001","0.002","0.004","0.008", "0.016")) +
  xlab(label = expression(paste(hat(theta)," per base"))) + ylab("") + 
  labs(title=expression(paste("GW: ", theta[w]," = open circles,", theta[pi]," = filled")), subtitle=subsample)
#annotate("text", label="Tw = open circles \n Tp = filled", x = 0.010, y=38)

#ggsave(filename = paste0("figs/thetas_Tw_Tp_",subsample,".png"), width = 8, height = 7, units = "in", dpi=200)

# Tw vs. Tpi (sierras only, col by reg)
datF_sierra <- filter(datF, region=="sierra")
ggplot() + theme_bw(base_size = 9) +
  geom_point(data=datF_sierra, aes(y = ID, x = Tw), color="gray40", size = 3, shape=21) +
  geom_point(data=datF_sierra, aes(y = ID, x = Tp, color=reg), size = 3, shape=16, show.legend = T) +
  xlab(label = expression(paste(hat(theta)," per base"))) + ylab("") + 
  labs(title=expression(paste("GW: ", theta[w]," = open circles,", theta[pi]," = filled")), subtitle=subsample)

#ggsave(filename = paste0("figs/thetas_Tw_Tp_sierras_",subsample, ".png"), width = 8, height = 7, units = "in", dpi=200)

# Tw vs. Tpi (coastal only, col by reg)
datF_coastal <- filter(datF, region=="coastal")
ggplot() + theme_bw(base_size = 9) +
  geom_point(data=datF_coastal, aes(y = ID, x = Tw), color="gray40", size = 3, shape=21) +
  geom_point(data=datF_coastal, aes(y = ID, x = Tp, color=reg), size = 3, shape=16, show.legend = T) +
  xlab(label = expression(paste(hat(theta)," per base"))) + ylab("") + 
  labs(title=expression(paste("GW: ", theta[w]," = open circles,", theta[pi]," = filled")), subtitle=subsample)

#ggsave(filename = paste0("figs/thetas_Tw_Tp_coastal_",subsample, ".png"), width = 8, height = 7, units = "in", dpi=200)


# SIERRA MEANS ------------------------------------------------------------

# summarize thetas
datF_sierra %>% group_by(reg) %>% 
  summarize(mean_Tw=mean(Tw),
            mean_Tp=mean(Tp)) %>% 
  ggplot(.) + geom_point(aes(x=reg, y = mean_Tw), color="gray40", size = 4, shape=21) +
  geom_point(aes(x=reg, y = mean_Tp, color=reg), size = 4, shape=16, show.legend = T) +
  xlab(label = expression(paste(hat(theta)," per base"))) + ylab("")


# TAJ D -------------------------------------------------------------------

datF_tD <- filter(datF, tD!=-Inf)

# all sites
ggplot() + theme_bw(base_size = 9) +
  geom_point(data=datF_tD, aes(y = ID, x = tD, color=site), size = 3, shape=16, show.legend = F) +
  xlab(label = "Tajima's D") + ylab("") +
  labs(title=paste0("Tajima's D at ", subsample))

ggsave(filename = paste0("figs/tD_",subsample,"_all.png"), width = 8, height = 7, units = "in", dpi=200)

# sierras sites
ggplot() + theme_bw(base_size = 9) +
  geom_point(data=datF_tD[datF_tD$region=="sierra",], aes(y = ID, x = tD, color=reg), size = 3, shape=16, show.legend = F) +
  xlab(label = "Tajima's D") + ylab("") +
  labs(title=paste0("Tajima's D at ", subsample))

ggsave(filename = paste0("figs/tD_",subsample, "_sierras.png"), width = 8, height = 7, units = "in", dpi=200)


# coastal sites
ggplot() + theme_bw(base_size = 9) +
  geom_point(data=datF_tD[datF_tD$region=="coastal",], aes(y = ID, x = tD, color=reg), size = 3, shape=16, show.legend = F) +
  xlab(label = "Tajima's D") + ylab("") +
  labs(title=paste0("Tajima's D at ", subsample))

ggsave(filename = paste0("figs/tD_",subsample, "_coastal.png"), width = 8, height = 7, units = "in", dpi=200)


# OTHER STATS -------------------------------------------------------------

# from here: http://rstudio-pubs-static.s3.amazonaws.com/150385_e616a4453cb5465f9736a7d5db2643a3.html


library(cowplot)
library(dplyr)
library(data.table)
library(tidyr)
library(rtracklayer)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# get data
teo<-fread("~/Desktop/ecl243_teo.thetas.graph.me") %>% mutate(taxa="teo") %>% select(3:15)
mz<-fread("~/Desktop/ecl243.thetas.graph.me") %>% mutate(taxa="mz") %>% select(3:15)
dat<-rbind(teo,mz) %>% mutate(pi=tP/nSites,Mb=WinCenter/1E6) 

# filter and plot
filter(dat,taxa=="mz") %>% 
  ggplot(aes(y=pi,x=Mb))+
  geom_point(color=cbPalette[1])+
  geom_smooth(method="loess",span=0.1,color=cbPalette[1])+
  ggtitle("maize nucleotide diversity")+
  xlab("Mb") + 
  ylab(expression(hat(theta[pi])))+
  theme(axis.text=element_text(size=14), plot.title=element_text(size=18),axis.title=element_text(size=18,face="bold")) 

filter(dat,taxa=="mz",nSites>1000) %>% 
  ggplot(aes(y=pi,x=Mb))+
  geom_point(color=cbPalette[1])+
  geom_smooth(method="loess",span=0.1,color=cbPalette[1])+
  ggtitle("maize nucleotide diversity ( > 1Kb sequence per window)")+
  xlab("Mb") + 
  ylab(expression(hat(theta[pi])))+
  theme(axis.text=element_text(size=14), plot.title=element_text(size=18),axis.title=element_text(size=18,face="bold")) 

filter(dat,taxa=="mz",nSites>1000) %>% 
  ggplot(aes(y=Tajima,x=Mb))+
  geom_point(color=cbPalette[3])+
  geom_smooth(method="loess",span=0.1,color=cbPalette[3])+
  ggtitle("maize Tajima's D")+
  xlab("Mb") + 
  ylab("Tajima D")+
  theme(axis.text=element_text(size=14), plot.title=element_text(size=18),axis.title=element_text(size=18,face="bold")) 


#pi vs watterson
filter(dat,nSites>1000,taxa=="mz") %>% 
  mutate(watterson=tW/nSites) %>%
  ggplot(aes(y=watterson,x=pi))+
  geom_point()+
  ggtitle("maize nucleotide diversity")+
  ylab(expression(hat(theta[W])))+scale_y_log10()+
  xlab(expression(hat(theta[pi])))+scale_x_log10()+
  geom_smooth(method="lm")+
  theme(axis.text=element_text(size=14), plot.title=element_text(size=18),axis.title=element_text(size=18,face="bold"),legend.position=c(0.8, 0.9)) 

#do some machinations, calculate ratio of pi values
merge(mz,teo,by="WinCenter",suffixes=c(".mz",".teo")) %>% 
  mutate(rat=(tP.mz/nSites.mz)/(tP.teo/nSites.teo),Mb=WinCenter/1E6) %>%
  filter(nSites.mz>1000,nSites.teo>1000) %>% 
  ggplot(aes(y=rat,x=Mb))+
  geom_point()+
  ggtitle("here be sweeps?")+
  ylab(expression(paste(theta[pi],"maize / ",theta[pi],"teosinte")))+
  xlab("Mb")+
  geom_smooth(method="loess",span=0.05)+
  theme(axis.text=element_text(size=14), plot.title=element_text(size=18),axis.title=element_text(size=18,face="bold"),legend.position=c(0.8, 0.9))



# GET GFF

chr10_genes<-import.gff(con="~/Desktop/Zea_mays.AGPv3.30.chromosome.10.gff3",format="gff3") %>% 
  as("data.frame") %>% 
  mutate(start=start/1E6,end=end/1E6) %>%
  filter(start<25,start>15,type=="gene")
chr10_exons<-import.gff(con="~/Desktop/Zea_mays.AGPv3.30.chromosome.10.gff3",format="gff3") %>% 
  as("data.frame") %>% 
  mutate(start=start/1E6,end=end/1E6) %>%
  filter(start<25,start>15,type=="exon")

# Plot with genes, zoomed in to our region

merge(mz,teo,by="WinCenter",suffixes=c(".mz",".teo")) %>% 
  mutate(rat=(tP.mz/nSites.mz)/(tP.teo/nSites.teo),Mb=WinCenter/1E6) %>%
  filter(nSites.mz>1000,nSites.teo>1000,Mb>17.5,Mb<18) %>% 
  ggplot(aes(y=rat,x=Mb))+
  geom_point()+
  ggtitle("here be sweeps?")+
  ylab(expression(paste(theta[pi],"maize / ",theta[pi],"teosinte")))+
  xlab("Mb")+
  scale_y_log10()+xlim(c(17.5,18))+
  geom_segment(aes(x = start, y = 0.08, xend = end, yend = 0.08),data=chr10_genes,colour=rgb(0,0,0,0.25),size=5)+
  geom_segment(aes(x = start, y = 0.08, xend = end, yend = 0.08),data=chr10_exons,colour=rgb(0,0,0,1),size=10)+
  theme(axis.text=element_text(size=14), plot.title=element_text(size=18),axis.title=element_text(size=18,face="bold"),legend.position=c(0.8, 0.9)) 

# candidates?

candidates<-merge(mz,teo,by="WinCenter",suffixes=c(".mz",".teo")) %>% 
  mutate(rat=(tP.mz/nSites.mz)/(tP.teo/nSites.teo),Mb=WinCenter/1E6) %>%
  filter(nSites.mz>1000,nSites.teo>1000,rat<quantile(rat,0.05)) 

mygene=filter(chr10_genes,start<17.8,start>17.7)$gene_id
print(paste('http://www.maizegdb.org/gene_center/gene/', mygene,sep=""))
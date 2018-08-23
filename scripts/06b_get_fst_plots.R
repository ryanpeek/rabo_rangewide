# Plot/model fst vs dist data

library(tidyverse)

# Load Data ---------------------------------------------------------------

# load fst & dist data
load("data_output/fst/final_fst_25k_n3.rda") 

# load metadat
metadat <- read_rds(path = "data_output/rapture_metadata_rabo_quant.rds")

# view summary of data:
metadat %>% group_by(HUC_8) %>% tally %>% print(n=Inf)


# Add EcoRegions and HUCs to Fsts Data ------------------------------------


# add ecoregion and HUC for each site/A/B
metadat_filt <- metadat %>% select(Locality:EcoRegion, lat:HU_8_NAME, county:NHD_Tot_DA_sqkm) %>% 
  mutate(Locality=tolower(Locality)) %>% arrange(Locality) %>% 
  distinct(Locality, .keep_all = T)

# fix trin with space
unique(metadat_filt$Locality)
metadat_filt$Locality<-gsub(pattern = "[[:space:]]", replacement = "-", x = metadat_filt$Locality)


# now join ecoregions:
fst_eco <- left_join(final_fst, metadat_filt[,c(1,4)], by=c("siteA"="Locality")) #site A
fst_eco <- left_join(fst_eco, metadat_filt[,c(1,4)], by=c("siteB"="Locality"))  %>% 
  rename(EcoRegion.A=EcoRegion.x, EcoRegion.B=EcoRegion.y) %>% 
  mutate_at(.vars = c("sitepair", "EcoRegion.A", "EcoRegion.B"), .funs = as.factor)

summary(fst_eco)

# MODEL -------------------------------------------------------------------

mod1 <- fst_eco %>% group_by(EcoRegion.B) %>% sample_n(25)

# model association between fst and dist, loci, and ecoregion

# https://cran.r-project.org/web/packages/brms/vignettes/brms_multivariate.html

library(brms)

fit1 <- brm(
  fst_weight ~ dist_km + no_obs_loci + (1|EcoRegion.B),
  data=fst_eco, chains=2, cores=2
)

# term (1|p|sitepair) indicates a varying intercept over sitepair. By writing |p| in between we indicate that all varying effects should be modeled as correlated

add_ic(fit1) <- "loo"
summary(fit1)
pp_check(fit1, resp = "fst_weight")
bayes_R2(fit1) # how much variation in response variables can be explained by model, Bayesian R^2


marginal_effects(fit1, "dist_km", resp = "fst_weight")



# DO SOME CODING OF VARIABLES ---------------------------------------------

# add healthy watershed metrics here and model against Fst
# ad permeabl, etc
# think about modeling after this step.



# PLOT --------------------------------------------------------------------

#fst_dist_pw$regType <- factor(fst_dist_pw$regType, levels = c("Unregulated","Bypass", "Hydropeaking"))
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# plotly version:

plotly::ggplotly(ggplot() + 
                   stat_smooth(data=fst_eco,
                               aes(y=fst_weight/(1-fst_weight), x=dist_km, color=EcoRegion.A, group=EcoRegion.A), 
                               method = "glm", show.legend = F, lty=2, lwd=1,
                               alpha=0.2, level = 0.89) +
                   geom_point(data=fst_eco, aes(y=fst_weight/(1-fst_weight), x=dist_km,
                                                    fill=EcoRegion.A, group=EcoRegion.A, label=sitepair),
                              show.legend = T, size=4, alpha=0.9) +
                   #scale_fill_manual("Regulation Type", values = c("Unregulated"=cbbPalette[3], "Bypass"=cbbPalette[2], "Hydropeaking"=cbbPalette[7])) +
                   #scale_shape_manual("Regulation Type", values=c(21,22,23))+
                   #scale_color_manual("Regulation Type", values = c("Unregulated"=cbbPalette[3], "Bypass"=cbbPalette[2], "Hydropeaking"=cbbPalette[7])) +
                   theme_bw(base_family = "Roboto Condensed")
)


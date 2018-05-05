
# libraries ---------------------------------------------------------------

library(here)
library(ggplot2)

# GET PCA -----------------------------------------------------------------

source("scripts/f_read_covar.R")

# source the function to plot
source(paste0(here(),"/scripts/f_read_covar.R"))

# load the metadata
load(file = paste0(here(), "/data_output/rapture06_metadata_sierra.rda"))

# set site/reads for bamlist/covar filepaths:
reads <- "25k_thresh"
sitea <- "all_rubyubfea_25k"
covarpath<- paste0(here(), "/data_output/angsd/", sitea, ".covMat")
site <- "all_rubyubfea"
bampath <- paste0(here(), "/data_output/bamlists/", site, "_", reads, ".bamlist")

# run function
(read_covar(covarpath, bampath, metadat, c(2,3),plotlyplot = TRUE))


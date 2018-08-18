
# GET FST Data ------------------------------------------------------------

# first we copy over from either farmer or using scp:
# scp -P 2022 rapeek@agri.cse.ucdavis.edu:/home/rapeek/projects/rangewide/pop_gen/results_fst/all_global_fst.txt .


# Load Libraries ----------------------------------------------------------

library(tidyverse)
library(stringr)


# GET DATA ----------------------------------------------------------------

dat <- read_csv("data_output/fst/all_global_fst.txt", col_names = "alldat")

# skip 1 line of every 3 rows           
dat <- dat[[-seq(from=1, to=10599, by=3),]


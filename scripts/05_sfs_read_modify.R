# load R functions
library(here)
library(viridis)

source(here("scripts/05a_sfs_demography_functions.R"))

# 1) READ AND PLOT OBSERVED GENETIC DATA (2D-SFS) ------------------------

# read in files
(sfsfiles <- tools::file_path_sans_ext(list.files("data_output/sfs/", pattern="\\.2dsfs$")))

# write a loop (cuz I'm too lazy to deal with purrr)
for(i in seq_along(sfsfiles)){
  file2d <- sfsfiles[i]
  file2dname <- gsub(x = file2d, pattern = "[.]", replacement = "_")
  pop1 <- substr(file2d, 1, 8)
  pop2 <- substr(file2d, 10,17)
  # read in file
  sfs2d <- scan(paste0("data_output/sfs/",file2d,".2dsfs", sep=""), quiet=T)
  # transpose and convert to matrix (no indiv in each sfs pop * 2 + 1)
  sfs2d = t(matrix(sfs2d, nrow=31, ncol=31))
  
  # Plot 2DSFS to see everything makes sense
  colv <- viridis(40,direction = -1)
  #colv<-rainbow(n=40,start=0, end =0.3)

  pdf(file = paste0(here(),"/figs/2dsfs_n15_", pop1,"_",pop2,".pdf"), width = 6, height= 4)
  n1 <- nrow(sfs2d)-1
  n2 <- ncol(sfs2d)-1
  image(x=seq(0,n1), y=seq(0,n2), z=-log10(sfs2d), xlab=pop1, ylab=pop2, col=colv)
  dev.off()
}

# 2) GET FOLDED SPECTRUM & FST ---------------------------------------------

# set up fst tibble: 
fsts <- list()

# SET INDIVIDUALS AND CHROMOS
for(i in seq_along(sfsfiles)){
  file2d <- sfsfiles[i]
  file2dname <- gsub(x = file2d, pattern = "[.]", replacement = "_")
  pop1 <- substr(file2d, 1, 8)
  pop2 <- substr(file2d, 10,17)
  # read in file
  sfs2d <- scan(paste0("data_output/sfs/",file2d,".2dsfs", sep=""), quiet=T)
  # transpose and convert to matrix (no indiv in each sfs pop * 2 + 1)
  sfs2d = t(matrix(sfs2d, nrow=31, ncol=31))

  # number of individuals
  r.nr_ind <- (nrow(sfs2d)-1)/2
  c.nr_ind <- (ncol(sfs2d)-1)/2
  
  # number of chromosomes
  r.nr_chrom <- nrow(sfs2d)-1
  c.nr_chrom <- ncol(sfs2d)-1
  
  # set to NA entries in the matrix where the sites is not a SNP
  sfs2d[1,1] <- NA
  sfs2d[(2*r.nr_ind+1), (2*c.nr_ind+1)]=NA

# fold
  sfs2dfold <- fold2DSFS(sfs2d)
  
  # format for fastcoalsim2
  sfs2dfold[is.na(sfs2dfold)]<-0 # replace all NAs with zero
  
  # add col/row headers (sequentially for each col)
  fastcolnames <- sprintf("d0_%02d", 1:ncol(sfs2dfold))
  fastrownames <- sprintf("d1_%02d", 1:nrow(sfs2dfold))
  
  colnames(sfs2dfold) <- fastcolnames
  rownames(sfs2dfold) <- fastrownames
  write.table(sfs2dfold, paste0("data_output/sfs/", file2dname, "_jointMAFpop1_0.obs"), sep='\t', row.names=T, col.names=T)
  
  # plot the folded spectrum
  n1 <- nrow(sfs2dfold)-1
  n2 <- ncol(sfs2dfold)-1
  
  colv <- viridis(40,direction = -1)
  
  png(file = paste0(here(),"/figs/2dsfs_n15_folded_", pop1,"_",pop2,".png"), width = 6, height= 4, units = "in", res = 150)
  n1 <- nrow(sfs2d)-1
  n2 <- ncol(sfs2d)-1
  image(x=seq(0,n1), y=seq(0,n2), z=-log10(sfs2dfold), xlab=pop1, ylab=pop2, col=colv)
  dev.off()
  fst_val <- doFST(sfs2d)
  tmp <- list(fst_estim=fst_val, sitepair=file2dname)
  fsts[[i]] <- tmp
  #return(fsts)
}

# convert to dataframe
fsts2 <- do.call(rbind, lapply(fsts, data.frame))

# write out
readr::write_csv(fsts2, path = "data_output/sfs/2DSFS_derived_fst.csv")

# 3) FIX FOLDED 2DSFS OBS FILE -----------------------------------------------------
# need to add headers / rownames 
# 1 observation
# d0_0    d0_1    d0_2    d0_3    d0_4    d0_5    d0_6    d0_7    d0_8    d0_9    d0_10   d0_11   d0_12



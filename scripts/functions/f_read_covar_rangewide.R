# function to read in covar and plot 

read_covar_range <- function(covmat, # @path to covar file 
                       bamlist, # @path to bamlist
                       metadata, # the name of pre-loaded metadata that contains "Seq" colname with full seq
                       pcs, # the PC to plot as c(1, 2)
                       colvar, # the value for the color
                       shapevar, # the col name for shape: Pop (River), ecoreg, HUC6, HU_8_NAME
                       plotlyplot=TRUE){
  
  if(!require(viridis)) { install.packages("viridis"); require(viridis, warn.conflicts = F)}
  if(!require(ggplot2)) { install.packages("ggplot2"); require(ggplot2, warn.conflicts = F)}
  if(!require(dplyr)) { install.packages("dplyr"); require(dplyr, warn.conflicts = F)}
  if(!require(plotly)) { install.packages("plotly"); require(plotly, warn.conflicts = F)}
  
  # read covar file
  covar <- read.table(covmat, stringsAsFactors = F)
  
  # get bamlist for spp
  bams <- read.table(bamlist, stringsAsFactors = F, header = F) # read in the bamlist
  bams$V2 <- sub('\\..*$', '', basename(bams$V1)) # remove the path and file extension
  annot <- inner_join(bams, metadata, by=c("V2"="Seq")) %>% select(-V1) # join with the metadata
  
  # Eigenvalues
  eig <- eigen(covar, symm=TRUE)
  eig$val <- eig$val/sum(eig$val)
  PC <- as.data.frame(eig$vectors)
  colnames(PC) <- gsub("V", "PC", colnames(PC))
  
  # these are all local variables that can be added...ultimately up to personal pref for plot
  PC$Pop <- factor(annot$River)
  PC$Locality <- factor(annot$Locality) 
  PC$ID <- factor(annot$SampleID)
  PC$HU_12_NAME <- factor(annot$HU_12_NAME)
  PC$HU_8_NAME <- factor(annot$HU_8_NAME)
  PC$HUC6 <- factor(annot$HUC_6)
  PC$LabID <- factor(annot$LabID)
  PC$ecoreg <- factor(annot$EcoRegion)
  
  # set up plot
  
  # PC's:
  pc1 <- pcs[1]
  pc2 <- pcs[2]
  
  # Title: (% explained)
  title <- paste("PC",pc1," (",signif(eig$val[pc1], digits=3)*100,"%)"," / PC", pc2," (",signif(eig$val[pc2], digits=3)*100,"%)", sep="",collapse="")
  
  # plotly
  if(plotlyplot){
    
      
    (plotpca <- ggplotly(p = 
               ggplot(data=PC, 
                      aes_string(x=paste0("PC",pc1), 
                                 y=paste0("PC",pc2), 
                                 color=colvar, 
                                 #shape=shapevar, # can be Pop, HU_8_NAME, HUC6, ecoreg etc 
                                 text=quote((paste0("ID: ", ID, 
                                                    " <br> River: ",
                                                    Locality))))) +
               geom_point(size=4, alpha=0.8) + 
               theme_bw(base_size = 9) +
               scale_color_viridis_d() + 
               ggtitle(paste0(title)))
    )
    return(plotpca)
  
  cat("All Finished! Available in current dataframe...\n")  
  } else {
    
    (ggpca <- ggplot(data=PC, aes_string(x=paste0("PC",pc1), y=paste0("PC",pc2),
                               color=quote(Locality), shape=quote(Pop),
                               text=quote((paste0("ID: ", ID, " <br> River: ", Locality))))) +
    
    geom_point(size=4, alpha=0.8) +
    #ggforce::geom_mark_circle(aes(group=watershed))+
    theme_bw(base_family = "Roboto Condensed") +
    #theme(legend.position="bottom") +
    scale_color_viridis_d("River") + 
    ggtitle(paste0(title)))
    
    return(ggpca)
  }

}


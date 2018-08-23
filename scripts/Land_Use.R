 #First find and read in Land Use data 
#https://www.mrlc.gov/nlcd11_data.php

library(raster)
#setwd("D:/Dropbox/R_data/")

#r <- raster("nlcd_2011_landcover_2011_edition_2014_10_10/nlcd_2011_landcover_2011_edition_2014_10_10/nlcd_2011_landcover_2011_edition_2014_10_10.img")
#https://www.mrlc.gov/nlcd11_leg.php
r <- raster("nlcd_2011_landcover_2011_edition_2014_10_10/nlcd_2011_landcover_2011_edition_2014_10_10.img")

#get to CA or WA 
#download from: http://www.ofm.wa.gov/pop/geographic/tiger.asp
#wa <- shapefile("Dis_Prop/Kelp/kelp_share/state10/state10.shp")
wa <- shapefile("data/state10/state10.shp")
wa <- spTransform(wa, crs(r))
plot(r)
plot(wa, border="red", add=T)
#Clip to WA #use cropNLCD file 

nr <- raster("WA-NLCD.tif")
#table(getValues(nr))

freq(nr)
hist(nr)


raster(rr, "Dis_Prop/Kelp/kelp_share/data/WA_nlcd.tif")

#Use legend to aggregate categories 
#https://www.mrlc.gov/nlcd11_leg.php

#then download watershed shapefile
require(rgdal)

# The input file geodatabase
setwd("/home/ktied/kelp/")
fgdb <- "data/WBD_WA_SP.gdb"

# List all feature classes in a file geodatabase
subset(ogrDrivers(), grepl("GDB", name))
fc_list <- ogrListLayers(fgdb)
print(fc_list)
#we want fourth level or fifth level - is that 8 or 10?
# Read the feature class
fc <- readOGR(dsn=fgdb,layer="WBDHU10")
fc <- spTransform(fc, crs(r))
#For each watershed polygon extract from nr 
#Divide into watershed  
v <- extract(nr, fc, df=TRUE)
v.counts <- lapply(v,table)
saveRDS(v, "extract_lcrast.RDS")
#Aggregate into Ag/Urban/Wild

#s <- merge(fc, v)
#Then divide into watershed  

#then loop through for each year in series 


#Get fancy if this works out -- ag types etc/pop density 

#Fire as well 
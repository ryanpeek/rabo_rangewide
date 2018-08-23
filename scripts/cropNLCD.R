# Cropping really big rasters
library(raster)
library(gdalUtils)

# Make sure you are in the dir with the data
setwd("nlcd_2011_landcover_2011_edition_2014_10_10/")
# Download GADM level 1 for US
us <- getData("GADM", country="USA", level=1)
# Pick a state, subset
aoi <- "CA"
#aoi <- "WA"
#subset and match proj of nlcd EPSG:5070
#aoishp <- spTransform(us[us$HASC_1=="US.WA",],CRS("+init=epsg:5070"))
aoishp <- spTransform(us[us$HASC_1=="US.CA",],CRS("+init=epsg:5070"))
cutlyr <- paste0(aoi,"5070")
cutshp <- paste0(cutlyr,".shp")
shapefile(aoishp, cutshp)

input <- "nlcd_2011_landcover_2011_edition_2014_10_10.img"
output <- paste0(aoi,"-NLCD.tif")
# crop by state
# mask by state
gdalwarp(srcfile=input
         ,dstfile=output
         ,of="GTiff"
         ,multi=TRUE
         ,co=c("COMPRESS=DEFLATE","PREDICTOR=1", "ZLEVEL=6")
         ,ot="Byte"
         ,cutline=cutshp
         ,cl=cutlyr
         ,crop_to_cutline = TRUE
         ,overwrite=TRUE)
# save new raster

nr <- raster(output)
#table(getValues(nr))

freq(nr)
hist(nr)


require(rgdal)

# The input file geodatabase
fgdb <- "C:/path/to/your/filegeodatabase.gdb"

# List all feature classes in a file geodatabase
subset(ogrDrivers(), grepl("GDB", name))
fc_list <- ogrListLayers(fgdb)
print(fc_list)

# Read the feature class
fc <- readOGR(dsn=fgdb,layer="some_featureclass")

#For each watershed polygon extract from nr 

v <- extract(r, polys, df=TRUE)
v
# mean for each polygon - probably wouldn't work 
unlist(lapply(v, function(x) if (!is.null(x)) freq(x, na.rm=TRUE) else NA ))

##########Blah blah blah#####
#Alternate
r <- raster(input)
#maybe don't write this to the file 
q <- crop(r, extent(cutshp), filename= output, 
     options = c("COMPRESS=DEFLATE","PREDICTOR=1", "ZLEVEL=6"), datatype="Byte")


g <- mask(q, cutshp) #or save as a file 
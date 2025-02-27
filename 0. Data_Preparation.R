# 0. Data Preparation
rm(list=ls())
gc()
#.rs.restartR()
options(java.parameters = "-Xmx5g") # increase the memory space for jave before loading any package
options("rgdal_show_exportToProj4_warnings"="none") # Silence packages updates warnings

setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set the working directory to the directory in which the code is stored
td<-tempdir()
dir.create(td,showWarnings = FALSE)
#
#
###//\/\/\/\/><\/\////\/\/\/\/\/\/\///\\\\\\//\/\/\/\////////////////////////></////##-#
##                        Single species SDM pipeline                               ##-#  
###///\/\/\/><\/\/\////\/\/\//\/\/\/\/\/\///\\\\\\\\\..\.\\\\.\\\\\\\\\\\\\\\><\\\\\##-#
#'
#'
# 0. Load the packages----
list.of.packages<-c("sf","terra","tidyverse","geodata","rJava","viridis")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages,require,character.only=TRUE)
rm(list.of.packages,new.packages)

# 0.1 Some project parameters----
# CRS and spatial standarization
crs_p <- "ESRI:54012" # This is an equal area projection system

# 0.1 Load the needed functions----
functions<-"./Functions" %>% list.files(recursive = FALSE,pattern = ".R$",full.names = TRUE)
lapply(functions,function(x) source(x))

# 1. Download the environmental data ----
# a. Defining the Study area----
# Dowload and format the spatial data
"./Data/Spatial" %>% dir.create(recursive = TRUE,showWarnings = FALSE) # Create the output folder

study_area <- geodata::world(resolution=5,path="./Data/Spatial") # get the country borders https://gadm.org/
class(study_area) # Spatvector object (terra compatible)
study_area <- study_area %>% st_as_sf() # Transform into an sf object for processing and manipulation

class(study_area) ; par(mfrow=c(2,1));plot(study_area %>% st_geometry(),main="WGS84",col=viridis::plasma(nrow(study_area)))
study_area <- study_area %>% st_transform(crs=crs_p) ; plot(study_area %>% st_geometry(),main="Eckert IV equal area",col=viridis::viridis(nrow(study_area)))

# Delimit our study area (in this case we are interested in Africa)
iso_doces_by_region <- read.csv("./Data/Complementary_info/ISO_codes.csv") # World Iso Codes with other spatial informaiton, from -> https://github.com/lukes/ISO-3166-Countries-with-Regional-Codes/blob/master/all/all.csv
head(iso_doces_by_region)

ISO_Africa <- iso_doces_by_region %>% filter(region=="Africa") # Select the iso codes of Africa
head(study_area$GID_0) ; head(ISO_Africa$alpha.3)

study_area <- study_area %>% filter(GID_0 %in% ISO_Africa$alpha.3)
Africa_pol <- study_area %>% st_union() # Combine all the polygons into a single spatial object

par(mfrow=c(1,2))
plot(study_area %>% st_geometry(),col=viridis::mako(n=nrow(study_area)),main="Country Selection") ; 
plot(Africa_pol %>% st_geometry(),col="firebrick",main="Merged Polygons")

# Export the spatial information
"./Data/Spatial/study_area" %>% dir.create(recursive = T,showWarnings = F) # Create the exit folder
st_write(study_area,paste("./Data/Spatial/study_area","Africa_countries.shp",sep="/"),append = FALSE)
st_write(Africa_pol,paste("./Data/Spatial/study_area","Africa_pol.shp",sep="/"),append = FALSE)

# b. Dowload the climatic environmental information----
# We want this to be relatively fast, so lets take the 2.5 min (degrees) resolution
bio_clim <- "./Data/Spatial/BioClim" ; dir.create(bio_clim,recursive = TRUE,showWarnings = FALSE) # Create the exit folder
geodata::worldclim_global(res=2.5,path=bio_clim,var="bio") # Download the data for the whole world

# Dowload spatial information for an African region
iso_region <- ISO_Africa %>% filter(sub.region=="Northern Africa") %>% dplyr::select("alpha.3") %>% unlist() %>% unname()
lapply(iso_region,function(w) geodata::worldclim_country(country=w,res=2.5,path=bio_clim,var="bio")) # Download the data for the whole world

# c. Format the data to meet the study area requirements ---- 
# The spatial information is in different CRS and projection systems, we need to reproject the environmental information to follow our
# standards. We can do this in two ways:
# 1. Provide a CRS or CRS description to the terra::project function
# 2. Create a reference raster that follows the CRS and extension we want to use 
#
#  Method 1: Since the initial raster is too big, we are going to crop it to match the size of our study area
bioClim_G <- list.files("./Data/Spatial/BioClim/wc2.1_2.5m", pattern=".tif$",full.names = T)
bioClim_G <- lapply(bioClim_G,rast) ; bioClim_G <- do.call("c",bioClim_G)

bioClim_A <- bioClim_G %>% crop(study_area %>% st_transform(crs="EPSG:4326") %>% vect(),mask=T) # EPSG:4326 is the identifier of WGS 84 -> https://epsg.io/4326
bioClim_A_Equal <- bioClim_A %>% project(y=crs_p) # Project the data into the desire CRS

# Method 2: Create a reference raster and use it to reproject the data
ref_rast <- rast(x=study_area,ncols=ncol(bioClim_A_Equal),nrows=nrow(bioClim_A_Equal),crs=crs_p,vals=NA)#,extent=terra::ext(bioClim_G))
bioClim_A_Equal2 <- bioClim_G %>% project(y=ref_rast) %>% mask(study_area %>% vect())

# This method requires less steps and in some ways is more efficient, however it requires more information such as the dimensions of the final
# raster object and the spatial context (crs, extension, and location) for the translocation of values. In this example this context is given 
# by the study area polygon but despite this and using the same number of columns and rows as the already rescalled and reprojected raster (bioClim_A_Equal)
# there is still a discrepancy between the resolution of the transformed rasters. So this method is only consistent if we generate a detailled reference raster
# e.g
bioClim_A_Equal3 <- bioClim_G %>% project(y=bioClim_A_Equal) %>% mask(study_area %>% vect())

lt<-layout(matrix(c(1,1,1,1,2,3,4,5),nrow=2,ncol=4,byrow=T))
layout.show(lt)
plot(bioClim_G$wc2.1_2.5m_bio_1,main="Original data") ; 
plot(bioClim_A$wc2.1_2.5m_bio_1,main="Cropped OR data")
plot(bioClim_A_Equal$wc2.1_2.5m_bio_1,main="Cropped and projected",col=viridis::magma(200),ces=0.7)
plot(bioClim_A_Equal2$wc2.1_2.5m_bio_1,main="Projected ref-rast",col=viridis::inferno(200),ces=0.7)
plot(bioClim_A_Equal3$wc2.1_2.5m_bio_1,main="Projected ref-rast\nPrecise",col=viridis::rocket(200),ces=0.7)

par(mfrow=c(1,1))

# d. What if our study area is made of multiple raster files----
# In many occassions remote sensing data is organized in tiles or segments. Sometimes we can work with the individual tiles, but in most occasions we need to
# fuse/merge the spatial data to meet our study region
list.dirs(bio_clim) # Get the BioClim country data
# bioClimReg <- list.files("./Data/Spatial/BioClim/wc2.1_country",pattern=".tif$",recursive = T,full.names = T) %>% rast()

# We cannot combine all the layers into a single raster stack, like the other data, since each object has different extends.
# We are going to load each raster individually and stored into a list:
pol_region <- study_area[study_area$GID_0 %in% iso_region,] 

tiles <- list.files("./Data/Spatial/BioClim/wc2.1_country",pattern=".tif$",recursive = T,full.names = T)
tiles_r <- lapply(tiles,rast)

# Transform into a SpatRasterCollection
tiles_r <- sprc(tiles_r) # the sprc transform our list into a spat raster colection
bioClimReg_comb <- mosaic(tiles_r)# Combines the tiles and if they partially overlap calculates the mean for those overlapping pixels

bioClimReg_comb <- bioClimReg_comb %>% project(y=crs_p) # Reproject the data to meet the same CRS
bioClimReg_comb <- bioClimReg_comb %>% mask(pol_region %>% vect()) # Use the regional geometries to mask the data

# Display the data
bioClimReg_comb[[1]] %>% plot(col=viridis::viridis(200), main="Combined raster",axes=F,box=F)
pol_region %>% st_geometry() %>% plot(add=T)

# e. Get some elevation data----
geodata::elevation_global(res=2.5,path="./Data/Spatial/Elevation")
elev <- "./Data/Spatial/Elevation" %>% list.files(full.names = TRUE,pattern=".tif$",recursive = T) %>% rast()

greyFun<-colorRampPalette(c("grey90","grey30"))
elev %>% plot(col=greyFun(100))

# Same as with the Bioclim data we are only interested in our study area. We are going to use the already adapted Bioclim data
# to re-sample (adapt dimension, resolution,  and distribution of cells)
elev_equal <- elev %>% project(y=ref_rast) %>% mask(Africa_pol %>% vect())
elev_equal %>% plot(col=greyFun(100))

elev_equal <- elev_equal %>% resample(y=ref_rast) # We are making sure that the cells of the raster are 

# With the elevation data we can calculate other metrics
terrain_lyrs <- terra::terrain(elev_equal,v=c("slope","aspect","roughness")) # Calculate the slope, aspect and terrain roughness
hill_shade <- terra::shade(slope=terra::terrain(elev_equal,v=c("slope"),unit="radians"),
                                                aspect=terra::terrain(elev_equal,v=c("aspect"),unit="radians")
                           ) # Calculate the hillshade map for plotting porpuses. Need to specify the units of slope and aspect as radians

# f. Get some future climatic data for forecasting----
clim_change_route <- "./Data/Spatial/BioClim/Scenarios" ; clim_change_route %>% dir.create(showWarnings = FALSE,recursive = TRUE)
scenarios<-c("126","245","370","585")

for(i in scenarios){
  save_route<-paste(clim_change_route,i,sep="/") ; save_route  %>% dir.create(recursive = TRUE,showWarnings = FALSE)
  geodata::cmip6_world(model="MPI-ESM1-2-HR",ssp=i,time="2061-2080",var="bioc",res=2.5,path=save_route)
}

# Process the spatial information
scenarios_route <- clim_change_route %>% list.files(pattern=".tif$",recursive = TRUE,full.names = TRUE)

for(i in scenarios){
  # Load the scenarios
  x<-scenarios_route[scenarios_route %>% grepl(pattern=i)] %>% rast()
  names(x) <- names(x) %>% gsub(pattern=paste0("_MPI-ESM1-2-HR_ssp",i,"_2061-2080"),replacement="") %>% gsub(pattern="bioc_",replacement="bio_")
  
  # Adapt the layer
  x <- x %>% project(y=ref_rast) %>% mask(Africa_pol %>% vect())
  x <- x %>% resample(y=ref_rast)
  
  # Export the information
  exit_lyr <- paste("./Data/Spatial/Processed/Scenarios",paste0("spp",i),sep="/")
  dir.create(exit_lyr,showWarnings = FALSE,recursive = TRUE)
  x %>% writeRaster(paste(exit_lyr,paste0("Bio_Spp_",i,".tif"),sep="/"),overwrite=T)
  
  plot(x$wc2.1_2.5m_bio_1)
}


# 2. Load and process the IUCN range information----
# The IUCN Red-List of Threatened species contains spatial data regarding the range and distribution for hundred of thousands 
# of species. We are going to load a small sample of this spatial information for the species we are going to model.
# a. Get the names of the species selected for the analysis (These are their names as they appear int the IUCN)

Species_analysis <- c(#"Cephalophus dorsalis","Cephalophus natalensis","Cephalophus rufilatus",  
                      #"Cephalophus silvicultor",
                      "Eidolon helvum")
                      #,"Epomophorus gambianus",                     
                      #"Gorilla beringei","Epomophorus pusillus","Pan troglodytes",
                      #"Rousettus aegyptiacus")

# Get the Spatial data from a local folder containing the IUCN spatial information (skip)
route_IUCN_pols <- "./Data/SP_Info/IUCN_range" %>% list.files(pattern=".shp$",recursive=T,full.names = TRUE)
head(route_IUCN_pols)

IUCN_pols <- lapply(Species_analysis, function(x) grep(pattern = x,route_IUCN_pols,value=T))
IUCN_pols<-IUCN_pols[-8] # We don't have spatial information for one species!

# Load and process the species ranges
IUCN_pols <- lapply(IUCN_pols,st_read)
IUCN_pols <- do.call("rbind",IUCN_pols) # We can combine all the spatial features into a single object
IUCN_pols <- IUCN_pols %>% st_transform(crs=crs_p) # IUCN polygons use a different crs than the project, we need to re-project the geometries

IUCN_pols # The sf package format the spatial data in a way that resembles a data.frame with a special object or colum named geometry, that contains the spatial information
IUCN_pols %>% st_drop_geometry() # If we drop the geometry we end up with a data.frame

# diplay the range data
plot(Africa_pol %>% st_geometry())
IUCN_pols %>% st_geometry() %>% plot(col=viridis::viridis(n=nrow(IUCN_pols))%>%adjustcolor(alpha.f = 0.15),add=T)

# 3. Export all the processed information----
# a. Study area and range information----
"./Data/Spatial/Processed/Vector" %>% dir.create(recursive = T,showWarnings = F)

Africa_pol %>% st_write(paste("./Data/Spatial/Processed/Vector","StdArea.shp",sep="/"),append=F)
study_area %>% st_write(paste("./Data/Spatial/Processed/Vector","CountryPols.shp",sep="/"),append=F)
pol_region %>% st_write(paste("./Data/Spatial/Processed/Vector","pol_region.shp",sep="/"),append=F)

# b. Species ranges----
# IUCN_pols %>% st_write(paste("./Data/Spatial/Processed/Vector","IUCN_ranges.shp",sep="/"),append=F)

# c. Environmental information----
"./Data/Spatial/Processed/Raster" %>% dir.create(recursive = T,showWarnings = F)

bioClim_A_Equal2 %>% writeRaster(paste("./Data/Spatial/Processed/Raster","Bioclim.tif",sep="/"),overwrite=TRUE)
elev_equal %>% writeRaster(paste("./Data/Spatial/Processed/Raster","Elevation.tif",sep="/"),overwrite=TRUE)
terrain_lyrs %>% writeRaster(paste("./Data/Spatial/Processed/Raster","Terrain.tif",sep="/"),overwrite=TRUE)
hill_shade %>% writeRaster(paste("./Data/Spatial/Processed/Raster","Hill_shade.tif",sep="/"),overwrite=TRUE)


# 3. Display the spatial information----
#~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~
# Some parameters
colors <- colorRampPalette(c("navy","skyblue","tomato")) ; colors<-colors(1000)
results_r <- "Results/Figures/" ; results_r %>% dir.create(recursive=TRUE,showWarnings = FALSE)

png(paste(results_r,paste0("Mock_map",".png"),sep="/"),width=15,height = 15,res=600,units="cm")
par(bg="#fffbe3ff")

# We are going to use the hill shade as the base layer for our map
hill_shade %>% 
  plot(axes=F,col=greyFun(1000) %>% adjustcolor(alpha.f = 0.5),legend=F,
       plg = list(loc = "bottom", size=c(0.5,1), title = "Temperature"))

# Overlay the temperature data and add the legend at the bottom of the map
bioClim_A_Equal2$wc2.1_2.5m_bio_1 %>% 
  plot(axes=F,col=colors %>% adjustcolor(alpha.f = 0.5),
       plg = list(loc = "bottom", size=c(0.5,1), title = "Temperature"),add=T)

# Add some extra spatial information
plot(pol_region %>% st_geometry(),add=TRUE,border="gold")
Africa_pol %>% st_geometry() %>% plot(border="black",add=TRUE)

plot(IUCN_pols %>% filter(SHAPE_AREA==max(SHAPE_AREA)) %>% 
       st_geometry(),add=TRUE,pch=19,cex=0.15,border="firebrick4")

# Legend
legend("left",legend=c(paste("Distribution of\n", 
                       IUCN_pols %>% st_drop_geometry() %>% 
                         filter(SHAPE_AREA==max(SHAPE_AREA)) %>% 
                           dplyr::select("BINOMIAL")),
                              "Study area"),col=c("firebrick4","gold"),pch=0,cex=1.2,bty="n")

dev.off()

# End of the script

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
##                        Process species spatial data                              ##-#  
###///\/\/\/><\/\/\////\/\/\//\/\/\/\/\/\///\\\\\\\\\..\.\\\\.\\\\\\\\\\\\\\\><\\\\\##-#
#'
#'
# Load the packages needed for the analysis
# 0. Load the packages needed ----
list.of.packages<-c("dplyr","sf","data.table","terra","tidyr","parallel","purrr",
                    "data.table","doParallel","rJava","rstudioapi","raptr")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages,require,character.only=TRUE)
rm(list.of.packages,new.packages)

# 0.1 Load the functions to run the analysis ----
functions<-"./Functions" %>% list.files(recursive = FALSE,pattern = ".R$",full.names = TRUE)
lapply(functions,function(x) source(x))

# 0.2 Other parameters ----
# CRS and spatial standarization
crs_p <- "ESRI:54012" # This is an equal area projection system

# 1. Load the species spatial information ---- 
# a. List of species
Species_list <- c("Cephalophus dorsalis","Cephalophus natalensis","Cephalophus rufilatus",  
                                      "Cephalophus silvicultor","Eidolon helvum","Epomophorus gambianus",                     
                                      "Gorilla beringei","Epomophorus pusillus","Pan troglodytes",
                                      "Rousettus aegyptiacus")

# b. Species ranges
Sp_ranges <- "./Data/Spatial/Processed/Vector/IUCN_ranges.shp" %>% st_read()

# c. Study area
stdy_area <- "./Data/Spatial/Processed/Vector/StdArea.shp" %>% st_read()

# 2. Get the spatial information for the species:
# We are going to dowload the species distribution data, this is the presence points, directly from the servers of GBIF 
# using the functions built in the rgbif package. However, since we are going to collect spatial information spanning several decades,
# we need to check the taxonomic information of our species in order to see their potential synonyms.
#
# a. Taxonomic checking ----
# We are going to use the retrive_syns function. This is a wrapper function of taxize package that automatically search multiple taxonomic repositories and 
# returns a unified taxonomi for each species along with a list of references and accepted synonyms. 

# We need an API key to access the IUCN red list information
options(iucn_redlist_key="hA36bu9uSrSsxPexjAmtyw8ocf9opEyTM9tF") # We need an api token from the IUCN to access the IUCN information
                                                                 # https://apiv3.iucnredlist.org/
spp_tax <- lapply(Species_list,retrieve_syns)
tax_info <- lapply(spp_tax,function(l) return(l[["TaxDat"]])) %>% rbindlist()

sp_list <- list()

for(w in 1:nrow(tax_info)){
    sp_y <- tax_info[w,]
    sp_y <- c(sp_y$Or_name,sp_y$ITIS_name,sp_y$IUCN_syn,sp_y$ITIS_syn,sp_y$ITIS_name) %>% na.omit()
    
    sp_list[[w]] <- lapply(sp_y,strsplit,split=";") %>% unlist() %>% unique()
}

write.csv(tax_info,paste("./Data/Sp_Info","Taxonomic_information.csv",sep="/"))
tax_info <- read.csv(paste("./Data/Sp_Info","Taxonomic_information.csv",sep="/"))

#
# b. Retrieve the Gbif data ----
# We are going to use the Dowload_gbif function to dowload the species spatial records from Gbif. This is a wrapper around
# the rgbif package. The function launch multiple calls to the Gbif servers to download records for each year and month of the total
# for the specified time-span. This way we can skip the hard limit that rgbif impose on the size of the records that can be downloaded at
# each time.
"./Data/SP_info/Sp_records" %>% dir.create(showWarnings = F,recursive = T)

for(i in tax_info$Or_name){
  # Download spatial records
    gbif_r <- try(rgbif::occ_search(scientificName = i,year = "2000,2025",limit=100000),silent=F)
    gbif_r[["data"]] %>% as.data.frame() %>% write.csv(paste("./Data/SP_info/Sp_records",paste0(i,".csv"),sep="/"),row.names=F)
    
    print(i)
  }

# c. Clean the spatial information ----
sp_routes <- list.files("./Data/SP_Info/Sp_records",full.names = TRUE,pattern = ".csv$")
"./Data/SP_Info/Sp_Records/Clean" %>% dir.create(showWarnings = FALSE,recursive = TRUE) # Exit route

for(i in 1:length(sp_routes)){
  # sp name
  sp_01 <- sp_routes[i] %>% basename() %>% gsub(pattern=".csv",replacement="")
  
  # Load the data and transform into a sf object
  y_points <- sp_routes[i] %>% read.csv() #
  y_points <- y_points[!y_points$decimalLongitude %>% is.na(),]
  
  # Filter the points using coordinate cleaner:
    # We need to get the coordinates and other spatial information into
    # longitude and latitude in order for the function to correctly classify and 
    # remove the erroneous points
  
  if(TRUE %in% c(Sp_ranges$BINOMIAL %in% sp_01)){
        rX <- Sp_ranges[Sp_ranges$BINOMIAL %in% sp_01,] # get the range data
        rX <- rX %>% st_union() %>% st_transform(crs="EPSG:4326")
        
        # Check the polygon
        if(!rX %>% st_is_valid()) rX <- NULL
  
      }else{
        rX <- NULL
      }
  
  clean.p <- Prepare_points(y_points, range_sp=rX, xy.c=c("decimalLongitude","decimalLatitude"), b.width=1.5, crs.r="EPSG:4326")
  
  # With the clean points, transform and adapt to meet the rest of the spatial data
  clean.p <- clean.p %>% st_as_sf(coords=c("decimalLongitude","decimalLatitude"),crs="EPSG:4326") # Transform into an sf object
  clean.p <- clean.p %>% st_transform(crs=crs_p) # Re-project the CRS of the spatial object
  clean.p <- cbind(clean.p %>% st_drop_geometry(),st_coordinates(clean.p))# Load the coordinates
  
  # Export the data
  clean.p %>% write.csv(paste("./Data/SP_Info/Sp_Records/Clean",paste0(sp_01,".csv"),sep="/"),row.names = FALSE)
  }

# 3. Let display the spatial information ----
# Load all the point information
all_points<-lapply(list.files("./Data/SP_Info/Sp_Records/Clean",pattern = ".csv$",full.names = TRUE),read.csv) %>% rbindlist(fill=TRUE) %>% as.data.frame()
all_points<-all_points[,-c(1)] %>% st_as_sf(coords=c("X","Y"),crs=crs_p)

plot(stdy_area %>% st_geometry())
plot(all_points["species"],pch=19,add=T)

# End of the script

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
##                              SDM - Future forecasting                            ##-#  
###///\/\/\/><\/\/\////\/\/\//\/\/\/\/\/\///\\\\\\\\\..\.\\\\.\\\\\\\\\\\\\\\><\\\\\##-#
#'
#'
# Load the packages needed for the analysis
# 0. Load the packages needed ----
list.of.packages<-c("dplyr","sf","data.table","terra","tidyr","parallel","purrr","MASS","nnet","biomod2",
                    "data.table","doParallel","rJava","rstudioapi","raptr","dismo","randomForest","caret")

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
# a. List of species----
Species_list <- list.files("./Data/SP_Info/Sp_records/Clean",full.names = TRUE)

# b. Species spatial information----
IUCN_ranges <- "./Data/Spatial/Processed/Vector/IUCN_ranges.shp" %>% st_read()
Stdy_area <- "./Data/Spatial/Processed/Vector/StdArea.shp" %>% st_read()

# c. Load the present information----
env_lyrs <- lapply("./Data/Spatial/Processed/Raster" %>% list.files(patter=".tif$",full.names = TRUE),rast)
env_data <- do.call("c",env_lyrs[-c(3)])

plot(env_data[[c(1:3,21:23)]],axes=F)

# d. Load the future scenarios data
terrain.lyrs <- do.call("c",env_lyrs[-c(1,3)])

fut_env <- lapply("./Data/Spatial/Processed/Scenarios" %>% list.files(patter=".tif$",full.names = TRUE,recursive = TRUE),rast)
names(fut_env)<-basename("./Data/Spatial/Processed/Scenarios" %>% list.files(patter=".tif$",full.names = TRUE,recursive = TRUE))

fut_env <- lapply(fut_env,function(x) c(x,terrain.lyrs))

panel(c(fut_env$Bio_Spp_126.tif$wc2.1_2.5m_bioc_1,
          fut_env$Bio_Spp_245.tif$wc2.1_2.5m_bioc_1,
            fut_env$Bio_Spp_370.tif$wc2.1_2.5m_bioc_1,
              fut_env$Bio_Spp_585.tif$wc2.1_2.5m_bioc_1),col=viridis::inferno(200))

# 2. load the saved models----
models<-readRDS(list.files(paste("./Results/SDM_models","Eidolon helvum",sep="/"),pattern=".rds$",full.names = TRUE))
     
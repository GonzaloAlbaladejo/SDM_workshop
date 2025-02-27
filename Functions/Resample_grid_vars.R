
########################################################
# This function takes a list of .tif or other compatible raster object locations and 
# adapt the objects to a particular extend, resolution and CRS. Layers are stacked and saved as
# and individual spatial object into a results folder
########################################################
#
resample.rast<-function(x,# Reference raster (should be a spatrast object or matrix)
                        res=NULL, # If x is missing we can built the reference raster using the resolution (res), extend(ex) and crs (crs.r)
                        ex=c(-180, 180, -60, 90), # if x is a matrix, the extent of the raster must be specified, by default it takes the value of the whole globe
                        crs.r="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", # Coordinate reference system, default WGS84 (lat~lon)
                        y, # filepath to the layers we want to re-size/sample
                        paralell=FALSE, # Use multiple when there is more than one object to resize? default=FALSE
                        cores=2,      # If paralell==TRUE, how many cores do you want to commit? default==2
                        route.r= NULL, # route to save the temporal re-sampled layers, by default it would create a temporal folder in the working environment with the name of r_resample
                        results.r=NULL, # Route to save the results, by default it create a ResultsRast folder in the working environment
                        mask=FALSE, # Should the final rasters need to be masked? default FALSE if TRUE a wrld_simpl dataset is usded to mask the rasters
                        mask_r=NULL # polygon to use as mask, only if mask=TRUE
                        ){
   # Load the required libraries to run the function
   # a. Load the packages needed----
   list.of.packages<-c("terra","doParallel","dplyr","foreach","sf")
   
   new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
   if(length(new.packages)) install.packages(new.packages)
   
   lapply(list.of.packages,require,character.only=TRUE)
   rm(list.of.packages,new.packages)
   
   # b. Creates the Temporal file folder and the folder to allocated the results
   if(is.null(route.r)){Tdir<-paste(getwd(),"r_resample",sep="/")}else{Tdir<-route.r}
   if(!dir.exists(Tdir)){dir.create(Tdir)}
   
   if(is.null(results.r)){Resdir<-paste(getwd(),"ResultsRast",sep="/")}else{Resdir<-results.r}
   if(!dir.exists(Resdir)){dir.create(Resdir)}
   
   # Tmp dir---disposable
   Tmp_dir<-paste(getwd(),"junk",sep="/")
   if(!dir.exists(Tmp_dir)){dir.create(Tmp_dir)}
   
   # c. Resample the selected layers
   # c.1 Reference raster, loaded or create it
   if(missing(x)){
      ref.rast<-rast(resolution=res,crs=crs.r,extent=ex,vals=NA)
   }else{
      if(is.matrix(x)){
         ref.rast<-rast(x,crs=crs.r,extent=ex)
      }
      if(class(x)=="SpatRaster"){
         ref.rast<-x
      }}
   
   # d. Resample the layers
   # d.1 If paralell processing is set to true and the number of layers is greater than 1
   if(length(y)>1 & paralell==TRUE){
      
      if(cores>detectCores(logical = FALSE)){
         cores<-detectCores(logical = FALSE)-1   
      }
      
      # Need to write the reference rast into disk to distribute it to the workers
      terra::values(ref.rast)<-NA
      terra::writeRaster(ref.rast,
                         filename=paste(Tmp_dir,paste0("ref_rast.tif"),sep="/"),
                         overwrite=TRUE)
      
      # Create the cluster structure
      cl<-makeCluster(cores)
      registerDoParallel(cl)
      
      foreach(k=1:length(y)) %dopar% {
         list.of.packages<-c("terra","dplyr")
            new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
               if(length(new.packages)) install.packages(new.packages)
         
         lapply(list.of.packages,require,character.only=TRUE)
         rm(list.of.packages,new.packages)
         
         ref.rast2<-paste(Tmp_dir,paste0("ref_rast.tif"),sep="/") %>% rast() # Load the reference raster
         r.rast<-rast(y[k]) # load the data to resample
         
         # d.2 Configure Terra to load the processing into the disk rather than the ram memory
         terraOptions(memfrac=0.9,todisk=TRUE,tempdir = Tmp_dir) # The object y will be allocated in the temporal folder of 
         terra::resample(x=r.rast,y=ref.rast2,filename=paste(Tdir,paste0(names(r.rast),".tif"),sep="/"),overwrite=TRUE)
         gc()
         
         print("SUCCESS!!")
      }
   stopCluster(cl)
   
   }else{
      
      if(length(y)>1){
         
         for(i in 1:length(y)){
            r.rast<-y[i] %>% rast()
            terraOptions(memfrac=0.9,todisk=TRUE,tempdir = Tmp_dir) # The object y will be allocated in the temporal folder of 
            terra::resample(x=r.rast,y=ref.rast,filename=paste(Tdir,paste0(names(r.rast),".tif"),sep="/"),overwrite=TRUE)
            print("SUCCESS!!")
         }
         gc()
      
         }else{
      
      terraOptions(memfrac=0.9,todisk=TRUE,tempdir = Tdir)
      r.rast<-rast(y)
      terra::resample(x=r.rast,y=ref.rast,filename=paste(Tdir,paste0(names(r.rast),".tif"),sep="/"),overwrite=TRUE)
      print("SUCCESS!!")
      }
   }
   
   # e. Unify all the re-sample layers into a single stack
   f.r<-list.files(Tdir,pattern = ".tif$",full.names = TRUE,recursive = FALSE) # put always the full route!!!
   resampled.r<-rast(f.r)
   
   if(mask==TRUE){
      x.mask<-mask_r %>% vect() # need to transform the sf object into a spatial vector object
      resampled.r<-terra::mask(resampled.r,mask=x.mask,touches=TRUE)   
      }
   
   terra::writeRaster(resampled.r,
                      filename=paste(Resdir,paste0("Resample_rast.tif"),sep="/"),
                      names=basename(f.r %>% gsub(pattern = ".tif$",replacement ="")),
                      overwrite=TRUE) # Export the full results as a stack of layers
   
   unlink(Tdir,recursive = TRUE)
   unlink(Tmp_dir,recursive = TRUE)
   print("ALL DONE!")
   return(basename(f.r %>% gsub(pattern = ".tif$",replacement =""))) # these are the names of the layers
}

# End of the function

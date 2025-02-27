#
#################################################################################
#       dowload species GBIF point distribution for a given set of species      #-//\/\/\/\/\/\
#################################################################################
#
# This script automatically dowload the Gbif data for a list of species on an specified area
#
#
# 
Dowload_gbif<-function(sp_list, # (Character) List of species from which to dowload spatial information
                       initial_date, # (Numeric/year) By default the function will dowload 500 records for each month of the year, from the year specified till present
                       # n_cores=2, # (Numeric) Number of cores commited to the processing (only when sp_list>1)
                       exit_route, # (Character) Route to store the dowloaded information
                       area=NULL, # (character) Searches for occurrences inside a polygon in Well Known Text (WKT) format. A WKT shape written as either
                       gadm_codes=NULL, # (character) The gadm id of the area occurrences are desired from. https://gadm.org/.
                       locality=NULL, # If iterating around different polygons, the name of the regions or polygon
                       n_records=150000 # Maximun number of records to retrieve
                       ){
  
  # 0.Load the required libraries----
  list.of.packages<-c("sp","rgbif","raster","data.table","dplyr","doParallel","parallel")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  rm(list.of.packages)

  # Set the dates for the searching
  present_date<-Sys.time() %>% format("%Y") %>% as.numeric()
  years=c(initial_date:as.numeric(present_date))
  
  # Create the output folder for the data
  if(!dir.exists(exit_route)){dir.create(exit_route,showWarnings = FALSE)}

  if(length(sp_list)<2){
    
    for (k in 1:length(years)){
      print(paste0("year=",k))
    
    for (j in 1:12){
      
      points <- NULL
      t_11 <- 1
      while(is.null(points) && t_11 <= 1000) {

        points<-try(occ_search(scientificName =sp_list,
                               hasCoordinate=TRUE,
                               year=years[k],
                               month=j,
                               hasGeospatialIssue=FALSE,
                               limit=n_records,
                               geometry=area,
                               gadmGid=gadm_codes),
            silent = FALSE)
        t_11 <- t_11 + 1
      }
      rm(t_11)
      # 
      # points<-occ_search(scientificName =sp_list,
      #                   hasCoordinate=TRUE,
      #                   year=years[k],
      #                   month=j,
      #                   hasGeospatialIssue=FALSE,
      #                   limit=n_records,
      #                   geometry=area,
      #                   gadmGid=gadm_codes)

      if(is.null(points$data)){
        next
      }else{
        if(!exists("y")){
          y<-points$data
        }else{
          y<-data.table::rbindlist(list(y,points$data),fill=TRUE)
        }}
      rm(points)
      print(paste(k,j,sep="-"))
    }
  }
  gc()
  
  if(exists("y")){
    
    if(is.null(sp_list)){
      sp_list<-paste("All_records",locality)
    }
      
     write.csv(y,paste(exit_route,paste0(sp_list,".csv"),sep="/"),row.names = FALSE)
      rm(y)
    }
  
}else{
  # Create a temporal folder to host the information until its combined
  Temp_tax<-paste(getwd(),"Temp_tax",sep="/")
  dir.create(Temp_tax)
  # Run a loop in parallel to extract the information
  # if the number of cores is greather than the number of species,
  # set n_core to the same length of the species
  # if(n_cores > length(sp_list)){
  #   n_cores<-length(sp_list)
  # }
  # 
  # cl<-makeCluster(n_cores)

  #foreach(f=1:length(sp_list)) %dopar% {
  # # load the library to the workers
  # 
  #   list.of.packages<-c("sp","rgbif","raster","data.table","dplyr","doParallel","parallel")
  #   lapply(list.of.packages,require,character.only=TRUE)
  #   rm(list.of.packages)

  for(f in 1:length(sp_list)){
    t1<-Sys.time()
    sp<-sp_list[f]

    for (k in 1:length(years)){
      print(paste0("year=",k))

      for (j in 1:12){

        points <- NULL
        t_11 <- 1

        while(is.null(points) && t_11 <= 1000) {

          points<-try(occ_search(scientificName =sp,
                                         hasCoordinate=TRUE,
                                         year=years[k],
                                         month=j,
                                         hasGeospatialIssue=FALSE,
                                         limit=n_records,
                                         geometry=area,
                                         gadmGid=gadm_codes),
              silent = FALSE)
          t_11 <- t_11 + 1
        }
        rm(t_11)

        # points<-occ_search(scientificName =sp,
        #                    hasCoordinate=TRUE,
        #                    year=years[k],
        #                    month=j,
        #                    hasGeospatialIssue=FALSE,
        #                    limit=n_records,
        #                    geometry=area,
        #                    gadmGid=gadm_codes)

        if(is.null(points$data)){
          next
        }else{
          if(!exists("y")){
            y<-points$data
          }else{
            y<-data.table::rbindlist(list(y,points$data),fill=TRUE)
          }}
        rm(points)
        print(paste(k,j,sep="-"))
      }
    }
    gc()

    if(exists("y")){
      write.csv(y,paste(Temp_tax,paste0(sp,".csv"),sep="/"),row.names = FALSE)
      rm(y)
    }
   }
  #stopCluster(cl)
  # Combine the files from the temporal folder, export it to the relevant folder and 
  # erase the temporal folder
    points_spp<-lapply(list.files(Temp_tax,pattern = "csv",full.names = TRUE),fread)
    points_spp<-points_spp %>% rbindlist(fill=TRUE)
    
    write.csv(points_spp,paste(exit_route,paste0(sp_list[1],".csv"),sep="/"),row.names = FALSE) # Export the point information under the original name
    unlink(Temp_tax,recursive=TRUE) # Erase the temporal data
  
  }
 return("GBIF data downloaded") 
}

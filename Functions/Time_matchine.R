##\/\/\/\/\/?|/\/?\/|??|??\/\//\/\/\/\/\/\/\/\/\/\/\/\/\\\\\\\\\\\\\\/\/\/\/\////////////////\/\/\/\/\//
#                                                                                                     //
#                                         Time_MatchIne                                               //
#         To coordinate spatio temporal environmental data with species time referenced point data    //           
#                                                                                                     //    
##\/\/\\/\/\/\/\/\/?|/\//\/\/\//\/\/\/\/\/////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\/\/\/\/\/\/\///

Time_matchine<- function(x, #[[RASTER]] Raster stack containing the formatted time and spatial environmental data
                         y, # [[SF]] An sf object containing the spatial [POINT] presence data for the species or group of species of interest
                         y.t = "eventDate", # [[Character]] Field in the sf object 
                         time.i = "year", # [Character] It can take the form of years, months or days
                         jump.i = 1, # [numeric] The jump in the intervals
                         id.p = "key",# [Numerical] the id that identifies every single species record
                         continuous_data = "9999-01-01",
                         
                         # Sampling the temporal data accordingly?
                         bk_sampling = TRUE, # [[LOGICAL]] Should a background sampling be implemented
                         sampling_type = "Random", # [[CHARACTER]] What type of background sampling needs to be implemented:(options)Random,BwData,BwData_inv,and EnvBK. For more information see the `backgroundPOINTS` function documentation
                         sampling_area, # [[SF-Spatial]] sf spatial object delimiting the sampling area for the backgroud points
                         number_points = 10000, # [[Numerica]] Number of ramdom points to sample at each time jump. This will also be the final number of Bk points returned by the function (a random sample of n bk_points will be extracted)
                         coord.xy=c("X","Y") # [[Internal]] Internal parameter for the transformation of extracted data
){
  
  # 0. Load the needed packages and functions----
  list.of.packages<-c("dplyr","sf","terra","lubridate","data.table")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  rm(list.of.packages,new.packages)
  
  # a. Get the initial details----
  start_date <- time(x) %>% min()
  end_date <- time(x[[time(x)!=continuous_data]]) %>% max()
  
  # b. Configure the point data----   
  y.index <- y[,colnames(y) %in% y.t] %>% st_drop_geometry() %>% is.na()
  y <- y[!y.index,]
  
  # check the date format, it should be y-m-d format
  y.index <- lapply(y[,colnames(y) %in% y.t] %>% st_drop_geometry() %>% unlist(),function(w) strsplit(w,split="-") %>% unlist() %>% length())
  y.index <- do.call("c",y.index)
  y <- y[y.index == 3,] ; rm(y.index)  
  
  time.points <- y[,colnames(y) %in% y.t] %>% st_drop_geometry() %>% unlist() %>% unname() %>% as.Date()
  
  # c. Subset the data based on the type and interval specified
  t.min <- time.points %>% min() 
  t.max <- time.points %>% max()
  
  t.seq <- seq(t.min,t.max,by=paste(paste0("+",jump.i),time.i))
  
  # Set up the final stage of the time sequence
  if(time.i=="year") t.max2 <- (t.seq %>% max()) %m+% years(jump.i)
  if(time.i=="month") t.max2 <- (t.seq %>% max()) %m+% months(jump.i)
  if(time.i=="day") t.max2 <- (t.seq %>% max()) %m+% days(jump.i)
  
  if(length(time.points)==0){
    stop("No valid `time` in `y` data. Format should be years-month-day")
  }
  
  # Extract the values of the raster's for the spatial points  
  r.values <- list()  
  
  for(i in 1:length(t.seq)){
    
    # Select the observations based on the time period
    if(i == length(t.seq)){
      y.index <- time.points >= t.seq[i] & time.points < t.max2 
      
    }else{
      y.index <- time.points >= t.seq[i] & time.points < t.seq[i+1]
    }
    
    # Jump to the next case if there is no information on the points or the raster data
    if(sum(y.index,na.rm=T)==0) next()
    #if(!t.seq[i] %in% time(x)) next()
    
    # Subset the points 
    r.y <- y[y.index,]
    
    # Subset the raster
    t.index<-time(x) >= t.seq[i] & time(x) < t.seq[i+1]
    if(sum(t.index[!is.na(t.index)])==0){next} 
    
    if(i == length(t.seq)){
      r.x <- x[[time(x) >= t.seq[i] & time(x) < t.max2]]
      
    }else{
      r.x <- x[[time(x) >= t.seq[i] & time(x) < t.seq[i+1]]]
    }  
    
    # If the jump paramter is greather than 1 we need to aggregate the raster information----
    if(jump.i > 1){
      
      col_groups <- split(seq_along(names(r.x)), names(r.x))
      
      # Apply the mean to each group
      aggregated_df <- lapply(col_groups, function(idc) {
        if (length(idc) > 1) {
          terra::mean(r.x[[idc]], na.rm = TRUE)
        } else {
          r.x[[idc]]
        }
      })
      r.xp <- rast(aggregated_df)
      
    }else{
      r.xp <- r.x
    }
    
    
    # Extract the information          
    if(bk_sampling==TRUE){
      
      # Generate the background points
      if(sampling_type %in% c("Random","BwData","BwData_inv")){
        bk_points <- backgroundPOINTS(presence = r.y,
                                      background_n = number_points,
                                      TrainTest = 1,
                                      range_samp = sampling_area,
                                      weights.p = sampling_type)$Train
      }
      
      # For the environmental background type
      if(sampling_type %in% "EnvBK"){
        # The selection of variables can have an impact on the shape and distribution of the PCA scores, therefore we are going to 
        if(nrow(r.y)<6){
          warning("EnvBK: Minimun numbers os presence points lower than threshold for EnvBK calculations. Sampling bk points at random")
          bk_points <- backgroundPOINTS(presence = r.y,
                                        background_n = number_points,
                                        TrainTest = 1,
                                        range_samp = sampling_area,
                                        weights.p = "Random")$Train
          
        }else{
          bk_points <-  bk_env(p.points=r.y,
                               env_var=r.xp,
                               n_bk=number_points,
                               density_pc=TRUE)[["points"]] %>% st_cast("POINT")
        }
      }
 
      #
      env_bk <- terra::extract(x=r.xp,y=bk_points %>% vect()) %>% cbind(st_coordinates(bk_points))
      env_bk <- env_bk %>% dplyr::mutate(ID.x=paste(i,1:nrow(env_bk),sep="-"),
                                         time_step=paste(t.seq[i],t.seq[i+1],sep=":"),
                                         presence=0,
                                         .before=1)
      
      env_vals <- terra::extract(x=r.xp,y=r.y %>% vect()) %>% cbind(st_coordinates(r.y))
      env_vals <- env_vals %>% dplyr::mutate(ID.x=r.y[,colnames(r.y) %in% id.p] %>% st_drop_geometry() %>% unlist(),
                                             time_step=paste(t.seq[i],t.seq[i+1],sep=":"),
                                             presence=1,
                                             .before=1)
      
      env_vals <- rbind(env_vals,env_bk)
      
    }else{
      env_vals <- env_vals %>% mutate(presence=1,.before=0) %>% dplyr::mutate(coordinates=st_coordinates(bk_points))
    }
    
    # Store the data   
    r.values[[i]]<-env_vals
    
    rm(y.index,r.x,r.xp,r.y,env_vals)
    print(t.seq[i])
  }
  
  r.values <- r.values[!lapply(r.values,is.null)%>%unlist()]
  r.values <- r.values %>% data.table::rbindlist(use.names = TRUE,fill = TRUE)
  
  # Get the observations with complete variable configuration
  r.values <- r.values[r.values %>% complete.cases(),]
  
  if(nrow(r.values)<number_points){
    warning("Time_matchine - Check the data\nNUMBER OF RECORDS LOWER THAN BACKGROUND SAMPLE! Check variables spatial and temporal domains!")
  }
  
  # Sample the background data to shape the final dataset
  index_bk <- sample(x=r.values %>% filter(presence==0) %>% dplyr::select(ID.x)%>%unlist(),size=number_points)
  r.values <- r.values %>% filter(presence==1|ID.x%in%index_bk)
  
  # For the non-temporal variables
  if(continuous_data %in% time(x)){ 
    r.values <- cbind(r.values,terra::extract(x = x[[time(x) %in% continuous_data]],r.values %>% st_as_sf(coords=coord.xy) %>% vect(),ID=FALSE))
  }
  
  return(r.values)
}

# End of the Function


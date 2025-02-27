#
###//\/\/\/\/\/\/\////\/\/\/\/\/\/\//\/\/\/\//\////\/\/\////---_-__---\/\\\\\\/////\/\/\/\/\/\##-#
#                       Run the prediction and extract some testing information
###///\/\/\/\/\/\/\////\/\/\/\\/\/\/\/\/\/\/\/\/\\\\\/\/\/\//\/\/\/\/\/\/\/\___-\-\-\-\/\/\///##-#
#
predict_MAXENT<-function(x, # route to the model location
                        spp, # name for the ouput files
                        new_vals, # rast stack contaning the variables of interest
                        var_list, # List contaning the variables to use for the predictions
                        results_route= paste(getwd(),"Results_maxent",sep="/"), # Route to the folder containing the model results
                        ref_rast=NULL, # If predictions are meant to be stacked, do you want to adapt them to a reference raster?
                        Cut=TRUE, # If predictions are meant to be restricted to a known geometry, add its route here (it needs to conform to a vect/terra compatible format)
                        buffer=NULL, # Do you want to apply a buffer to increase the area over which to run the predictions? (only if cut is not NULL)
                        report=FALSE, # Creates a short report with some model information
                        mod_name # name for the output model
                        ){
  
  # Create the output folders
  dir.create(results_route,showWarnings = FALSE)
  tif_route<-paste(results_route,spp,sep="/") ; dir.create(tif_route,showWarnings = FALSE)
  #test_route<-paste(results_route,"Test",sep="/") ; dir.create(test_route,showWarnings = FALSE)
  
  # Load the model and spatial information
  load(x)
  modX<-mx_mod$Train_model  
  area_m<-mx_mod$Area_mod
  
  # Crop the prediction variables at the adequate size
  if(Cut==TRUE){
    print(paste0("Warning! Distance units are in",
        st_crs(area_m, parameters = TRUE)$units_gdal,"!!")
        )
    
  if(is.null(buffer)){
    px<-terra::crop(new_vals,area_m)
    
  }else{
    area_m<-r.buff<-sf::st_buffer(area_m,dist=buffer) #%>% sf::st_union()
    px<-terra::crop(new_vals,area_m)
  
    }
    px<-terra::mask(px,area_m %>% vect())  
    
  }
  
  # Select the variables and run the predictions
  # for(w in 1:length(var_list)){
  #       if(w==1){
  #         xy<-terra::predict(modX,px[[var_list[[w]]]])
  #         
  #           }else{
  #             py<-px[[var_list[[w]]]]
  #               names(py)<-var_list[[1]] # Set the names of the variables to be the same as in the model, needs to adjust this
  #                 xy<-terra::predict(modX,py)
  #         }
    
  
    xy<-terra::predict(modX,px[[var_list]])
  
  
    # Does the raster need to be transform?
    if(!is.null(ref_rast)){
      xy<-raster::resample(x=xy,y=ref_rast)
    }
    
      # Save the predicted Tiff into a folder
      terra::writeRaster(xy %>% mask(area_m %>% vect()),
                          filename = paste(tif_route,paste0(spp,"_",mod_name,".tif"),sep="/"),
                          overwrite=TRUE)
      
        # Export the results of the sensitivity into an RMarkdown htlm file
        # if(report==TRUE){
        # date<-Sys.Date()
        # rmarkdown::render(input = "./Base_code/Scripts/SDM_portion/ShortReport.Rmd",clean=TRUE, 
        #                   output_file = paste(test_route,spp,sep="/"))
        # }

  # Extract some metrics to run basic model filtering and selection
  mod.details<-data.frame(model=spp,
                          Presence=mx_mod$Train_model@presence %>% nrow(),
                          Absence=mx_mod$Train_model@absence %>% nrow(),
                          AUC=mx_mod$Test@auc,
                          Kappa=mx_mod$Test@cor)
  return(mod.details)
}

#
# End of the Function
#

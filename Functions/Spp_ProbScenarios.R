#/\/\/\//\/\/\/\/\/\\\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/////\/\//\/\/\/\///\/\/\/\#
# FUNCTION TO CALCULATE THE DIFFERENT SCENARIOS OF SPECIES OCCURRENCE        #
#///\/\/\/\/\/\/\\\\/\/////\/\/\/\/\/\//\/\/\\\\\\//////\/\/\/\/\/\/\/\/\/\/\#
#
#
Spp_ProbScenarios<-function(route_dis, # route to the predictions folder
                            plot=TRUE, # Should results be ploted
                            all.p=TRUE # Should all scenarios be caltulated and returned?
                            ){
  # 0. Load the packages needed ----
  list.of.packages<-c("dplyr","terra","parallel")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  rm(list.of.packages,new.packages)
  
  # 0.1 Custom functions to run the calculations----
  # a. Extract the information from the species and disease folders
  model_objs<-function(xy){
    
    # Check the disease folder and extract the information
      y_rast<-list.files(xy,pattern = ".tif$",full.names = TRUE) %>% rast() # 
      names(y_rast)<-basename(list.files(xy,pattern = ".tif$")) %>% gsub(pattern=".tif$",replacement="")
      
      return(y_rast)
    }
  
  # b. Calculate the probabilities of: a) having at least one host in a given pixel; ----
    # 1) 1 to n host for pixel; and c) all hosts
    # 2) Prob of having at least one Host per pixel (X1 or X2 or Xi...Xn)= P(X1)+P(X2)+P(Xi...Xn)-P(X1*X2) 
    # - P(X1*Xi...Xn) + (P(X1)*P(X2)...P(Xn)) 
    cp <- function(species_probs,all.p=all.p,...){ # I need to optimise the probability calculations
      
      if(all(is.na(species_probs))){
        return(NA)
      
        }else{
      # 0.1 Replace the NA values with 0
        species_probs[is.na(species_probs)]<-0
          
      # a. Get the probability of x species being present at the same time---- 
        ev <- do.call(expand.grid,replicate(length(species_probs),0:1,simplify=FALSE))
        
        pe <- apply(ev,1,function(x) prod(species_probs*(x==TRUE)+(1-species_probs)*(x==FALSE)))
        y<-tapply(pe,rowSums(ev),sum) # this is the probabylity of finding exacly 1, 2,...n species
        
          # b. To calculate the probability of at least 1 events happening;---- 
          # we need to 1-p(0)-p(1...n) the exact probabilities
          y_p<-c(y[1],1-(y[1]),y[length(y)]) # Probability of at least one event, probability of all events
          
            # b.1 Probability of more than one event happening, when y>2----
            if(length(species_probs)>2){
              for(i in 2:c(length(y)-2)){
                y_p<-c(y_p,1-(sum(y[1],y[2:i]))) # probability of at least x events
                }
          
                # c. name the probabilities  
                names(y_p)<-c("No_spp","At_least_1_spp","All_spp",
                              paste0("At_least_",2:(length(y)-2),"_spp"))
                # return(y_p)
              
              }else{
              # c. name the probabilities
              names(y_p)<-c("No_spp","At_least_1_spp","All_spp")
              return(y_p)
              }
            }
        }
  
  # 1. Check the disease folder and extract the information ----
  # For each disease, extract the host and vector predictors, if we have more than one model for species average 
  if(length(route_dis)==0){
    paste0("No prediction information for ",basename(route_dis))
    
    }else{
    if(length(route_dis)==1){
      # There is only one file per specie
      spp_ProbDist<- model_objs(route_dis) # stack the layers
      names(spp_ProbDist)<-names(spp_ProbDist)
      
      # There is only 1 species
      if(terra::nlyr(spp_ProbDist)==1){
        # The probability is equal to the probability of the species
        prob_spp<-spp_ProbDist
        
      }else{
        prob_spp<-terra::app(spp_ProbDist,fun=cp,cores=1)
      }
      
    }else{
      # The probability of a host being presence is a combination of the probabilities of the different species being or not present in each
      # pixel.
      spp_ProbDist<-lapply(route_dis,function(x) model_objs(x) %>% mean()) %>% rast() # stack the layers
      prob_spp<-terra::app(spp_ProbDist,fun=cp)
      
      }
    }
  
  # 3. Display the information ----
  if(plot==TRUE){
    plot(prob_spp,main=names(prob_spp))
    }
  
  # 4. Return the scenarios ----
    return(prob_spp)
}

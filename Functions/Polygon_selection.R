#
# Area function: select the subpolygons of a spatial polygon data.frame by area and returns 
# spatialpolygon dataframe withouth the polygons that falls below an area threshold
#
select_area<-function(x, # Spatial polygon data.frame with the areas we want to filter/select
                      area_p # Area threshold in m2
                      ){

  # 1. Add the polygon ID to the SpatialPolygonDataset
  x@data$ID<-lapply(x@polygons,function(p) p@ID) %>% unlist()
  
  # 2. Extract and process the SpatialPolygons from the SpatialPlygonsDataFrame----
  y<-polygons(x)
  
  calc_area<-function(j,ref_x=x){ # returns the area in m2
    
    if(class(j)=="Polygons"){
      j %>% list() %>% SpatialPolygons(proj4string=CRS(proj4string(ref_x))) %>%
        raster::area()  
        
    }else{
        j %>% list() %>% Polygons(ID="1") %>% list() %>% 
        SpatialPolygons(proj4string=CRS(proj4string(ref_x))) %>%
        raster::area()
    }
  }
  
  to_remove<-NULL
  
  for(i in 1:length(y@polygons)){
     if(length(y@polygons[[i]]@Polygons)>1){ # if we have multiple polygons
        # area_pols<-lapply(y@polygons[[i]]@Polygons,function(x) slot(x,"area")) %>% unlist()   
          
       area_pols<-lapply(y@polygons[[i]]@Polygons,calc_area) %>% unlist() # Calculate the areas
       
     if(all(area_pols>area_p)){# Are all the polygons greather than the treshold?
          next()
        }else{
          
   # Are all polygons below the treshold?
     if(!TRUE %in% c(area_pols>area_p)){ # Then remove the whole SpatialPolygon
        # y@polygons[[i]]<-NULL
        # y@plotOrder[[i]]<-NULL
       to_remove<-c(i,to_remove)
        
        }else{ 
        # Remove the individual subpolygons
          y@polygons[[i]]@Polygons<-y@polygons[[i]]@Polygons[c(area_pols>area_p)] # select the polygons over the area threshold
          y@polygons[[i]]@plotOrder<-1:length(y@polygons[[i]]@Polygons) # remove the plot_order of removed polygons
          
        # Remove the comments of the discarded polygons
          comment_p<-strsplit(slot(y@polygons[[i]],"comment"),split = " ") %>% unlist()
          comment_p<-comment_p[c(area_pols>area_p)]
            
          comment(y@polygons[[i]])<-paste(comment_p,collapse=" ")
        }
          
      }} else {
        # area_pols<-y@polygons[[i]] %>% slot(name="area")
        area_pols<-y@polygons[[i]] %>% calc_area()# Calculate the areas
        
        if(area_pols<area_p){
          # y@polygons[i]<-NULL
          # y@plotOrder[i]<-NULL
          to_remove<-c(i,to_remove)
          }
        }  
      }
  
  # 3. Remove the polygons that don't meet the threshold
  # Are they polygons to remove?
  if(is.null(to_remove)){
    return(x)
  }else{
    
    # Are they polygons left?
    if(length(to_remove)==length(y@polygons)){
      
      return(NULL)
      
      }else{
  
  # Remove the polygons that doesn't pass the cut    
  y@polygons<-y@polygons[-c(to_remove)]
  
  #comment_list<-sapply(slot(y,"polygons"), function(pp) slot(pp, "plotOrder"))
  y@plotOrder<-y@plotOrder[-to_remove]
  
  # 4. Remove the NA polygons from the SpatialPolygons object 
  y@polygons<-y@polygons[!is.na(y@polygons)]
  y@plotOrder<-c(1:length(y))
  
  #comment_list<-sapply(slot(y,"polygons"), function(pp) slot(pp, "comment"))

  # 5. Add the comments and the plot order into the remaining polygons
  # for(l in 1:length(y)){
  #     comment(y@polygons[[l]])<-comment_list[[l]]
  #     y@polygons[[l]]@plotOrder<-l
  #     }
  
  # 6. Add the data.frame of the remaining polygons into the filtered object
  pid <- sapply(slot(y, "polygons"), function(x) slot(x, "ID"))
  
  pols_df<-x@data[x@data$ID %in% pid,]
  row.names(pols_df)<-pols_df$ID
  
  y<-sp::SpatialPolygonsDataFrame(y,pols_df,match.ID = TRUE)
  
  return(y)
      }
    }
  }


# check some things
# Ahh<-function(p){
# 
#   ww<-list()
#   for(k in 1:length(p@polygons)){
#   x<-p@polygons[[k]]@Polygons %>% length()
#   y<-slot(p@polygons[[k]],"comment") %>% strsplit(split=" ") %>% unlist() %>% length()
#   w<-slot(p@polygons[[k]],"plotOrder") %>% length()
#   
#   ww[[k]]<-data.frame(x=x,y=y,w=w,equal=c(x==y & y==w & w==x))
#   }
#     return(ww %>% rbindlist())
# 
# }
# 
# p<-Ahh(y)

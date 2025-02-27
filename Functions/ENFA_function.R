###########################################################################.
#  Colection of Functions to Run &                                        #
#                                 Plot the results from an ENFA analysis  #
###########################################################################.
#
#
# Main Function to run the ENFA analysis
# Adapted from the ENFA ade4 package function - see https://github.com/cran/adehabitatHS/blob/master/R/enfa.r for the original function
#
ENFA_function<-function(data, # Data.frame containing the environmental information with no NAs
                        presence_index, # an index of length equal to nrow(data) showing the presence values for each row (0=absence; 1-n= number of presence)
                        n_speciation_axes=2, # number of speciation axes to preserve, it should be in the range of 1 to ncol(data) (default=1)
                        row.weights= rep(1, nrow(data))/nrow(data), # weights for the rows of the data.frame
                        col.weights= rep(1, ncol(data))# weight for the variables of the data.frame
                        ){
  
  # 0. Load the needed packages
  list.of.packages<-c("tidyr","terra","sf","data.table","dplyr")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  rm(list.of.packages,new.packages)
  
  ## Adapted from the ENFA ade4 package function - see https://github.com/cran/adehabitatHS/blob/master/R/enfa.r for the original function
  # Presence vector (0,1)
  prb <- presence_index
  presence_index <- presence_index/sum(presence_index) # Proportion of presences from the total in each cell
  
  # Proportion of total row weigths; if lw=rep(1, nrow(rast_vals$tab[,-1]))/nrow(rast_vals$tab[,-1]) their sum =1
  row.weights <- row.weights/sum(row.weights)

  # Convert the data.frame into matrix and subtract the sum of weighted rows from each variable 
  Z <- as.matrix(data)
  n <- nrow(Z)
  
  f1 <- function(v) sum(v * row.weights)
  
  center <- apply(Z, 2, f1) 
  Z <- sweep(Z, 2, center)
  Ze <- sweep(Z, 2, sqrt(col.weights), "*") # multiply by the square route of the col.weights
  
  # Ze is the environmental matrix that contains the data for the whole environment or study area
  
  ## Build the data matrices for the species presence pixels/points 
  # Inertia matrices S and G
  DpZ <- apply(Ze, 2, function(x) x*presence_index) # weighted species or inertia matrix
  
  ## Marginality computation
  # Marginality vector coordinates
  mar <- apply(Z,2,function(x) sum(x*presence_index)) # Calculate the marginality over the Row 
                                                      # weighted presence matrix. If weights are equal 
                                                      # to 1 is the same as with the inertia matrix
  # Weighted marginality vector coordinates
  me <- mar*sqrt(col.weights)
  
  # Crossproduct of the environmental and inertia matrizes
  Se <- crossprod(Ze, DpZ) 
  Ge <- crossprod(Ze, apply(Ze,2,function(x) x*row.weights))
  
  ## Computation of S^(-1/2)
  eS <- eigen(Se)
  kee <- (eS$values > 1e-9)     ## keep only the positive values
  S12 <- eS$vectors[,kee] %*% diag(eS$values[kee]^(-0.5)) %*% t(eS$vectors[,kee])
  
  ## Passage to the third problem
  W <- S12 %*% Ge %*% S12
  x <- S12 %*% me
  b <- x / sqrt(sum(x^2))
  
  ## Eigenstructure of H
  H <- (diag(ncol(Ze)) - b%*%t(b)) %*% W %*% (diag(ncol(Ze)) - b%*%t(b))
  s <- eigen(H)$values[-ncol(Z)]
  # barplot(s)
  
  ## Number of eigenvalues
  
  if (n_speciation_axes <= 0 | n_speciation_axes > (ncol(Ze) - 1)){
    n_speciation_axes <- 1}
  
  ## coordinates of the columns on the specialization axes
  co <- matrix(nrow = ncol(Z), ncol = n_speciation_axes + 1)
  tt <- data.frame((S12 %*% eigen(H)$vectors)[, 1:n_speciation_axes])
  ww <- apply(tt, 2, function(x) x/sqrt(col.weights))
  norw <- sqrt(diag(t(as.matrix(tt))%*%as.matrix(tt)))
  co[, 2:(n_speciation_axes + 1)] <- sweep(ww, 2, norw, "/")
  
  ## coordinates of the columns on the marginality axis
  m <- me/sqrt(col.weights)
  co[, 1] <- m/sqrt(sum(m^2))
  
  ## marginality
  m <- sum(m^2)
  
  ## Coordinates of the rows on these axes (values of marginality and specificity for the environmental matrix)
  li <- Z %*% apply(co, 2, function(x) x*col.weights)
  
  ## Output
  co <- as.data.frame(co) # vector coordinates
  li <- as.data.frame(li) # row.values
  names(co) <- c("Marginality", paste("Specialization", (1:n_speciation_axes), sep = ""))
  row.names(co) <- dimnames(data)[[2]]
  names(li) <- c("Marginality", paste("Specialization", (1:n_speciation_axes), sep = ""))

  ## Calculate the environmental suitability (Mahalanobis distance)
  Zli <- li[, 1:(n_speciation_axes + 1)]
  f1 <- function(x) rep(x, prb)
  
  Sli <- apply(Zli, 2, f1)
  m_cord <- apply(Sli, 2, mean) # Specificity and marginality coordinates
  cov <- t(as.matrix(Sli)) %*% as.matrix(Sli)/nrow(Sli)
  maha <- data.frame(Mahalanobis.Dist = mahalanobis(Zli, center = m_cord, cov = cov))
  
  return(list(marginality=m,
              row.weights=row.weights,
              col.weights=col.weights,
              presence_index=presence_index,
              marginality_coordiantes=me,
              coordinates_axis=co,
              marginality_specificity_vals=li,
              niche_centroid_coordinates=m_cord,
              prediction=maha))
  }

###########################
# Complementary functions
###########################
#
# Plot enfa results
plot_enfa<-function(mar, # Marginality vector
                    spc, # Specialization vector
                    m, # Niche centroid
                    sp_rec, # Species records index
                    plot_sp=TRUE, # should we plot the results
                    pts=FALSE, # Should all points be plotted
                    col_p=c("grey","cyan4") # colors for the environmental and spp records/points
){
  
  # Load the needed packages
  list.of.packages<-c("tidyr","terra","sf","data.table","dplyr")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  rm(list.of.packages,new.packages)
  
  # Calculate environmental and species variability area (convex hull over the hyperspace)----
  points_index<-sp_rec>0
  
  hull_complete<-chull(matrix(c(mar,spc),ncol=2))
  X_sp <-matrix(c(mar,spc),ncol=2)[points_index,]
  
  hull_sp<-chull(X_sp)
  
  coord_env<-matrix(c(mar,spc),ncol=2)[c(hull_complete,hull_complete[1]),]       
  coord_sp<-X_sp[c(hull_sp,hull_sp[1]),]
  
  env.pol<-st_polygon(list(coord_env))
  spp.pol<-st_polygon(list(coord_sp))
  
  ENFA_pol <- st_sf(data.frame(ID=1:2,Type=c("Env","Spp"),geom=st_sfc(env.pol,spp.pol)))
  
  # Some metrics----
  ENFA_pol$Area <- st_area(ENFA_pol)
  #ENFA_pol %>% st_centroid() %>% plot(add=T)
  
  prop_spp_area <- (ENFA_pol$Area[2]/ENFA_pol$Area[1])*100 
  
  if(plot_sp){
    # Build the base plot----
    par(mar=c(5,5,4,4))
    
    plot(mar[!points_index],spc[!points_index],
         pch=19,col=NA,bty="n",
         # col=ifelse(pts,col_p[1] %>% adjustcolor(alpha.f = 0.3),NA),
         # col="grey85" %>% adjustcolor(alpha.f = 0.3),
         xlab="Marginality",ylab="Specificity")
    
    axis(1,col="white") ; axis(2,col="white")

    
    # Add the points
    if(pts){
      points(mar[!points_index],spc[!points_index],pch=19,
             col=col_p[1] %>% adjustcolor(alpha.f = 0.5))
      
      points(mar[points_index],spc[points_index],pch=19,
             col=col_p[2] %>% adjustcolor(alpha.f = 0.75))
    }
    
    # Add the lines
    plot(ENFA_pol %>% st_geometry(),add=TRUE,col=col_p %>% adjustcolor(alpha.f = 0.25),border=col_p,lwd=2)
    
    # Add the Niche centroid
    points(x=m[1],y=m[2],cex=2,col="black",bg="white",pch=21)
    abline(h=0,v=0,lty=3)
    
    # Add the legend
    legend("topright",legend=c("Environment","Species domain","Species niche centroid"),cex=0.8,

           pch=c(15,15,21),bg=c(NA,NA,"white"),col=c(col_p,"black"),pt.cex=c(1,1,1.5),bty="n")
  }
  
  return(list(centroid=m,
              areas=ENFA_pol,
              prop.overlap=prop_spp_area))
  }

#######################################
# Function to evaluate ENFA results   #----
#######################################
#
performance_ENFA<-function(EnFa_p, # model to test
                            index_vals, # index of presence and absence data
                            threshold_dist = 100 # number of classes at which to evaluate the distance accuracy
){
  
  y.acc<-data.frame(test=NA)
  threshold_seq <- seq(from=range(EnFa_p)[1],to=range(EnFa_p)[2],length.out=threshold_dist)
  threshold_seq <- threshold_seq[-c(1,length(threshold_seq))]
  
  
  for(l in 1:length(threshold_seq)){
    
    # Create confusion matrix
    Conf.Mat<-table(predict=ifelse(EnFa_p >= threshold_seq[l], 1, 0),
                    real=index_vals)
    
    if(nrow(Conf.Mat)<2){
      
      A<-list(Accuracy=NA,
              TypeI=NA,
              TypeII=NA,
              TNR=NA,
              TPR=NA,
              kappa=NA)
      
      b <- A %>% unlist() %>% as.data.frame() 
      b$test<- names(A)
      
      y.acc <- merge(y.acc,b,by="test",all=TRUE) ; colnames(y.acc)[ncol(y.acc)] <- threshold_seq[l]
      next()
    }
    
    Tn<-Conf.Mat[1,1] ; Fn<-Conf.Mat[1,2]
    Fp<-Conf.Mat[2,1] ; Tp<-Conf.Mat[2,2]
    
    
    A<-list(Accuracy=(Tn+Tp)/(Tn+Tp+Fp+Fn),
            TypeI=(Fp)/(Tn+Tp+Fp+Fn),
            TypeII=(Fn)/(Tn+Tp+Fp+Fn),
            TNR=(Tn)/(Tn+Fp),
            TPR=(Tp)/(Tp+Fn),
            kappa=(2*(Tp*Tn-Fn*Fp))/((Tp+Fp)*(Fp+Tn)+(Tp+Fn)*(Fn+Tn))
            )
    
    b <- A %>% unlist() %>% as.data.frame() 
    b$test<- names(A)
    
    y.acc <- merge(y.acc,b,by="test",all=TRUE) ; colnames(y.acc)[ncol(y.acc)] <- threshold_seq[l]
    
    #rm(A,b)
  }
  
  if(length(y.acc)<2){
    print("No performance information for the Model")
    return(NULL)
  }else{
    y.acc<-y.acc[!is.na(y.acc$test),] ; rownames(y.acc)<-y.acc$test
    y.acc<-y.acc[,-1] %>% as.matrix()
    threshold_seq %in% as.numeric(colnames(y.acc))
    
    
    return(list(accuracy=y.acc,dist_vals=threshold_seq))
  }
}

# End of the function
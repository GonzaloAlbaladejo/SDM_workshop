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
##                        Species distribution modelling                              ##-#  
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

# c. Load the environmental information----
  env_lyrs <- lapply("./Data/Spatial/Processed/Raster" %>% list.files(patter=".tif$",full.names = TRUE),rast)
  env_data <- do.call("c",env_lyrs[-c(3)])
  
  plot(env_data[[c(1:3,21:23)]],axes=F)

# 1. Some more preparations----
# Species Distribution Models use a combination of presence and absence/background data to establish the relationships between
# species spatial records and environmental gradients. Once defined, these relationships would define the species potential/realize or 
# realize niche of the species and allow us to extrapolate this niche into space. The configuration of environmental variables, the 
# distribution and number of presence points and the selection of background or pseudo-absence data would have an impact on our niche 
# modelling.
# 
# a. Background and absence data
# We are going to generate 4 different types of pseudo-absence and background data for the Straw-coloured fruit bat (Eidolon helvum).
# For this we are going to use the backgroundPOINTS function. This function allow us to distribute spatial points around our study area
# using 4 different methods:
#
# a.1 Prepare the background and pseudo-absence information----
# GBIF.org (25 February 2025) GBIF Occurrence Download  https://doi.org/10.15468/dl.gvsrkj
# https://doi.org/10.15468/dl.gvsrkj (GBIF DOI)
  sp_points <- Species_list[Species_list %>% grepl(pattern="Eidolon helvum")]
  
  spp <- sp_points %>% basename() %>% gsub(pattern=".csv",replacement="") # load the points and extract the species name
  sp_points <- sp_points %>% read.csv() %>% st_as_sf(coords=c("X","Y"),crs=crs_p) # load and transform our points into an sf object
  
  rX <- IUCN_ranges %>% filter(BINOMIAL==spp) %>% st_union()

# For the pseudo-absence data we are going to draw a buffer of 100 km around our presence data, so presence and absence data do not overlap
  dist_buff <- 50
  units(dist_buff)<-"km"

# We are going to define our sampling area
  to_crop <- sp_points %>% st_buffer(dist = dist_buff) %>% st_union()
  sample_area <- rX %>% st_buffer(dist = dist_buff) %>% st_difference(to_crop) %>% st_union()
  sample_area <- Stdy_area %>% st_intersection(sample_area)

# a.2 Generate the pseudo-absence data
# The weights.p argument on the backgroundPOINTS function controls how points are sampled into our study area
  bk_points <- backgroundPOINTS(presence = sp_points, # [sf points] distribution points of the species
                 background_n = nrow(sp_points), # [Numeric] number of background points
                  TrainTest=1, # [Numeric] proportion of data used for the testing and training of the models
                    range_samp=sample_area, # [Character] Route to the range information of species
                      weights.p="Random")

# display the information
  plot(Stdy_area %>% st_geometry(),col="grey95",xlim = st_bbox(rX)[c(1, 3)],
       ylim = st_bbox(rX)[c(2, 4)])
  plot(rX %>% st_geometry(),add=T,col="grey")
  plot(sample_area %>% st_geometry(),border="black",lty=3,add=T)
  
  plot(sp_points %>% st_geometry(),pch=19,col="gold1",cex=1,add=T)
  plot(bk_points$Train %>% st_geometry(),pch=19,col="#556b2f",add=T)
  
  legend("left",legend=c("Spp range","Sampling area","Presence","Absence"),
         pch=15,col=c("grey","black","gold1","#556b2f"),bty="n",xpd=T)

# Lets repeat the process for the rest of sampling methods
  bk_types <- c("Random","BwData")
  pseudo_abs <- lapply(bk_types,function(x) backgroundPOINTS(presence = sp_points, 
                                background_n = nrow(sp_points), TrainTest=1,
                                range_samp=sample_area, weights.p=x))

# Set the background observations
  sample_area_bk <- rX %>% st_buffer(dist = dist_buff) %>% st_union()
  sample_area_bk <- Stdy_area %>% st_intersection(sample_area_bk)
  
  bk_points <- lapply(bk_types,function(x) backgroundPOINTS(presence = sp_points, 
                                background_n = nrow(sp_points), TrainTest=1,
                                range_samp=sample_area_bk, weights.p=x))

# Display the information
  lt <- layout(matrix(c(1,2,3,4),ncol=2,nrow=2,byrow = F))
# layout.show(lt)

  for(i in 1:length(bk_types)){
    # a. Psudo-absences
    plot(Stdy_area %>% st_geometry(),col="grey95",xlim = st_bbox(rX)[c(1, 3)],
          ylim = st_bbox(rX)[c(2, 4)])
    
    plot(rX %>% st_geometry(),add=T,col="grey")
    plot(sample_area %>% st_geometry(),border="black",lty=3,add=T)
    
    plot(sp_points%>%st_geometry(), pch=19, col="gold1", cex=0.5,add=T)
    plot(pseudo_abs[[i]]$Train %>%st_geometry(), col="#556b2f", cex=0.25,add=T)
    
    mtext(bk_types[i],side=3,adj=0)
    legend("bottomleft",legend=c("Spp range","Sampling area","Presence","Absence"),
           pch=15,col=c("grey","black","gold1","#556b2f"),bty="n",xpd=T,cex=0.8)
    
    
    # b. bk_points
    plot(Stdy_area %>% st_geometry(),col="grey95",xlim = st_bbox(rX)[c(1, 3)],
       ylim = st_bbox(rX)[c(2, 4)])
    plot(rX %>% st_geometry(),add=T,col="grey")
    plot(sample_area_bk %>% st_geometry(),border="black",lty=3, add=T)
  
    plot(sp_points %>% st_geometry(),pch=19,col="gold1",cex=0.5, add=T)
    plot(bk_points[[i]]$Train %>% st_geometry(),col="tomato",cex=0.25, add=T)
  
    mtext(bk_types[i],side=3,adj=0)
    legend("bottomleft",legend=c("Spp range","Sampling area","Presence","Absence"),
         pch=15,col=c("grey","black","gold1","tomato"),bty="n",xpd=T,cex=0.8)
    }

# a.2 Extract the data from the environmental variables ----
  p_data <- env_data %>% terra::extract(y=sp_points %>% vect(),ID=F)

# add the type of record (1=presence, 0= psudo-absence/background)
  p_data <- p_data %>% mutate(record=1,.before=1) 

# Repeat the process for the background and psudo-absence data
  pseudo_abs <- lapply(pseudo_abs,function(x) env_data %>% terra::extract(y=x$Train %>% vect(),ID=F) %>% mutate(record=0,.before=1))
  names(pseudo_abs) <- bk_types
  
  bk_data <- lapply(bk_points,function(x) env_data %>% terra::extract(y=x$Train %>% vect(),ID=F) %>% mutate(record=0,.before=1))
  names(bk_data) <- bk_types

# We now have all the information for modelling the niche of our species and predict their distribution
#  
#/\/\/\//\/\/\/\/\/\/\/\//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\  
# 2. Niche modelling ----
#/\/\/\//\/\/\/\/\/\/\/\//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\  
#
# Species distribution models (SDM's) methods can be classified into three theoretical groups. The first one deals with situations in
# which only presence data is available for our study. This are usually termed presence only methods and are mostly use to create basic 
# species suitability maps, basic niche extrapolations, and species niche descriptions. Similar to presence only methods, presence 
# background methods use presence data with a combination of background or pseudo-absence points. These background points are meant to 
# offer a contrast against to which compare the presence data. Neither the background or pseudo-absence data reflect the real distribution 
# of the species but an attempt to designate areas with a low probability of being occupied by our species.   
#
# 2.1 Environmental Niche Factor Analysis (ENFA)----
# ENFA is a presence only method that analyze the multidimensional space form by the environmental information and the presence data to calculate
# the centroid of the species niche, the distance between that centroid and the rest of the environment (marginality), and the suitability
# of the environment
#
  ind_x <- env_data %>% crop(rX %>% vect) # We are only interested in the species range

# Rather than working directly with the raster data we are working with the values of the raster and the location of presence 
# data within that index
  r.dat <- rast_to_vect(ind_x) # Extract the information from the raster
  n.row.dat <- prod(r.dat[["dim"]])
  
  pres_index<-rep(0,times=n.row.dat)
  
  obs <- terra::extract(x=ind_x, y=sp_points %>% vect(), cells=T)$cell # Extract the index (cell id's) for the presence data
  pres_index[obs]<-1
  
  obs_index <- pres_index[-r.dat$index_missin]

# a. Run the ENFA analysis----
  ENFA.r <- ENFA_function(data = r.dat$tab[,!colnames(r.dat$tab) %in% "cell"], # Data.frame containing the environmental information with no NAs
                          presence_index = obs_index)

# Transform the ENFA_results into rasters for the export
  empty_rast <- ind_x[[1]] ; empty_rast[-is.na(empty_rast)] <- NA

# Get the different ENFA values
  maha <- empty_rast ; maha[r.dat$tab$cell] <- ENFA.r$prediction
  Marginality <- empty_rast ; Marginality[r.dat$tab$cell] <- ENFA.r$marginality_specificity_vals$Marginality
  Specialization <- empty_rast ; Specialization[r.dat$tab$cell] <- ENFA.r$marginality_specificity_vals$Specialization1

  ENFA_rast <- c(maha,Marginality,Specialization) ; names(ENFA_rast) <- c("Mahalanobis_dist","Marginality","Specificity")

# b. Get the rest of the ENFA parameters----
  par(mfrow=c(1,2))
  ENFA_extra<- plot_enfa(mar=ENFA.r$marginality_specificity_vals$Marginality, # Marginality vector
                       spc=ENFA.r$marginality_specificity_vals$Specialization1, # Specialization vector
                       m=ENFA.r$niche_centroid_coordinates, # Niche centroid
                       sp_rec=obs_index, # Species records index
                       plot_sp=T, # should we plot the results
                       pts=FALSE)
  
  plot(ENFA_rast$Mahalanobis_dist %>% log1p(),col=viridis::inferno(n=200) %>% rev(),axes=F,main="Distance to\nniche centroid")
    
# c. Exploratory analysis of the species niche
  ENFA.r$niche_centroid_coordinates
  ENFA.r$marginality_coordiantes # Marginality coordinates

#  
# 2.2 Regression base approaches - Generalized Linear Models #### ----
# 2.2.1 Generalized Linear Models (GLM)----
  
  # Let's split the data into train and test datasets (70-30% respectively)
  n.test <- (nrow(p_data)*0.3) %>% round(digits=0)
  i.test <- sample(1:nrow(p_data),size=n.test)
  
  data.train <- rbind(p_data[-c(i.test),],pseudo_abs$Random[-c(i.test),]) 
  data.test <- rbind(p_data[i.test,],pseudo_abs$Random[i.test,])
  
  data.train.b <- rbind(p_data[-c(i.test),],bk_data$Random[-c(i.test),]) 
  data.test.b <- rbind(p_data[i.test,],bk_data$Random[i.test,])
  
# a. Formulation ----
  # Saturated model
  sat.form <- paste0("record ~ ",paste(names(data.train)[-1],collapse = " + ")) %>% as.formula()
  
  glm.1 <- glm(sat.form,family = binomial,data=data.train)  
  summary(glm.1)
  
# b. How does the model performs ----
  par(mfrow=c(1,2))
  # McFadden's pseudo r-squared
  glm.R2.1<-with(summary(glm.1), 1 - deviance/null.deviance) # The adjustment is not great
  
  # Boice index
  cbi <- boyce_index(pred = terra::predict(env_data,glm.1,type="response"),obs=sp_points[i.test,] %>% st_coordinates())
  
  # Model performance
  y.glm <- performance_model(glm.1,dat=data.test[,-1],index_vals = data.test$record,threshold_seq = seq(0.01,0.99,by=0.01))
  TSS.threshold <- seq(0.01,0.99,by=0.01)[which.max(y.glm[rownames(y.glm) %in% "TSS",])]

# c. Check the model ----  
  y.glm <- cbind(Threshold=colnames(y.glm)%>%as.numeric(),y.glm %>% t()) %>% as.data.frame()
  
  plot(x=y.glm$Threshold,y=y.glm$Accuracy,
       xlab="Probability threshold",ylab="Metric Value",ylim=c(0,1),
       xlim=c(0,1),axes=F,pch=NA)
  
  axis(1,labels=seq(0,1,by=0.1),at=seq(0,1,by=0.1)) ; axis(2,las=2,labels=seq(0,1,by=0.1),at=seq(0,1,by=0.1))
  col.r <- viridis::viridis(n=ncol(y.glm))
  
  for(f in 2:ncol(y.glm)){
    lines(y=y.glm[,f],x=y.glm$Threshold,lty=1,col=col.r[f]%>%adjustcolor(alpha.f = 0.75),xpd=T,lwd=1.2)  
   }
  
  legend("bottomleft",legend=colnames(y.glm)[-1],col=col.r,lty=1,lwd=2,horiz=F,bty="o",xpd=T,ncol = 2,cex=0.5)
  abline(v=TSS.threshold,lty=3)
  mtext(side=3,adj=0,"Model performance")
  
  mtext(side=3,adj=1,paste0("Boyce_index=",cbi$Boyce%>%round(digits=2),
                            "\n","r-squ=",glm.R2.1%>%round(digits=2)),col="grey22",cex=0.8)
  
 # d. Model reduction and optimization ----
  # d.1 AIC step selection ----
  step.glm1<-MASS::stepAIC(glm.1,direction=c("both"))
  step.glm1$model # Returns the dataset with the final set of variables
  
  glm.2 <- glm(record~.,data=step.glm1$model,family = binomial)
  summary(glm.2)  
  
  # d.2 Model performance---- 
  # McFadden's pseudo r-squared
  glm.R2.2<-with(summary(glm.2), 1 - deviance/null.deviance) # The adjustment is not great
  
  # Boice index
  cbi.2 <- boyce_index(pred = terra::predict(env_data,glm.2,type="response"),obs=sp_points[i.test,] %>% st_coordinates())
  
  # Model performance
  y.glm2 <- performance_model(glm.2,dat=data.test[,-1],index_vals = data.test$record,threshold_seq = seq(0.01,0.99,by=0.01))
  TSS.threshold.2 <- seq(0.01,0.99,by=0.01)[which.max(y.glm2[rownames(y.glm2) %in% "TSS",])]
  
  # d.3 Check the model ----  
  y.glm2 <- cbind(Threshold=colnames(y.glm2)%>%as.numeric(),y.glm2 %>% t()) %>% as.data.frame()
  
  plot(x=y.glm2$Threshold,y=y.glm2$Accuracy,
       xlab="Probability threshold",ylab="Metric Value",ylim=c(0,1),
       xlim=c(0,1),axes=F,pch=NA)
  
  axis(1,labels=seq(0,1,by=0.1),at=seq(0,1,by=0.1)) ; axis(2,las=2,labels=seq(0,1,by=0.1),at=seq(0,1,by=0.1))
  col.r <- viridis::viridis(n=ncol(y.glm2))
  
  for(f in 2:ncol(y.glm2)){
    lines(y=y.glm2[,f],x=y.glm2$Threshold,lty=1,col=col.r[f]%>%adjustcolor(alpha.f = 0.75),xpd=T,lwd=1.2)  
  }
  
  legend("bottomleft",legend=colnames(y.glm2)[-1],col=col.r,lty=1,lwd=2,horiz=F,bty="o",xpd=T,ncol = 2,cex=0.5)
  abline(v=TSS.threshold.2,lty=3)
  mtext(side=3,adj=0,"Model performance")
  
  mtext(side=3,adj=1,paste0("Boyce_index=",cbi.2$Boyce%>%round(digits=2),
                            "\n","r-squ=",glm.R2.2%>%round(digits=2)),col="grey22",cex=0.8)
  
  # d.4 Compare the models
  glm.x.1 <- terra::predict(env_data,glm.1,type="response")
  glm.x.2 <- terra::predict(env_data,glm.2,type="response")
  
  anova(glm.1,glm.2) # There is no difference between the two models
  
  par(mfrow=c(1,1))
  panel(c(glm.x.1,glm.x.2),col=viridis::viridis(150),main=c("glm.1","glm.2"))  
  plot(glm.x.2-glm.x.1,main=("glm.2-glm.1"),col=viridis::magma(100),axes=F)
 
# On paper the GLM is doing a good job on classifying presence and absence data, however the map shows 
# some odd patterns of distribution, and clear artifacts (althoght these might be cause by the prediction 
# of the habitat of the species way far away from the study area/range).
# With more time and tweaking, GLM's GAM's and other Regression base approaches are robust tools for the 
# analysis of species niches.

# 2.3. Machine learning - Classification methods ----
# 2.3.1. Random forest ----
  # Preparing the data
  rf_data <- data.train[complete.cases(data.train),]
  rf_data.test <- data.test[complete.cases(data.test),]
  
  rf_data$record <- rf_data$record %>% as.factor() # Need to transform our response variable into a factor to force the model to run a classification
  rf_data.test <- rf_data.test$record %>% as.factor()                                                
  
  # Run a simple random forest classification model
  rf.1 <- randomForest(sat.form,data=rf_data,ntree=1000,importance=TRUE)
  rf.x.1.res <- predict(env_data,rf.1,type="response")
  rf.x.1.prob <- predict(env_data,rf.1,type="prob")  
  
# 2.3.1.a Model performance and reduction ----
  rf.1 # The confusion matrix is the most important piece of information here!
  v.imp <- importance(rf.1) # Importance of the variables
  v.imp <- v.imp[v.imp[,4] %>% order(decreasing = T),] # Lets take the first 13 variables
  
  rf_data.2 <- rf_data[,colnames(rf_data) %in% c("record",rownames(v.imp)[1:13])]
  rf.2 <- randomForest(record~.,data=rf_data.2,ntree=1000,importance=TRUE)
  rf.2 ; importance(rf.2)
  rf.2
  
  rf.x.2.res <- predict(env_data,rf.2,type="response")
  rf.x.2.prob <- predict(env_data,rf.2,type="prob")  
  
  # Display the models
  panel(c(rf.x.1.res,rf.x.2.res),col=c("skyblue","tomato"),main=c("Rf.1","Rf.2"))
  panel(c(rf.x.1.prob$X1,rf.x.2.prob$X1),col=viridis::viridis(100),main=c("Rf.1","Rf.2"))
  app(catalyze(c(rf.x.1.res,rf.x.2.res)),fun=sum,na.rm=T,cores=2) %>% plot(axes=F,col=c("grey","gold","green4"))
  points(sp_points,col="tomato")
  
# RandomForest is a fast and easy way to fit species distribution models. As other SDM methods RF is sensitive to the number of trees
# we use, the number and quality of our observations and the type of background or psudoabsence data we are using.  
#
# 2.3.2. Artificial Neural Networks (ANNs)----
# a. We are going to fit a single-hidden layer neural network----
  # This first model is fitted using some random parameters of decay, interations (maxit), and number of hidden layers (size)
  set.seed(185)
  init.neural <- nnet(record~.,data=data.train,size=1, decay=0.3,maxit=200)
  
  # However these parameters can be optimized using a cross-validation procedure
  ControlCV <- trainControl(method = "cv",
                            number = 3)
  
  nnGrid <-  expand.grid(size = seq(1, 15, 1),
                         decay = seq(0, 0.5, 0.1))
  
  set.seed(185)
  nnetFit <- train(record ~ ., data = rf_data,
                   method = "nnet", maxit = 200,
                   tuneGrid = nnGrid,
                   trControl = ControlCV)
  
  print(nnetFit)
  nnet.params <- nnetFit$bestTune
  
  # Fit the final ANN with the selected parameters
  init.final <- nnet(record~.,data=data.train,size=nnet.params$size,
                     decay=nnet.params$decay,maxit=200)

# b. Check the performance of both models----
  par(mfrow=c(1,1))
  panel(c(predict(env_data,init.neural),predict(env_data,init.final)),main=c("ANN.1","ANN.2"),col=viridis::inferno(150))
  
  init.y.1 <- performance_model(init.neural,dat=data.test[,-1],index_vals = data.test$record,threshold_seq = seq(0.01,0.99,by=0.01))
  init.y.2 <- performance_model(init.final,dat=data.test[,-1],index_vals = data.test$record,threshold_seq = seq(0.01,0.99,by=0.01))
  
  plot(x=colnames(init.y.1)%>% as.numeric(),y=init.y.1[rownames(init.y.1)%in% "TPR",],
       type="l",col="gold2",ylim=c(0,1),xlim=c(0,1),lwd=2,lty=2,xlab="Probability threshold",ylab="Value")
  
  lines(x=colnames(init.y.1)%>% as.numeric(),y=init.y.1[rownames(init.y.1)%in% "TNR",],lwd=2,lty=2,col="purple3")
  
  lines(x=colnames(init.y.2)%>% as.numeric(),y=init.y.2[rownames(init.y.2)%in% "TPR",],lwd=2,lty=1,col="gold2")
  lines(x=colnames(init.y.2)%>% as.numeric(),y=init.y.2[rownames(init.y.2)%in% "TNR",],lwd=2,lty=1,col="purple3")
  
  legend("bottomleft",legend=c("TPR","TNR","init1","init2"),col=c("gold2","purple3","black","black"),lty=c(1,1,2,1))
  mtext("ANN performance",adj=0,side=3)
  
# 2.4. Maximum Entropy - MaxEnt ---- 
#' MaxEnt (Maximum Entropy) is a statistical modeling principle that selects the probability distribution with the highest entropy among all 
#' distributions that satisfy come given constraints. MaxEnt is exclusively use in SDM and niche modelling.
#' 
# By default MaxEnt requires more background data than the other methods, we are going to generate 10.000 random points
  bk_Max <- backgroundPOINTS(presence = sp_points,background_n = 10000, TrainTest=1,
                             range_samp=sample_area_bk, weights.p="Random")
  
  bk_Max <- rbind(p_data,cbind(record=0,terra::extract(env_data,bk_Max$Train %>% vect(),ID=F)))
  
  # Run the MaxEnt model
  Train.m <- dismo::maxent(x=bk_Max[,-1], p=bk_Max$record,
                           removeDuplicates=rm.dp)
  
  Train.m # We can access the model performance and other information by calling the model
  
  MaxEnt.mod <- predict(Train.m,env_data)
  Maxend.Pred <- plot(MaxEnt.mod)

# 3. Export the models
  
  
  
  
  
  
  list(ENFA=,
        GLM=glm.2,
       RandomForest=rf.2,
       ANN = init.final, 
       MaxEnt=Train.m)
  
#
# End of the script
#
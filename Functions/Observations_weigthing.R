#
###//\/\/\/\/\/\/\////\/\/\/\/\/\/\//\/\/\/\//\////\/\/\///////////////////\\\\\\##-#
#                 Create the background weigthed points for MAXENT algorithm
###///\/\/\/\/\/\/\////\/\/\/\/\/\/\\/\////\/\/\/\/\/\/\/\/\/\\\\\/\/\/\/\/\/\///##-#
#
# This functions generates background points and weighting vectors according to different methods  
# a) Background weighting generation based on a density kernell created from the presence data
# b) Weighting vector based on a distance difussion function based on the distribution or reference localities for the species presence data
# c) 
#

function(x,method=c("BackGround","Distance")){
  
  # Building a kernel density raster of the presence points used to create the background weighting data:
library(MASS)
coord <- occ_train@coords
k = kde2d(x=coord[,1], y=coord[,2], h=3, n=c(1760, 2405), lims = c(range(39.95569:60.00046), range(28.37629:43.04296)))
r1 = raster(k)
r.agg <- rasterize(agg.bound, r1)
rr <- mask(r1, r.agg)
plot(rr)

# Reclassify kernel density raster to 20 classes:
library(classInt)
class <- classIntervals(na.omit(sampleRegular(rr, 100000)), n=20, style="fisher")
rec1 <- reclassify(rr, class$brks)
plot(rec1, col=colorRampPalette(c("grey", "yellow", "red", "black"))(20))

# Generating 10000 background points randomly considering the probability distribution of the density raster
set.seed(123) 
bg.wt <- randomPoints(rec1, 10000, p=occ_train, prob=T)
bg.wt <- data.frame(bg.wt)
colnames(bg.wt) <- c("Longitude", "Latitude")

plot(agg.bound)
points(bg.wt, pch=19, cex=0.4, col="grey")
points(occ_train, pch=19, cex=1, col="red")

  
  
  
}


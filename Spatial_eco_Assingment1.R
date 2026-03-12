setwd("~/Desktop/GEOG71922/isashult-alt")

#Install packages & Load Libraries
install.packages(c("terra","sf","mapview","dismo","spatstat","cowplot","ggplot2","glmnet"))

library(terra) #for spatial data
library(sf) #data frames with geometry
library (mapview) #overlays an interactive arcgis map 
library(dismo) #species distribution modelling package 
library (cowplot) #plot theme for ggplot package
library (ggplot2) #data visualization package
library (glmnet) #fits lasso and regression modelling 
library (spatstat) 


##############################################################################
########## Cleaning the Data / Prep 
##############################################################################

# Loading my file into R studio
meles <- read.csv("melesmeles.csv")

# Explore Data set
head(meles)
str(meles)
dim(meles)
names(meles)

# Subset the data to only include points with complete coordinates
meles <- meles[!is.na(meles$Latitude),]

# Subset the data to only include confirmed occurrences
meles <- meles[!is.na(meles$Identification.verification.status),]

# Remove all points with uncertainty > 1000m
meles <- meles[meles$Coordinate.uncertainty_m<=1000,]


##########################################################################
########## Setting up the study extent
##########################################################################

# Create crs object
meles.latlong = data.frame(x=meles$Longitude,y=meles$Latitude)

# Use coordinates object to create our spatial points object
meles.sp = st_as_sf(meles.latlong,coords=c("x","y"),crs="epsg:4326")

# Load in polygon of the study area (central Scotland) 
scot = st_read('scotSamp.shp')

# Load in the land cover map (LCM) and clip to the polygon
LCM = rast("LCMUK.tif")

# Crop to the extent of the study area plus a little more 
LCM = crop(LCM,st_buffer(scot, dist= 1000))

# Aggregate LCM raster from 25 to 100m resolution
LCM = aggregate(LCM$LCMUK_1,fact=4,fun="modal")

# Project meles data
meles.sp=st_transform(meles.sp,crs(LCM))

# Crop the occurence points to the study area
melesFin = meles.sp[scot,]

# Mask the LCM to this boundary
LCM = crop(LCM,scot,mask=TRUE)

# Inspect
plot(LCM)
plot(melesFin$geometry,add=T)


##########################################################################
########## Setting up Land Cover Covariates
##########################################################################

# Access levels of the raster by treating them as categorical data 
LCM = as.factor(LCM$LCMUK_1)

# Create an vector object called reclass
reclass = c(0,1,rep(0,20))

# Combine with the LCM categories into a matrix of old and new values.
RCmatrix = cbind(levels(LCM)[[1]],reclass)

RCmatrix = RCmatrix[,2:3]

# Apply function to make sure new columns are numeric 
RCmatrix = apply(RCmatrix,2,FUN=as.numeric)

# Assign new values to LCM with our reclassification matrix
broadleaf = classify(LCM, RCmatrix)

#Create an vector object called reclassUrban 
reclassUrban = c(rep(0,19),1,1)

# Combine with the LCM categories into a matrix of old and new values.
RCmatrixUrban = cbind(levels(LCM)[[1]],reclass)

RCmatrixUrban = RCmatrixUrban[,2:3]

# Apply function to make sure new columns are numeric 
RCmatrixUrban = apply(RCmatrixUrban,2,FUN=as.numeric)

# Use the reclassify() function to assign new values to LCM with our reclassification matrix
urban = classify(LCM, RCmatrixUrban)


##########################################################################
########## Building a Spatial Weight's Matrix for Woodland
##########################################################################

# Our buffer analysis (earlier practical) suggests a bi modal response with a 
# strong signal at 100m radius, but an overall highest log likelihood value at 
# 900m. Assume the latter as characteristic scale of the response of meles meles 
# to broadleafed woodland for following section. We will  create a variable to 
# reflect this by generating a raster layer for which the value at each cell 
# gives the proportion of broadleaf woodland within a 900m radius. To 
# characterise the spatial distribution of this characteristic we build a spatial 
# weights matrix and entering this into a focal (neighbourhood) analysis 

# Get number of pixels needed to cover the 900m radius 
nPix = round(900/res(LCM)[1])

# Double this number (for the distance in two direction i.e. diameter) and 
# add one (for an odd numbered grid with central cell)
nPix = (nPix*2)+1

# Build weights matrix
weightsMatrix = matrix(1:nPix^2,nrow=nPix,ncol=nPix)

# Get focal cell 
x = ceiling(ncol(weightsMatrix)/2)
y = ceiling(nrow(weightsMatrix)/2)
focalCell = weightsMatrix[x,y]

# Return the index of the cell that is the focal cell in weightsMatrix
indFocal= which(weightsMatrix==focalCell,arr.ind = TRUE)

# Compute distances for populating matrix, window of analysis
distances = list()

for(i in 1:nPix^2){
  ind.i=which(weightsMatrix==i,arr.ind=T)
  diffX=abs(ind.i[1,1]-indFocal[1,1])*res(LCM)[1]
  diffY=abs(ind.i[1,2]-indFocal[1,2])*res(LCM)[1]
  
  dist.i=sqrt(diffX^2+diffY^2)
  distances[[i]]=dist.i
  
}

# Add distance values to the weights matrix
weightsMatrix[] = unlist(distances)

# Set cells outside search radius to NA
weightsMatrix[weightsMatrix>900] = NA

# Inspect
plot(rast(weightsMatrix)) 

# Normalize the weights matrix cells so distances sum to 1
weightsMatrixNorm = weightsMatrix
weightsMatrixNorm[!is.na(weightsMatrixNorm)]=1/length(weightsMatrixNorm[!is.na(weightsMatrixNorm)])

sum(weightsMatrixNorm,na.rm=T)

# Inspect
plot(rast(weightsMatrixNorm))

# # Sum neighborhood values from all surrounding cells 
# to apply our spatial weights matrix to the woodland layer
lcm_wood_900 = focal(broadleaf,w=weightsMatrixNorm,fun="sum")

plot(lcm_wood_900)


##########################################################################
########## Building a Spatial Weight's Matrix for Urban 
##########################################################################

# A seperate buffer analysis suggests a bi modal response with a 
# strong signal at 100m radius, but an overall highest log likelihood value at 
# 1500m. Assume the latter as characteristic scale of the response of meles meles 
# to urban land.

# Get number of pixels needed to cover the 1500m radius for the urban class 
nPixUrban = round(1500/res(LCM)[1])

# Next, double and add one
nPixUrban = (nPixUrban*2)+1

# Build weights matrix
weightsMatrixUrban = matrix(1:nPixUrban^2,nrow=nPixUrban,ncol=nPixUrban)

# Get focal cell 
x=ceiling(ncol(weightsMatrixUrban)/2)
y=ceiling(nrow(weightsMatrixUrban)/2)

focalCell = weightsMatrixUrban[x,y]

# Return the index of the focal cell in Urban_weightsMatrix.
indFocal= which(weightsMatrixUrban==focalCell,arr.ind = TRUE)

# Compute distances
distancesUrban = list()

for(i in 1:nPixUrban^2){
  ind.i=which(weightsMatrixUrban==i,arr.ind=T)
  diffX=abs(ind.i[1,1]-indFocal[1,1])*res(LCM)[1]
  diffY=abs(ind.i[1,2]-indFocal[1,2])*res(LCM)[1]
  
  dist.i=sqrt(diffX^2+diffY^2)
  distancesUrban[[i]]=dist.i
  
}

# Add distance values to the weights matrix
weightsMatrixUrban[] = unlist(distancesUrban)

# Set cells outside search radius to NA
weightsMatrixUrban[weightsMatrixUrban>1500] = NA

# Normalise the weights matrix 
weightsMatrixUrban[!is.na(weightsMatrixUrban)] = 1/length(weightsMatrixUrban[!is.na(weightsMatrixUrban)])

# Sum urban class from all surrounding cells
lcm_urban_1500 = focal(urban,w=weightsMatrixUrban,fun="sum")

plot(lcm_urban_1500)


#########################################################################
########## Bringing in Elevation 
##########################################################################

# Load file
demScot = rast('demScotland.tif')

# Ensure DEM file has the same resolution and origin as the other data so that 
# we can stack them together using the the resample() function. 
demScot = terra::resample(demScot,lcm_wood_900)

#inspect
plot(demScot)


#########################################################################
########## Combining the Predictor Variables, Presence Data & Coordinates
##########################################################################

# Stack the covariate file nlayers together
allEnv = c(lcm_wood_900,lcm_urban_1500,demScot)
names(allEnv) = c("broadleaf","urban","elev")

# Check 
names(allEnv)
plot(allEnv)

# Create background points
set.seed(11)

# Sample background - one point for every cell (9775)
back = spatSample(allEnv,size=2000,as.points=TRUE,method="random",na.rm=TRUE) 
back = back[!is.na(back$broadleaf),]
back = st_as_sf(back,crs="EPSG:27700")

# Get environmental covariates at presence locations
eP = terra::extract(allEnv,melesFin)

# Bind together the presence data 
Pres.cov = st_as_sf(cbind(eP,melesFin))
Pres.cov$Pres=1

# Remove the first ID field coloumn
Pres.cov = Pres.cov[,-1]

# Get coordinates for spatial cross-validation later
coordsPres = st_coordinates(Pres.cov)

# Create plain data frame version of Pres.cov for modelling
Pres.cov.df = as.data.frame(Pres.cov)[, c("broadleaf", "urban", "elev", "Pres")]

# Drop geometry column 
Back.cov = st_as_sf(data.frame(back,Pres=0))

# Get coordinates of background points for cross validation later
coordsBack = st_coordinates(back)

# Create plain data frame version of Back.cov for modelling
Back.cov.df = as.data.frame(Back.cov)[, c("broadleaf", "urban", "elev", "Pres")]

# Combine into a table
coords = data.frame(rbind(coordsPres,coordsBack))

# Assign coloumn names of DF
colnames(coords) = c("x","y")

# Stacks pres and absence data into a df
all.cov = rbind(Pres.cov, Back.cov)

# Add in coordinates for later spatial cross validation 
all.cov = cbind(all.cov,coords)

# Remove any rows with NA
all.cov = na.omit(all.cov)

# Drop the geomery coloum 
all.cov = st_drop_geometry(all.cov)


#########################################################################
########## GLM Model Construction
##########################################################################

# Specify model with all variables
glm.meles <- glm(Pres~broadleaf+urban+elev,binomial(link='logit'), data=all.cov)

# Predict 
prGLM = predict(allEnv,glm.meles,type="response")

# Plot
plot(prGLM, main='GLM, regression')

# Inspect
summary(glm.meles)

#########################################################################
########## GLM Model Evaluation
##########################################################################

# Split sample & K fold assessment of AUC metric across 5 models 
folds = 5

# Partition presence and absence data according to folds
kfold_pres <- kfold(Pres.cov, folds)
kfold_back <- kfold(Back.cov, folds)

# Create an empty list to hold our results 
eGLM <- list()
par(mfrow=c(2,3))

# Loop to iterate over folds
for(i in 1:folds){
  train <- Pres.cov.df[kfold_pres != i,] #for presence values, select all folds which are not 'i' to train the model
  test <- Pres.cov.df[kfold_pres == i,] #the remaining fold is used to test the model
  backTrain <- Back.cov.df[kfold_back != i,] #now for background values, select training folds
  backTest <- Back.cov.df[kfold_back == i,] #use the remainder for model testing
  dataTrain <- as.data.frame(rbind(train, backTrain)) #bind presence and background training data together
  dataTest <- as.data.frame(rbind(test, backTest)) #bind test data together
  glm_eval <- glm(Pres ~ broadleaf + urban + elev,
                  binomial(link="logit"),
                  data = dataTrain)   # Train GLM 
  eGLM[[i]] <- evaluate(p = dataTest[which(dataTest$Pres==1),],
                        a = dataTest[which(dataTest$Pres==0),],
                        glm_eval)   # Evaluate on held-out fold
  plot(eGLM[[i]], "ROC")
}

# Model evaluation results 
aucGLM <- sapply( eGLM, function(x){slot(x, 'auc')} )

# Calculate the mean AUC values 
mean(aucGLM)


#########################################################################
########## Expand on Variable Relationships in GLM 
##########################################################################

############ BROADLEAF response 
# Build new data frame based on mean of elev and urban but varying values for broadleaf. 
glmNewBroadleaf= data.frame(broadleaf=seq(0,max(all.cov$broadleaf),length=1000),
                  elev=mean(all.cov$elev),
                  urban=mean(all.cov$urban))


# Probability-scale predictions using type = "response" 
preds = predict(glm.meles, newdata = glmNewBroadleaf, type = "response", se.fit = TRUE)
glmNewBroadleaf$fit = preds$fit
glmNewBroadleaf$se = preds$se.fit

head(glmNewBroadleaf)

# Plot
ggplot(glmNewBroadleaf, aes(x = broadleaf, y = fit)) +
  
  geom_ribbon(data = glmNewBroadleaf, aes(y = fit, ymin = fit - 1.96 * se, ymax = fit + 1.96 * se),
              fill = "blue", alpha = 0.3) +
  geom_line(data = glmNewBroadleaf, aes(y = fit)) 


############ URBAN response
# Build new data frame based on mean of elev and broadleaf but varying values for urban. 
glmNewUrban=data.frame(urban=seq(0,max(all.cov$urban),length=1000),
                       elev=mean(all.cov$elev),
                       broadleaf=mean(all.cov$broadleaf))

# Use type = "response" for probability-scale predictions    
predUrban = predict(glm.meles, newdata = glmNewUrban, type = "response", se.fit = TRUE)
glmNewUrban$fit = predUrban$fit
glmNewUrban$se = predUrban$se.fit

# Plot
ggplot(glmNewUrban, aes(x = urban, y = fit)) +
  
  geom_ribbon(data = glmNewUrban, aes(y = fit, ymin = fit - 1.96 * se, ymax = fit + 1.96 * se),
              fill = "blue", alpha = 0.3) +
  geom_line(data = glmNewUrban, aes(y = fit)) 


######## ELEVATION response
# Build new data frame based on mean of broaleaf and urban but varying values for elev. 
glmNewElev=data.frame(elev=seq(0,max(all.cov$elev),length=1000),
                      urban=mean(all.cov$urban),
                      broadleaf=mean(all.cov$broadleaf))


# Use type = "response" for probability-scale predictions    
predElev = predict(glm.meles, newdata = glmNewElev, type = "response", se.fit = TRUE)
glmNewElev$fit = predElev$fit
glmNewElev$se = predElev$se.fit

# Plot
ggplot(glmNewElev, aes(x = elev, y = fit)) +
  
  geom_ribbon(data = glmNewElev, aes(y = fit, ymin = fit - 1.96 * se, ymax = fit + 1.96 * se),
              fill = "blue", alpha = 0.3) +
  geom_line(data = glmNewElev, aes(y = fit)) 


#########################################################################
########## Polynomial GLM Construction 
##########################################################################

# Glm model using the appropriate polynomial term for broadleaf
glm.meles.poly = glm(Pres ~ broadleaf + urban + poly(elev, 2), binomial(link = 'logit'), data = all.cov)

# Predict 
prGLMpoly = predict(allEnv,glm.meles,type="response")

# Plot
plot(prGLMpoly, main='GLM, polynomial')

# Inspect
summary(glm.meles.poly)


#########################################################################
########## Polynomial GLM Model Evaluation
##########################################################################

# Split sample & K fold assessment of AUC metric across 5 models 
folds = 5

# Partition presence and absence data according to folds
kfold_pres <- kfold(Pres.cov, folds)
kfold_back <- kfold(Back.cov, folds)

# Create an empty list to hold our results 
eGLMpoly <- list()
par(mfrow=c(2,3))

# Loop to iterate over folds
for(i in 1:folds){
  train <- Pres.cov.df[kfold_pres != i,] #for presence values, select all folds which are not 'i' to train the model
  test <- Pres.cov.df[kfold_pres == i,] #the remaining fold is used to test the model
  backTrain <- Back.cov.df[kfold_back != i,] #now for background values, select training folds
  backTest <- Back.cov.df[kfold_back == i,] #use the remainder for model testing
  dataTrain <- as.data.frame(rbind(train, backTrain)) #bind presence and background training data together
  dataTest <- as.data.frame(rbind(test, backTest)) #bind test data together
  glm_eval <- glm(Pres ~ broadleaf + urban + poly(elev, 2),
                  binomial(link="logit"),
                  data = dataTrain)   # Train GLM 
  eGLMpoly[[i]] <- evaluate(p = dataTest[which(dataTest$Pres==1),],
                        a = dataTest[which(dataTest$Pres==0),],
                        glm_eval)   # Evaluate on held-out fold
  plot(eGLMpoly[[i]], "ROC")
}

# Model evaluation results 
aucGLMpoly <- sapply( eGLMpoly, function(x){slot(x, 'auc')} )

# Calculate the mean AUC values 
mean(aucGLMpoly)


#########################################################################
########## Point Process Data Preperation
##########################################################################

# In the Spatstat package, predictor variables take the form of Im (image) objects. 
# converted from rasters. 

# Load 
library(spatstat) 

#Function to convert raster to images for spatstat. Takes one argument "im" that should be a raster object
raster.as.im = function(im) {
  r = raster::res(im)[1] #get the resolution (cell size of the raster)
  orig = ext(im)[c(1,3)]  #get the origin (bottom left corner of the raster/image)
  xx = orig[1] + seq(from=0,to=(ncol(im) - 1)*100,by=r) #set the coordinates of the columns which is just a series of number from zero increasing by 100 metres (the resolution of the raster) for every cell along the rows and columns.
  yy = orig[2] + seq(from=0,to=(nrow(im) - 1)*100,by=r) #set the coordinates of the columns
  mat=matrix(raster::values(im), ncol = ncol(im), 
             nrow = nrow(im), byrow = TRUE)[nrow(im):1, ] #now build a single matrix with the cell values and dimension we want - note that we reverse the rows in the matrix by setting nrow(im):1. This is just because spatstat reads images in reverse order to how raster object are usually organised (so the below code ensures our image is not upside-down when we come to analyse it)
  return(spatstat.geom::im(mat, 
                           xcol = xx, yrow = yy))
}

# Convert broadleaf object
broadleafIm<-raster.as.im(raster(allEnv$broadleaf))
# Convert urban object
urbanIm<-raster.as.im(raster(allEnv$urban))
# Convert elevation object
elevIm<-raster.as.im(raster(allEnv$elev))

# Inspect
names(allEnv)
print(allEnv)

# Create study window using the urbanIm image object as a template
window.poly <- as.owin(urbanIm)

# Inspect
plot(window.poly)

# Get coordinates from melesFin object for creating point pattern
melesCoords <- st_coordinates(melesFin)

# Create point pattern object 
pppMeles <- ppp(melesCoords[,1], melesCoords[,2], window = window.poly)

# Inspect object 
plot(allEnv$broadleaf)
plot(pppMeles, add=T)

# Use as.ppp() to remove points outside the window
pppMeles <- as.ppp(pppMeles)

# Re-scale from m to km
pppMeles <- rescale(pppMeles, 1000)

# Re-scale image objects to conform with the point pattern
broadleafIm <- rescale(broadleafIm, 1000)
elevIm <- rescale(elevIm, 1000)
urbanIm <- rescale(urbanIm, 1000)


#########################################################################
########## Exploratory Spatial Analysis 
##########################################################################

########## Ripley's K Test to check for clustering 
# Test using the envelope() function
Kcsr <- envelope(pppMeles, Kest, nsim=39, VARIANCE=T, nSD=2, global=TRUE)

# Inspect
plot(Kcsr, shade=c("hi","lo"), legend=T)

############ Quadrature Loop for determining how many dummy points are needed 
# for stable likelihood approximation? See where AIC change slows. 

# Find optimal quadrature scheme resolution
ndTry <- seq(100, 1000, by=100)

for(i in ndTry){
  Q.i <- quadscheme(pppMeles, method="grid", nd=i)
  fit.i <- ppm(Q.i ~ broadleafIm + elevIm + urbanIm)
  print(i)
  print(AIC(fit.i))
}

# Set quadrature scheme at optimal nd value identified above
Q <- quadscheme(pppMeles, method="grid", nd=900)

########## Rohohat Plots for looking at the raw relationship between Badger 
# intensity and each covariate to determine polynomial degree of ppm. 

# Plot 
plot(rhohat(pppMeles, broadleafIm))  
plot(rhohat(pppMeles, elevIm))       
plot(rhohat(pppMeles, urbanIm))      


#########################################################################
########## Point Process Modelling 
##########################################################################

# First PPModel: standard inhomogeneous Poisson process model, assumes spatial independence
# Construct simple model with polynomial terms assessed in rohohat plots
firstPPMod <- ppm(Q ~ poly(broadleafIm, 3) + poly(elevIm, 3) + poly(urbanIm, 2) + x + y)

# Ripley's K Test 
firstModEnv <- envelope(firstPPMod, Kest, nsim=39, VARIANCE=TRUE, nSD=2, global=TRUE)
plot(firstModEnv)

# Modified PPModel: Thomas Cluster Model
# Accounting for clustering by using a Matern process
thomasMod <- kppm(Q ~ poly(broadleafIm, 3) + poly(elevIm, 3) + poly(urbanIm, 2) + x + y, "Thomas")

# Ripley's K Test applying kinhom rather than Kest to account for inhomogeneity
thomasEnv <- envelope(thomasMod, Kinhom, nsim=39, VARIANCE=TRUE, nSD=2, global=TRUE)
plot(thomasEnv)

# Predict
prPPMod <- predict(thomasMod)

# Plot
plot(prPPMod)

# Convert plot to Raster
plot(rast(prPPMod)) 


#########################################################################
########## Evaluating the modified PPM  
##########################################################################

# here AUC = probability that a randomly-selected data point has higher predicted 
# intensity than a randomly-selected spatial location (now spatially integrated). 
# For ppm's this differs slightly from binomial models in that here it is the 
# intensity rather than the success of correctly 
# classifying (probability) a binary outcome that is being tested. 

library(spatstat)

# Evaluate 
plot(roc(thomasMod))

# Print AUC
auc.kppm(thomasMod)


#########################################################################
########## Comparing The GlM and PPM with Thomas Clustering Predictive Plots
##########################################################################

# Final Plots 
par(mfrow=c(1,2))
plot(prGLMpoly, main="Polynomial GLM prediction")
plot(rast(prPPMod), main="Thomas PPM prediction")




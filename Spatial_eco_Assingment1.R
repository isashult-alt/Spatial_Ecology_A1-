setwd("~/Desktop/GEOG71922/isashult-alt")

#Install packages & Load Libraries
install.packages(c("terra","sf","mapview","dismo","spatstat","cowplot","ggplot2","precrec","glmnet","maxnet","mlr"))

library(terra) #for spatial data
library(sf) #data frames with geometry
library (mapview) #overlays an interactive arcgis map 
library(dismo) #species distribution modelling package 
library (cowplot) #plot theme for ggplot package
library (ggplot2) #data visualization package
library (precrec) #performance analysis of binary classifications, makes ROC curves, versatile model evaluation
library (glmnet) #fits lasso and regression modelling 
library (maxnet) #maxent SDMs for glmnet
library (mlr) #machine learning



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

#Apply function to make sure new columns are numeric (here the "2" specifies that we want to apply the as.numeric function to columns, where "1" would have specified rows)
RCmatrixUrban = apply(RCmatrixUrban,2,FUN=as.numeric)

#Use the reclassify() function to asssign new values to LCM with our reclassification matrix
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

# Drop geometry column 
Back.cov = st_as_sf(data.frame(back,Pres=0))

# Get coordinates of background points for cross validation later
coordsBack = st_coordinates(back)

#Combine into a table
coords = data.frame(rbind(coordsPres,coordsBack))

# Assign coloumn names of DF
colnames(coords) = c("x","y")


#########################################################################
########## Creating Learner in the MLR Package for Binomial Model 
##########################################################################

# Combine pres and background
all.cov = rbind(Pres.cov,Back.cov)

# Add coordinates
all.cov = cbind(all.cov,coords)

# Remove any NAs
all.cov = na.omit(all.cov)

# Remove the geometry column
all.cov = st_drop_geometry(all.cov)

# Set task and make target variable categorical
task = all.cov
head(all.cov)

task$Pres=as.factor(task$Pres)

task = makeClassifTask(data = task[,c(1:4)], target = "Pres",
                       positive = "1", coordinates = task[,5:6])

# Use the make learner function to build the model approach for binomial mdel
# (logistic regression).
lrnBinomial = makeLearner("classif.binomial",
                          predict.type = "prob",
                          fix.factors.prediction = TRUE)


#########################################################################
########## Cross Validation of simple MLR Binaomial model 
##########################################################################

# The below code repeats 5-fold cross validation for both strategies, splitting 
# data into 5 chunks, training on 4, testing on 1, repeated 5 times = 25 models total.

# Set up resampling strategy: non-spatial
perf_levelCV = makeResampleDesc(method = "RepCV", predict = "test", folds = 5, reps = 5)

# Set up resampling strategy: spatial cross-validation
perf_level_spCV = makeResampleDesc(method = "SpRepCV", folds = 5, reps = 5) 

# Now we run the model re-sampling the data (task) according to the learner set up.
############## Binomial conventional cross validation 
cvBinomial = resample(learner = lrnBinomial, task =task,
                      resampling = perf_levelCV, 
                      measures = mlr::auc,
                      show.info = FALSE)

print(cvBinomial)

# Create Spatial Resampling Plots
plots = createSpatialResamplingPlots(task,resample=cvBinomial,
                                     crs=crs(allEnv),datum=crs(allEnv),color.test = "red",point.size = 1)

# Use  cowplot function to plot all folds out in a grid
cowplot::plot_grid(plotlist = plots[["Plots"]], ncol = 3, nrow = 2,
                   labels = plots[["Labels"]])


############# Binomial spatial cross validation
sp_cvBinomial = resample(learner = lrnBinomial, task =task,
                         resampling = perf_level_spCV, 
                         measures = mlr::auc,
                         show.info = FALSE)

print(sp_cvBinomial)

# Make partition plots
plotsSP = createSpatialResamplingPlots(task,resample=sp_cvBinomial,
                                       crs=crs(allEnv),datum=crs(allEnv),color.test = "red",point.size = 1)

# Use the cowplot function to plot all folds out in a grid
cowplot::plot_grid(plotlist = plotsSP[["Plots"]], ncol = 3, nrow = 2,
                   labels = plotsSP[["Labels"]])


#########################################################################
########## EXpand on GLM: Variable Relationshsips
##########################################################################

# Specify the model
glm.meles = glm(Pres~broadleaf+urban+elev,binomial(link='logit'), data=all.cov)

#Predict 
prGLM = predict(allEnv,glm.meles,type="response")

# Plot
plot(prGLM)

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

############# Specify the glm model using the appropriate polynomial term for broadleaf
glm.meles = glm(Pres ~ broadleaf + urban + poly(elev, 2), binomial(link = 'logit'), data = all.cov)

# Inspect
summary(glm.meles)


#########################################################################
##########Cross Validation of updated glm.meles model
##########################################################################

# Add squared elevation term to data to encode poly(elev,2) for MLR
all.cov$elev2 <- all.cov$elev^2

# Rebuild task with elev2 included
task.poly <- all.cov
task.poly$Pres <- as.factor(task.poly$Pres)

task.poly <- makeClassifTask(data = task.poly[, c("Pres", "broadleaf", "urban", "elev", "elev2")],
                             target = "Pres",
                             positive = "1",
                             coordinates = task.poly[, c("x", "y")])

# Learner remains the same as before
lrnBinomial = makeLearner("classif.binomial",
                          predict.type = "prob",
                          fix.factors.prediction = TRUE)

############## Conventional cross validation of updated model
cvBinomial.poly = resample(learner = lrnBinomial, task = task.poly,
                           resampling = perf_levelCV,
                           measures = mlr::auc,
                           show.info = FALSE)

print(cvBinomial.poly)

# Spatial resampling plots
plots.poly = createSpatialResamplingPlots(task.poly, resample = cvBinomial.poly,
                                          crs = crs(allEnv), datum = crs(allEnv),
                                          color.test = "red", point.size = 1)

cowplot::plot_grid(plotlist = plots.poly[["Plots"]], ncol = 3, nrow = 2,
                   labels = plots.poly[["Labels"]])

############## Spatial cross validation of updated model
sp_cvBinomial.poly = resample(learner = lrnBinomial, task = task.poly,
                              resampling = perf_level_spCV,
                              measures = mlr::auc,
                              show.info = FALSE)

print(sp_cvBinomial.poly)

# Spatial resampling plots
plotsSP.poly = createSpatialResamplingPlots(task.poly, resample = sp_cvBinomial.poly,
                                            crs = crs(allEnv), datum = crs(allEnv),
                                            color.test = "red", point.size = 1)

cowplot::plot_grid(plotlist = plotsSP.poly[["Plots"]], ncol = 3, nrow = 2,
                   labels = plotsSP.poly[["Labels"]])

# Compare AUC of simple vs polynomial model
cat("Simple Glm model - conventional CV AUC:", cvBinomial$aggr, "\n")
cat("Simple Glm model - spatial CV AUC:", sp_cvBinomial$aggr, "\n")
cat("Poly Glm model - conventional CV AUC:", cvBinomial.poly$aggr, "\n")
cat("Poly Glm model - spatial CV AUC:", sp_cvBinomial.poly$aggr, "\n")

# Both models show a drop in AUC when you move from conventional to spatial cross validation.
# Importantly,  some of the apparent predictive performance in conventional CV was inflated by spatial autocorrelation.
# The marginal improvement from adding the polynomial elevation term disappears under spatial cross validation. 
# The improvement seen in conventional CV was likely a minor overfitting artefact.
# Given that the polynomial model offers no real spatial predictive improvement and adds complexity, 
# you could make a reasonable argument for preferring the simpler linear model on the basis of parsimony. 


#########################################################################
########## Further Basic Evaluation of GLM model 
##########################################################################

install.packages("pROC") 
library(pROC) 

# Re-specify simple glm model
glm.meles = glm(Pres~broadleaf+urban+elev,binomial(link='logit'), data=all.cov)

# Creating object for ROC evaluation 
roc.meles <- roc(all.cov$Pres, fitted(glm.meles))

# Results
plot(roc.meles, col = "red")
auc(roc.meles)

# Find the threshold that maximises sensitivity + specificity
opt.threshold <- coords(roc.meles, "best", best.method = "closest.topleft")
print(opt.threshold)

# Calculate TTS value
TSS <- 0.8111753 + 0.6765 - 1
print(TSS)

# Model in oderate range (0.4–0.6). Model is better at detecting true presences 
# than correctly rejecting absences, better to slightly-overpredict than
# miss habitat, common for SDMs. 

############# Confusion matrix metrics
#install.packages("caret")
#library(caret)
# Use that threshold instead of 0.5 to classify predictions
#cpredicted.class <- ifelse(fitted(glm.meles) >= opt.threshold$threshold, 1, 0)
# Now run the confusion matrix with the optimal threshold
# confusionMatrix(as.factor(predicted.class), as.factor(all.cov$Pres), positive = "1")

#########################################################################
########## Point Process Modelling
##########################################################################
# In the Spatstat package, predictor variables take the form of Im (image) objects. 
# converted from rasters. 

#Function to convert raster to images for spatstat. Takes one argument "im" that should be a raster object
raster.as.im = function(im) {
  #get the resolution (cell size of the raster)
  r = raster::res(im)[1]
  #get the origin (bottom left corner of the raster/image)
  orig = ext(im)[c(1,3)]
  #set the coordinates of the columns which is just a series of number from zero increasing by 100 metres (the resolution of the raster) for every cell along the rows and columns.
  xx = orig[1] + seq(from=0,to=(ncol(im) - 1)*100,by=r)
  #set the coordinates of the columns
  yy = orig[2] + seq(from=0,to=(nrow(im) - 1)*100,by=r)
  
  #now build a single matrix with the cell values and dimension we want - note that we reverse the rows in the matrix by setting nrow(im):1. This is just because spatstat reads images in reverse order to how raster object are usually organised (so the below code ensures our image is not upside-down when we come to analyse it)
  mat=matrix(raster::values(im), ncol = ncol(im), 
             nrow = nrow(im), byrow = TRUE)[nrow(im):1, ]
  return(spatstat.geom::im(mat, 
                           xcol = xx, yrow = yy))
}

# Load spatstat
library(spatstat) #stat analysis of spatial point patterns and converting between raster and pixel image objects

# Convert broadleaf object
broadleafIm<-raster.as.im(raster(allEnv$broadleaf))
# Convert urban object
urbanIm<-raster.as.im(raster(allEnv$urban))
# Convert elevation object
elevIm<-raster.as.im(raster(allEnv$elev))

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

# Inspect
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
########## Ripley's K Test for pppMeles
##########################################################################
# This tests, for multiple radii, whether points within a certain radius are 
# more clustered (high values of K) or more uniform (low values of K) than would 
# be expected under complete spatial randomness.

# Test using the envelope() function
Kcsr <- envelope(pppMeles, Kest, nsim=39, VARIANCE=T, nSD=2, global=TRUE)

# Inspect
plot(Kcsr, shade=c("hi","lo"), legend=T)

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

# Plot increasing values for each environmental covariate against intensity (~density) of the point pattern
plot(rhohat(pppMeles, broadleafIm))  
plot(rhohat(pppMeles, elevIm))       
plot(rhohat(pppMeles, urbanIm))      


#########################################################################
########## Point Process Model Cont. 
##########################################################################
# Based on polynomial terms for the three covariates and the projected coordinates
firstPPMod <- ppm(Q ~ poly(broadleafIm, 3) + poly(elevIm, 3) + poly(urbanIm, 2) + x + y)

# Test whether model conforms to spatial randomness using the envelope() function
firstModEnv <- envelope(firstPPMod, Kest, nsim=39, VARIANCE=TRUE, nSD=2, global=TRUE)
plot(firstModEnv)

# Accounting for clustering by using a Matern process, apply Thomas cluster process 
thomasMod <- kppm(Q ~ poly(broadleafIm, 3) + poly(elevIm, 3) + poly(urbanIm, 2) + x + y, "Thomas")

# Create envelope based on simulations of inhomogeneous (clustering) point patterns
# Using Kinhom rather than Kest to account for inhomogeneity
thomasEnv <- envelope(thomasMod, Kinhom, nsim=39, VARIANCE=TRUE, nSD=2, global=TRUE)
plot(thomasEnv)

# Create mapped output for meles predictions
prPPMod <- predict(thomasMod)
plot(prPPMod)

# Convert plot to Raster
plot(rast(prPPMod)) 


#########################################################################
########## Evaluating the PPM with Thomas CLustering 
##########################################################################

# here AUC = probability that a randomly-selected data point has higher predicted 
# intensity than a randomly-selected spatial location (now spatially integrated). 
# for ppm's this differs slightly from binomial models (models with 2 outcomes) 
# in that here it is the intensity rather than the success of correctly 
# classifying a binary outcome that is being tested. 

# Detach pROC temporarily to avoid masking spatstat's roc function
detach("package:pROC", unload=TRUE)

# Now spatstat's roc will be used correctly
thomasROC <- roc(thomasMod)
plot(thomasROC)

# AUC
auc.kppm(thomasMod)

# AIC Comparison: Poisson vs Thomas Cluster Model
# ????
  
  
  
#########################################################################
########## Comparing The Glm and PPM with Thomas Clustering Predictive Plots
##########################################################################

# Final Plots 
par(mfrow=c(1,2))
plot(prGLM, main="GLM prediction")
plot(rast(prPPMod), main="Thomas PPM prediction")

# Normalise both to 0-1 scale for fair visual comparison
glm_norm <- (prGLM - minmax(prGLM)[1]) / (minmax(prGLM)[2] - minmax(prGLM)[1])

ppm_rast <- rast(prPPMod)
ppm_norm <- (ppm_rast - minmax(ppm_rast)[1]) / (minmax(ppm_rast)[2] - minmax(ppm_rast)[1])

# Plot side by side on same scale
par(mfrow=c(1,2))
plot(glm_norm, main="GLM prediction (normalised)", range=c(0,1))
plot(ppm_norm, main="Thomas PPM prediction (normalised)", range=c(0,1))

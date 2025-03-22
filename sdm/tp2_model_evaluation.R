# think: which data to use for training and which for validation?

library(biomod2)
library(terra)
library(ecospat)
library(PresenceAbsence)
library(ltm)
library(modEvA)
library(ROCR)
library(randomForest)
library(ggplot2)
library(caret)
library(boot)
library(pROC)
library(splines)
library(gam)
library(gtools)
library(viridis)
library(reshape2)
library(sf)
library(plyr)
library(dplyr)

data("ecospat.testData")
sp0bs <- ecospat.testData[,c("long", "lat", "Veronica_chamaedrys")]
head(sp0bs)

coord <- sp0bs[,1:2]

# project coordinates in crs of climate maps
coord<-sf_project(from=st_crs("EPSG:21781"),to= st_crs("EPSG:2056"),pts=coord)

# load and extract bioclim data
raster_files <- mixedsort(list.files("sdm/data_tp1/covariates", pattern = "bio", full.names = T))
bioclim <- rast(raster_files)

# change names of rasters
names(bioclim) = sub("_.*", "", names(bioclim))

# remove coordinates that have no climatic information
spData<-na.omit(data.frame(coord,Veronica_chamaedrys=spObs$Veronica_chamaedrys,
                           terra::extract(bioclim,coord)))


# run models for which to evaluate performance

# glm
GLM <- glm(formula = Veronica_chamaedrys ~ bio7 + I(bio7^2) + I(bio6 ^2),
           family = binomial, data = spData)
GLM.pred <- predict (GLM, newdata = spData[,c("bio7", "bio6", "bio15")],
                     type = "response")

# gam
GAM <- gam(Veronica_chamaedrys ~ s(bio7) + s(bio6) + s(bio15),
           data = spData, family = "binomial")
GAM.pred <- predict(GAM, newdata = spData[,c("bio7", "bio6", "bio15")],
                    type = "response")
# RF
rf = randomForest (x = spData[,c("bio7", "bio6", "bio15")],
                   y = as.factor(spData$Veronica_chamaedrys),
                   ntree = 1000)
rf.pred = predict(rf, type = "prob") [,2] # extract second column - probabilities

# final data set with all predictions
ObsPA <- spData$Veronica_chamaedrys
plotID <- 1:nrow(spData)
EvalData <- data.frame(cbind(plotID, ObsPA, GLM.pred, GAM.pred, RF.pred))
colnames(EvalData) <- c ("plotID", "ObsPA", "GLM", "GAM", "RF")
head(EvalData)

# 2 - measuring model accuracy
# 2.1 calibration 
par(oma = c(2, 2, 0, 0), mar = c(2, 2, 2, 1), mfrow = c(2, 2),
    cex = 0.7, cex.lab = 1.4, mgp = c(2, 0.5, 0))
for (mod in 1:3) {
  calibration.plot(EvalData, which.model = mod, color = mod + 1,
                   xlab = "", ylab = "", N.bins=10)
  mtext("Predicted Probability of presence", side = 1, line = 1,
        cex = 1.4, outer = TRUE)
  mtext("Proportion of observed presence", side = 2, line =-4,
        cex = 1.4, outer = TRUE)
}

data(SIM3DATA)
calibration.plot(SIM3DATA, which.model = 1, color = 5,
                 xlab = "", ylab = "", N.bins = 10,
                 main = "Observed vs Predicted (Bad calibration")
?calibration.plot
# useful!

# 2.1 discrimination
# contingency table for 1 glm model and threshold 0.5
table(EvalData$GLM>0.5, EvalData$ObsPA, dnn = c("Prediction", "Observation"))
accu <- presence.absence.accuracy(EvalData, which.model = 1, threshold = 11,
                                  st.dev = FALSE)
accu[,-c(1, 2)] <- signif(accu[,-c(1, 2)], digits = 2)
accu

pred.prev <- predicted.prevalence(EvalData, threshold = 11)
pred.prev[, 2:5] <- round(pred.prev[, 2:5], digits = 2)
pred.prev

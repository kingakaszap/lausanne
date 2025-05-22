# start on 06/05/2025
# sdm for current and future projections of Craterellus_tubaeformis in the Vaud alps

# Prep -----
library(biomod2) # V.4.2-4
library(mgcv)
library(ade4) # V.1.7-22
library(MASS) # V.7.3-60
library(gridExtra) # V.2.3
library(rpart) # V.4.1.21
library(randomForest) # V.4.7-1.1
library(terra) # V.1.7-55
library(ggplot2) # V.3.4.4
library(corrplot) # V.0.92
library(caret) # V.6.0-94
library(gtools) # V.3.9.4
library(reshape2) # V.1.4.4
library(gam) # V.1.22-2
library(ecospat) # V.4.0.0
library(sf) # V.1.0-14
library(viridis) # V.0.6.4

# load data
pa_data <- read.table ("sdm/speciesData/Soil_Microorganisms/Fungi/Craterellus_tubaeformis.txt", row.names = NULL)
head(pa_data)
pa_data <- pa_data[, -1]
View(pa_data)
names(pa_data)

coord <- pa_data[,1:2]
head(coord)

# environmental data
raster_files <- mixedsort(list.files(path = "sdm/Covariate_currentTime",
                                     full.names = TRUE))
env_variables <- rast(raster_files)
# crs do not match

# fixing the coordinate system issue and playing around with rasters ----
crs(env_variables)

# Get the CRS of each raster
raster_crs_list <- sapply(raster_files, function(f) crs(rast(f)))

# Check which ones are different
unique_crs <- unique(raster_crs_list)
length(unique_crs)

# Take the first CRS as reference
ref_crs <- raster_crs_list[[1]]

# Compare each one to reference
which(sapply(raster_crs_list, function(crs) crs != ref_crs))

# Define the target CRS (EPSG:2056)
target_crs <- "EPSG:2056"

# Reproject each raster to the target CRS
reprojected_rasters <- lapply(raster_files, function(file) {
  r <- rast(file)
  if (crs(r) != target_crs) {
    r <- project(r, target_crs)
  }
  return(r)
})

# Check the CRS of the first raster after reprojection
crs(reprojected_rasters[[1]])

env_variables <- rast(reprojected_rasters)
 
# is this ok????
# hope so. Lol

cat(crs(env_variables))

# change the names of the rasters
names(env_variables) = sub("_.*", "", names(env_variables))

# are the rasters aligned? yes fab -----
# Check if all rasters have the same resolution and extent

# Check if all rasters have the same extent and resolution
same_extent <- all(sapply(env_variables, function(r) all(ext(r) == ext(env_variables[[1]]))))
same_resolution <- all(sapply(env_variables, function(r) all(res(r) == res(env_variables[[1]]))))

if (same_extent && same_resolution) {
  print("Rasters are aligned.")
} else {
  print("Rasters need alignment.")
}


names(env_variables)

# plot climatic maps -----
plot(env_variables, col =rev(viridis::inferno(100, direction = -1)))
# too many variables to show them properly

# extract environmental data values???
env_values <- na.omit(data.frame(extract(env_variables, coord)))
head(env_values)
# idk if this is ok, also try without reprojecting the coordinates in the first step

dim(env_values)
str(env_values)

# from first practical:
# now we have predictor data points corresponding to
# species data points
# (how ??)
# select which predictor to include in the model
# usually try: relevant for species, not correlates

# prelim analysis of correlation between predictors
# general rule: min 10 occurrences for 1 exp variable
# cor threshold: usually 65-75%

# create correlation matrix -----
df.cor.envdata <- data.frame(round(cor(env_values), 3))
head(df.cor.envdata)
head(env_values)
# remove id column - hope this doesnt cause issues idk
env_values <- env_values %>% dplyr::select(-ID)
View(env_values)
corrplot.mixed(cor(env_values), number.cex= 0.5, tl.cex = 0.5, tl.pos = "lt")
# this looks wrong lol

# select least correlated variables based on threshold of 0.7

colnames(df.cor.envdata) [-c(findCorrelation(as.matrix(df.cor.envdata), cutoff =0.70))]

cor.clust <- hclust(as.dist(1-abs(df.cor.envdata)))
# plot these distances as a tree
par(mar = c(1,4,3,1))
plot(cor.clust, main = "Cluster of the correlations among variables",
     ylab = "Height as 1-abs(cor)", xlab = "")
# threshold value
abline(h = 0.3, lty = 2, col = "red", lwd = 2)

# and now is the point where I can select the variables, in theory

nrow(pa_data)

str(env_values)

cor(env_values[,c( "bio13", "bio15", "f", "r", "forestaggr.1", "slope")])

cor(env_values[,c("bio1", "bio15", "forestaggr", "slope", "r", "f",  "twi25ss", "srad", "hillshade")])
# h, n remove
# maybe change forestaggr - see which is the best predictor (but 500 focal too correlated)
# keep bio13, drop bio1 - mabye max T?
# srad ok - but is it rly reasonable?

# try to see if there is relationship?? idk

# modelling prep----
species_data <- na.omit(data.frame(coord, presence=pa_data$resp,
                             extract(env_variables, coord)))
head(species_data)
View(species_data)

# try to start modelling ----
# glm -----
# Simple GLM with linear predictors
glm1 <- glm(presence ~ bio15 + gdd3 + f +r+ forestaggr + slope  + srad + hillshade,
            data=species_data, family="binomial")
summary(glm1)
# r is the only significant predictor

# Simple GLM with quadratic predictors
glm2 <- glm(presence ~ I(bio15^2)+I(gdd3^2)+I(f^2)  +I(r^2)  +I(forestaggr^2)+I(slope^2)
             +I(srad^2) + +I(hillshade^2),
            data=species_data, family="binomial")
summary(glm2)
# same
anova(glm1, glm2, test = "Chisq")

# A bit more complex GLM : both linear and quadratic predictors
glm3  <- glm(presence ~ bio15 + I(bio15^2) + gdd3 + I(gdd3^2) + f + I(f^2) + 
                              r + I(r^2) + forestaggr + I(forestaggr^2) + 
                              slope + I(slope^2)   +
             srad +I(srad^2)+ hillshade + I(hillshade^2),
             data = species_data,family=binomial)
summary(glm3)
# idk i think this would be too much, too many predictors
# lol this is so much worse...

# Centering the continuous variables
species_data$gdd3_centered <- scale(species_data$gdd3, center = TRUE, scale = FALSE)
env_variables$gdd3_centered <- scale(env_variables$gdd3, center = TRUE, scale = FALSE)
# Create the squared term for the centered gdd3 variable
species_data$gdd3_squared <- species_data$gdd3_centered^2
env_variables$gdd3_squared <- env_variables$gdd3_centered^2
# let this go bc fucks up the projections


# Re-run the model with the centered gdd3 and its squared term
glm_final <- glm(presence ~ bio15 + gdd3 + I(gdd3^2) + f + 
                   r + forestaggr + slope + srad + hillshade,
                 data = species_data, family = binomial)
summary(glm_final)

anova(glm1, glm3, test = "Chisq")

vif(glm_final)
# ok moving on ----
# predictions of relationship between variables and p/a in glm -----
var.names <- c("bio15", "gdd3", "f", "r", "forestaggr", "slope",  "srad", "hillshade")
# fix other values
medians<- apply(species_data[, var.names],2,median)

medians_table <- data.frame(sapply(medians, function(x)rep(x, 100)))
summary(medians_table)

# create gradient for bio13
foc.var <- spData[,"bio13"]
var.new <- seq(min(foc.var), max(foc.var), length = 100)
# using 100 points
new.data <- medians_table
new.data[,"bio13"] <- var.new # put the gradient in the table of medians
summary(new.data)
head(new.data)

# ask the 3 glm-s to predict on the new ("fake") dataset
pred.glm1 <-predict(glm1,newdata=new.data,type="response")
pred.glm2<-predict(glm2,newdata=new.data,type="response")
pred.glm_final<-predict(glm_final,newdata=new.data,type="response")

pred.glob <- c(pred.glm1, pred.glm2)

tmp1<-cbind(Occ.prob=pred.glob,Env.val=var.new)
tmp2 <-data.frame(cbind(Algorithm = rep(c("GLM1", 'GLM2'), each = 100), 
                        Var.name = "bio13"))
resp1.glm<-ggplot(cbind(tmp1,tmp2), aes(x=Env.val,y=Occ.prob,color=Algorithm,
                                        linetype=Algorithm)) +
  geom_line(lwd=1) +
  scale_color_manual(values=c("#007991","#439A86","#BCD8C1"))+
  facet_wrap(~Var.name,ncol=2,scale="free_x")+
  theme_classic()+
  labs(x="Value",y="Occurence probability")
resp1.glm

# plotting rspc for all  vars
out1.glm = NULL
for(i in 1:(length(var.names))) {
  foc.var <- species_data[,var.names[i]]
  new.data <- medians_table
  var.new <- seq(min(foc.var), max(foc.var), length = 100)
  new.data[,i] <- var.new
  pred.glm1 <- predict (glm1, newdata = new.data, type = "response")
  pred.glm2 <- predict (glm2, newdata = new.data, type = "response")
 # pred.glm3 <- predict(glm3, newdata = new.data, type = "response")
  pred.glob <- c(pred.glm1, pred.glm2)
  tmp1 <- cbind(Occ.prob= pred.glob, Env.val=var.new)
  tmp2 <- data.frame(cbind(Algorithm = rep(c("GLM1", "GLM2"), each = 100),
                           Var.name = var.names[i]))
  out1.glm <- rbind(out1.glm, cbind(tmp1, tmp2))
}
# ones that look like without effect: bio15 (kinda), forestaggr, slope, twi25ss

resp1.glm<-ggplot(out1.glm,aes(x=Env.val,y=Occ.prob,color=Algorithm,
                               linetype=Algorithm)) +
  geom_line(lwd=1) +
  scale_color_manual(values=c("#007991","#439A86","#BCD8C1"))+
  facet_wrap(~Var.name,ncol=2,scale="free_x")+
  theme_bw()+
  labs(x="Value",y="Occurence probability")+
  theme_bw()
  print(resp1.glm)

ggsave("sdm/report/plots/glm_curves_initial_v2.png", dpi = 600)

# response curves for the final glm only ----

out_final <- NULL

for(i in 1:length(var.names)) {
  foc.var <- species_data[, var.names[i]]
  
  # Create new dataset for predictions
  new.data <- medians_table
  var.new <- seq(min(foc.var), max(foc.var), length.out = 100)
  new.data[, i] <- var.new
  
  # Predict using glm_final
  pred <- predict(glm_final, newdata = new.data, type = "response")
  
  # Store output
  tmp <- data.frame(
    Occ.prob = pred,
    Env.val = var.new,
    Var.name = var.names[i],
    Algorithm = "GLM_final"
  )
  out_final <- rbind(out_final, tmp)
}

# Plot
resp_glm_final <- ggplot(out_final, aes(x = Env.val, y = Occ.prob)) +
  geom_line(color = "#007991", lwd = 1) +
  facet_wrap(~Var.name, ncol = 2, scales = "free_x") +
  theme_bw() +
  labs(x = "\nValue", y = "Occurrence probability\n")

print(resp_glm_final)

ggsave("sdm/report/plots/responses_glm.png", dpi = 600)

# slope looks unnecessary
# can keep these minus slope for now, but also play around with other variables later
# also maybe bio13 not as important??
# maybe bio1 instead of bio13

# final model: only include quadratic term (as well as linear) for gdd. 

# stepwise parameter selection -----

glmStart <- glm(presence~1, data = species_data, family = binomial)
glm.formula<-bm_MakeFormula("presence",
                            species_data[,c("bio15", "gdd3", "f", "r", "forestaggr", "slope", "twi25ss", "srad", "hillshade")],
                            "quadratic",interaction.level=1)
# stepwise selection with AIC
glmModAIC <- stepAIC(glmStart, glm.formula,
                     direction = "both",
                     trace = TRUE,
                     k = 2)

# only keeps gdd and r

# Stepwise selection with BIC
glmModBIC<-stepAIC(glmStart,
                   glm.formula,
                   direction="both",
                   trace=FALSE,
                   k=log(nrow(species_data)))
glmModBIC
?stepAIC

out2.glm=NULL
for(i in 1:length(var.names)){
  foc.var<-species_data[,var.names[i]]
  new.data<-medians_table
  var.new<-seq(min(foc.var),max(foc.var),length=100)
  new.data[,i]<-var.new
  pred.glm1<-predict(glm1,newdata=new.data,type="response")
  pred.glm2<-predict(glm2,newdata=new.data,type="response")
#  pred.glm3<-predict(glm3,newdata=new.data,type="response")
  pred.glmAIC<-predict(glmModAIC,newdata=new.data,type="response")
  pred.glmBIC<-predict(glmModBIC,newdata=new.data,type="response")
  pred.glob <-c(pred.glm1,pred.glm2,pred.glmAIC,pred.glmBIC)
  tmp1<-cbind(Occ.prob=pred.glob,Env.val=var.new)
  tmp2<-data.frame(cbind(Algorithm=rep(c("GLM1","GLM2","glmAIC","glmBIC"),
                                       each=100),Var.name=var.names[i]))
  out2.glm<-rbind(out2.glm,cbind(tmp1,tmp2))
}
resp2.glm<-ggplot(out2.glm,aes(x=Env.val,y=Occ.prob,color=Algorithm,
                               linetype=Algorithm)) +
  geom_line(linewidth=1)+
  scale_color_manual(values=c("#222E50", "#007991", "#439A86", "#BCD8C1",
                              "#E9D985"))+
  facet_wrap(~Var.name,ncol=2,scale="free_x")+
  theme_bw()+
  labs(x="Value",y="Occurence probability")+
  theme_bw()
print(resp2.glm)

# for bio13 the curve changes !! it becomes a hump !!

# recursive partitioning ----
RP <- rpart(presence ~bio13 + bio15 + f + forestaggr.1 + r,
            data = species_data,
            control=rpart.control(xval=1000), method = "class")

par(mfrow=c(1,1), mar=c(0,0,0,0))
plot(RP, uniform=F, margin=0.1, branch=1)
text(RP)

# random forest-----
rf_final = randomForest(x = species_data[,c("bio15", "f", "r", "forestaggr", "gdd3", "hillshade",  "slope", "srad")],
                  y = as.factor(species_data$presence),
                  ntree = 1000,
                  importance =TRUE)
RF.pred = predict(rf_final, type = "prob")[,2]

importance(rf_final)
?importance


out.rf <- NULL
for (i in 1:length(var.names)){
  foc.var <- species_data[,var.names[i]]
  new.data <- medians_table
  var.new <- seq(min(foc.var), max(foc.var), length = 100)
  new.data[,i] <- var.new
  pred.rf <- predict(RF, newdata = new.data, type = "prob")[,2]
  tmp1<-cbind(Occ.prob=pred.rf,Env.val=var.new)
  tmp2<-data.frame(cbind(Algorithm=rep(c("RF"),each=100),Var.name=var.names[i]))
  out.rf<-rbind(out.rf,cbind(tmp1,tmp2))
}
resp.rf<-ggplot(out.rf, aes(x=Env.val,y=Occ.prob)) +
  geom_line(size=1, color = "#007991")+
  scale_color_manual(values=c("#007991"))+
  facet_wrap(~Var.name,ncol=2,scale="free_x")+
  theme_bw()+
  labs(x="\nValue",y="Occurence probability\n")+
  theme_bw()
  print(resp.rf)
  
  ggsave("sdm/report/plots/rf_final.png", dpi = 600)


# gam with stepwise selection ----
gamStart<-gam::gam(presence~1,data=species_data,family=binomial)


gam_df2 <- gam::step.Gam(gamStart,
                         list(
                           ~1 + gam::s(bio13, 2),
                           ~1 + gam::s(bio15, 2),
                           ~1 + gam::s(f, 2),
                           ~1 + gam::s(forestaggr.1, 2),
                           ~1 + gam::s(r, 2)
                         ),
                         trace = TRUE, direction = "both")

gam_df2$anova

# wit df 4
gam_df4 <- gam::step.Gam(gamStart,
                         list(
                           ~1 + gam::s(bio13, 4),
                           ~1 + gam::s(bio15, 4),
                           ~1 + gam::s(f, 4),
                           ~1 + gam::s(forestaggr.1, 4),
                           ~1 + gam::s(r, 4)
                         ),
                         trace = TRUE, direction = "both")
gam_df4$anova

gam::step.Gam


##with df = 8.8 (df can be a decimal number)
gam_df8 <- gam::step.Gam(gamStart,
                         list(
                           ~1 + gam::s(bio13, 8),
                           ~1 + gam::s(bio15, 8),
                           ~1 + gam::s(f, 8),
                           ~1 + gam::s(forestaggr.1, 8),
                           ~1 + gam::s(r, 8)
                         ),
                         trace = TRUE, direction = "both")
gam_df8$anova

# if you increase df, allow more flexibility BUT risk overfitting.
# chatgpt suggests using mcv:gam instead of this.


# response curves
out.gam <- NULL
#prepare median table
# var.names <- c("bio13","bio15","f", "forestaggr.1", "r")
medians <- apply(species_data[,var.names],2,median)
medians_table <- data.frame(sapply(medians,function(x)rep(x,100)))
#get predictions for each variable gradient for each model
for(i in 1:length(var.names)){
  foc.var <- species_data[,var.names[i]]
  new.data <- medians_table
  var.new <- seq(min(foc.var), max(foc.var), length=100)
  new.data[,i] <- var.new
  pred.gam_df2 <- predict(gam_df2, newdata=new.data, type="response")
  pred.gam_df4 <- predict(gam_df4, newdata=new.data, type="response")
  pred.gam_df8 <- predict(gam_df8, newdata=new.data, type="response")
  pred.glob <- c(pred.gam_df2,pred.gam_df4,pred.gam_df8.8)
  tmp1 <- cbind(Occ.prob=pred.glob,Env.val=var.new)
  tmp2 <- data.frame(cbind(Algorithm=rep(c("GAM_df2","GAM_df4","GAM_df8.8"),each=100),
                           Var.name=var.names[i]))
  out.gam <- rbind(out.gam,cbind(tmp1,tmp2))
}

# plot results
resp.gam <- ggplot(out.gam, aes(x=Env.val, y=Occ.prob, color=Algorithm,
                                linetype=Algorithm)) +
  geom_line(linewidth=1) +
  scale_color_manual(values=c("#007991", "#439A86", "#BCD8C1"))+
  facet_wrap(~Var.name,ncol=2,scale="free_x")+
  theme_bw()+
  labs(x="Value",y="Occurence probability")+
  theme_bw()
print(resp.gam)

# rising value of df -> model can take more complex shape
# sometimes simple models better, somtimes increasing complexity helps
# but overfitting, beware

# final gam model ----
library(mgcv)
gam_final <- gam(presence ~ s(bio15) + s(gdd3) + s(f) + s(r) + 
                   s(forestaggr) + s(slope) + 
                   s(srad) + s(hillshade),
                 data = species_data, family = binomial)
summary(gam_final)

# Compute medians for all variables
medians <- apply(species_data[, var.names], 2, median, na.rm = TRUE)
medians_table <- as.data.frame(sapply(medians, function(x) rep(x, 100)))

# Generate predictions across the gradient of each variable
out.gam <- NULL

for (i in 1:length(var.names)) {
  foc.var <- species_data[[var.names[i]]]
  var.new <- seq(min(foc.var, na.rm = TRUE), max(foc.var, na.rm = TRUE), length = 100)
  
  new.data <- medians_table
  new.data[[var.names[i]]] <- var.new
  
  pred.gam <- predict(gam_final, newdata = new.data, type = "response")
  
  tmp1 <- data.frame(Occ.prob = pred.gam, Env.val = var.new)
  tmp2 <- data.frame(Algorithm = "GAM_final", Var.name = var.names[i])
  out.gam <- rbind(out.gam, cbind(tmp1, tmp2))
}

resp.gam <- ggplot(out.gam, aes(x = Env.val, y = Occ.prob)) +
  geom_line(color = "#007991", linewidth = 1) +
  facet_wrap(~Var.name, ncol = 2, scales = "free_x") +
  theme_bw() +
  labs(x = "\nValue", y = "Occurrence probability\n")
print(resp.gam)

ggsave("sdm/report/plots/gam_final_plots.png", dpi = 600)

# model evaluation ----
# libraries ----
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

# run models for which to evaluate performance

# glm
#GLM <- glm(formula = presence ~ bio13 +  I(bio13^2) + (bio15) + r,
#          family = binomial, data = species_data)
GLM.pred <- predict (GLM, newdata = species_data[,c("bio13", "bio15", "r")],
                     type = "response")

# gam
GAM <- gam(presence ~ s(bio13) + s(bio15) + s(r),
           data = species_data, family = "binomial")
GAM.pred <- predict(GAM, newdata = species_data[,c("bio13", "bio15", "r")],
                    type = "response")
# RF
rf = randomForest (x = species_data[,c("bio13", "bio15", "r")],
                   y = as.factor(species_data$presence),
                   ntree = 1000)
rf.pred = predict(rf, type = "prob") [,2] # extract second column - probabilities


ObsPA <- species_data$presence
plotID <- 1:nrow(species_data)
EvalData <- data.frame(cbind(plotID, ObsPA, GLM.pred, GAM.pred, RF.pred))
colnames(EvalData) <- c ("plotID", "ObsPA", "GLM", "GAM", "RF")
head(EvalData)

# final models and predictions for evaluating perf----
glm_final
GLM.pred =predict(glm_final, newdata = species_data [,c("bio15",
 "gdd3", "srad","r", "f", "hillshade", "slope", "forestaggr")], 
 type = "response")

gam_final
GAM.pred =predict(gam_final, newdata = species_data [,c("bio15",
            "gdd3", "srad","r", "f", "hillshade", "slope", "forestaggr")], 
                  type = "response")

rf_final
RF.pred = predict(rf_final, type="prob")[,2]

# final dataset with all predictions
ObsPA <- species_data$presence
plotID <- 1:nrow(species_data)
EvalData <- data.frame(cbind(plotID, ObsPA, GLM.pred, GAM.pred, RF.pred))
colnames(EvalData) <- c ("plotID", "ObsPA", "GLM", "GAM", "RF")
head(EvalData)


# measuring model accuracy ----
# calibration ----
# extent to which model correctly predicts presence

par(oma = c(2, 2, 0, 0), mar = c(2, 2, 2, 1), mfrow = c(2, 2),
    cex = 0.7, cex.lab = 1.4, mgp = c(2, 0.5, 0))
for (mod in 1:3) {
  calibration.plot(EvalData, which.model = mod, color = mod + 1,
                   xlab = "", ylab = "", N.bins=10)
  mtext("\nPredicted Probability of presence", side = 1, line = 1,
        cex = 1.4, outer = TRUE)
  mtext("Proportion of observed presence\n", side = 2, line =-4,
        cex = 1.4, outer = TRUE)
}

data(SIM3DATA)
calibration.plot(SIM3DATA, which.model = 1, color = 5,
                 xlab = "", ylab = "", N.bins = 10,
                 main = "Observed vs Predicted (Bad calibration")
?calibration.plot

# gam seems the best

# discrimination -----
# compare continuous predictions to the binary observations
# binarize predictions
# choosing threshold

table(EvalData$GAM>0.5, EvalData$ObsPA, dnn= c("Prediction", "Observation"))
accu <- presence.absence.accuracy(EvalData, which.model = 2, threshold = 11,
                                  st.dev = FALSE)
# good at predicting absences but not great at predicting presences.
# also not many presences in general - could be an issue.
accu[,-c(1, 2)] <- signif(accu[,-c(1, 2)], digits = 2)
accu

pred.prev <- predicted.prevalence(EvalData, threshold = 11)
pred.prev[, 2:5] <- round(pred.prev[, 2:5], digits = 2)
pred.prev

ecospat.meva.table(Pred = EvalData$GAM, Sp.occ = EvalData$ObsPA, th = 0.5)

kappa.max <- ecospat.max.kappa(Pred = EvalData$GAM, Sp.occ = EvalData$ObsPA)
head(kappa.max$table)

kappa.max$max.Kappa
# suggests moderate agreement
# above, this, predict presence, below this: predict absence

tss.max <- ecospat.max.tss(Pred = EvalData$GAM, Sp.occ = EvalData$ObsPA)
tss.max$max.TSS
# a better threshold bc not affected by prevalence.
head(tss.max$table)

tmpKappa<-NULL
tmpTSS<-NULL
for(i in 3:5){
  inter<-cbind.data.frame(ecospat.max.kappa(Pred=EvalData[,i],
                                            Sp.occ=EvalData$ObsPA)$table,
                          Model= colnames(EvalData)[i],
                          Evaluation.Metric="Kappa")
  tmpKappa<-rbind.data.frame(tmpKappa,inter)
  inter<-cbind.data.frame(ecospat.max.tss(Pred=EvalData[,i],
                                          Sp.occ=EvalData$ObsPA)$table,
                          Model= colnames(EvalData)[i],
                          Evaluation.Metric="TSS")
  tmpTSS<-rbind.data.frame(tmpTSS, inter)
}
colnames(tmpKappa)= colnames(tmpTSS)=c("Threshold","Value","Model",
                                       "Evaluation.Metric")
tmp<-rbind.data.frame(tmpKappa,tmpTSS)
ggplot(tmp,aes(x=Threshold,y=Value,color=Model))+
  geom_line(linewidth=1.1)+
  scale_color_discrete(type=c("#007991","#439A86","#BCD8C1"))+
  theme_bw()+
  facet_wrap(~Evaluation.Metric)

# gam performs best

# error threshold plots
data <- EvalData[1:5]
N.models <- ncol(data)- 2
par(oma=c(0,0,0,0), mar=c(4,1,1.5,1), mfrow=c(2,2), cex=0.7, cex.lab=1.4, mgp=c(2, 0.5,0))
for (mod in 1:N.models){
  error.threshold.plot(data,
                       which.model = mod,
                       color = TRUE,
                       add.legend = TRUE,
                       legend.cex = 0.9)
}


par(mfrow  = c(1,1))
auc.roc.plot(data, color=c("#007991", "#439A86", "#BCD8C1"), legend.cex=0.7, main="")

# should be towards the top left corner

# calculate optimal threshold using different methods
optimal.thresholds(EvalData, opt.methods = 1:12, req.sens=0.9, req.spec = 0.9, FPC = 2, FNC = 1)
# can stick to tss



# DATA PARTITIONING STRATEGIES  / validation ----
# unsure if correct bc some things changed compared to tp script

# Ensure the data is correct
fungus <- data.frame(bio15 = species_data$bio15,
                     gdd3 = species_data$gdd3,
                     srad = species_data$srad,
                     f = species_data$f,
                     r = species_data$r,
                     forestaggr = species_data$forestaggr,
                     slope = species_data$slope,
                     hillshade = species_data$hillshade,
                     presence = as.factor(species_data$presence))

# k-fold cross validation
# 1 - for glm
set.seed(123)
cv.error.10 <- rep(0, 10)

for (i in 1:10) {
  # Fit the original GLM model with the quadratic term for gdd3
  glm.fit <- glm(presence ~ bio15 + gdd3 + srad + f + r + forestaggr + slope + hillshade + I(gdd3^2),
                 family = "binomial", 
                 data = fungus)
  
  # Cross-validation on the same model
  cv.error.10[i] <- cv.glm(fungus, glm.fit, K = 10)$delta[1]
}

cv.error.10
mean(cv.error.10)
# 0.1715536

#  rf model
set.seed(123)
fungus$presence<-make.names(fungus$presence)
# idk what this does lol
# define training control
# cv of 10
train_control <- trainControl(method = "cv", number = 10, savePredictions = T,
                              summaryFunction = twoClassSummary, classProbs = T)
# train model
rf_model <- train(presence ~ ., data = fungus, trControl = train_control, method = "rf")

# get cv results
(cv_results<- rf_model$results)
# highest accuracy: 
# plot roc
selectedIndices <- model$pred$mtry == 3
ROC <- roc(as.numeric(model$pred$obs[selectedIndices]),
           as.numeric(model$pred$X0[selectedIndices]), auc = TRUE)
confMat <- caret::confusionMatrix(data = model$pred$obs[selectedIndices],
                                  reference = model$pred$pred[selectedIndices],
                                  mode= "everything")
confMat$overall
par(mfrow=c(1,1))
plot.roc(smooth(ROC))
ROC$auc
# 0.7351

# mean AUC(ROC) pf 0-77+-0.11, sensitivity of 0.9 and specificity of 0.37 
# not good at predicting absences
# accuracy of 0.77, (0.708  0.826), kappa 0.40

# on gam - ??

library(caret)
library(pROC)

set.seed(123)

# Define train control
train_control <- trainControl(method = "cv", number = 10,
                              summaryFunction = twoClassSummary,
                              classProbs = TRUE, savePredictions = TRUE)

# Train the GAM model
model_gam <- train(presence ~ bio15+
                   gdd3+
                   srad +
                   f +
                   r+
                   forestaggr +
                   slope +
                   hillshade,
                   data = fungus,
                   method = "gam",
                   trControl = train_control,
                   metric = "ROC",
                   family = "binomial")

model_gam
library(pROC)

# Use all predictions across folds
roc_obj <- roc(model_gam$pred$obs, model_gam$pred$X1)

# Plot the ROC
plot.roc(roc_obj)

# If you want the AUC:
auc(roc_obj)
# 0.7025

# If you want the max TSS:
TSS <- roc_obj$sensitivities + roc_obj$specificities - 1
max(TSS) # 0.3347

cv_results <- model_gam$results
cv_results  # Mean AUC, Sensitivity, Specificity for different mtry values


# other data partitioning methods -----
fungus$presence <- make.names(fungus$presence)

#-- Usual bootstrap
train_control <- trainControl(method="boot", number=5, savePredictions=T,
                              summaryFunction=twoClassSummary, classProbs=T)
model.boot <- train(presence ~ ., data=fungus, trControl=train_control,
                    method="rf")
ROC.boot <- roc(as.numeric(model.boot$pred$obs[selectedIndices]),
                as.numeric(model.boot$pred$X0[selectedIndices]), auc=TRUE)

#-- the 0.632 bootstrap estimator (Efron 1983)
train_control <- trainControl(method="boot632", number=5, savePredictions=T,
                              summaryFunction=twoClassSummary,classProbs=T)
model.boot632 <- train(presence ~ ., data=fungus,
                       trControl=train_control, method="rf")
ROC.boot632 <- roc(as.numeric(model.boot632$pred$obs[selectedIndices]),
                   as.numeric(model.boot632$pred$X0[selectedIndices]), auc=TRUE)
#-- optimism bootstrap estimator (Efron and Tibshirani, 1994)
train_control <- trainControl(method="optimism_boot", number=5, savePredictions=T,
                              summaryFunction=twoClassSummary,classProbs=T)
model.boot.optim <- train(presence ~ ., data=fungus,
                          trControl=train_control, method="rf")
ROC.boot.optim <- roc(as.numeric(model.boot.optim$pred$obs[selectedIndices]),
                      as.numeric(model.boot.optim$pred$X0[selectedIndices]),auc=TRUE)
#-- 5-fold cross validation
train_control <- trainControl(method="cv", number=5, savePredictions=T,
                              summaryFunction=twoClassSummary, classProbs=T)
model.cv <- train(presence ~ ., data=fungus, trControl=train_control,
                  method="rf")
ROC.cv <- roc(as.numeric(model.cv$pred$obs[selectedIndices]),
              as.numeric(model.cv$pred$X0[selectedIndices]), auc=TRUE)
#-- repeated split-sample CV where 75% of the data is used for calibrating and 25% for the evaluation
# You can change the ratio with the parameter p
train_control <- trainControl(method="LGOCV", number=5, savePredictions=T,
                              summaryFunction=twoClassSummary, classProbs=T, p= 0.75)
model.LGOCV <- train(presence ~ ., data=fungus, trControl=train_control,
                     method="rf")
ROC.LGOCV <- roc(as.numeric(model.LGOCV$pred$obs[selectedIndices]),
                 as.numeric(model.LGOCV$pred$X0[selectedIndices]), auc=TRUE)
#-- Plot the results
plot.roc(smooth(ROC.boot),col=1)
plot.roc(smooth(ROC.boot632),add=T,col=2)
plot.roc(smooth(ROC.boot.optim),add=T,col=3)
plot.roc(smooth(ROC.cv),add=T,col=4)
plot.roc(smooth(ROC.LGOCV),add=T,col=5)
legend("bottomright", c("boot","boot632","boot.optim","cv","LOGCV"),
       col=1:5, lty=1, inset=0.01)

# ???

# comparison of different algorithms----
nCV <- 20 # number of cross validation
nRow <- nrow(species_data)
Test_results <- as.data.frame(matrix(0,ncol=nCV,nrow=3,
                                     dimnames=list(c("GLM","GAM","RF"),
                                                   NULL)))
# array to store predicted habitat suitability
Pred_results <- array(0, c(nRow, 3, nCV),
                      dimnames=list(seq(1:nRow),
                                    c("GLM","GAM","RF"),
                                    seq(1:nCV)))
for(i in 1:nCV){
  # separate the original data in one subset for calibration and another for evaluation.
  a <- bm_SampleBinaryVector(obs = species_data$presence, ratio=0.7)
  calib <- species_data[a$calibration,]
  # Center the gdd3 variable
  calib$gdd3_centered <- scale(calib$gdd3, center = TRUE, scale = FALSE)
  # Create the squared term
  calib$gdd3_squared <- calib$gdd3_centered^2
  
  
  eval <- species_data[a$validation,]
  ### GLM ###
  # taken out stepwise
  # Use your final GLM model here instead of the stepwise model
  glm_model <- glm(presence ~ bio15 + gdd3 + I(gdd3^2) + f + r + forestaggr +
                     slope + srad + hillshade, 
                   data = calib, family = binomial,
                   control = glm.control(maxit = 100))
  
  # prediction on the evaluation data and evaluation using the AUC approach
  Pred_test <- predict(glm_model, eval, type = "response")
  Test_results["GLM", i] <- as.numeric(auc(roc(eval$presence, Pred_test)))
  
  # prediction on the total dataset
  Pred_results[,"GLM",i] <- predict(glm_model, species_data, type="response")
  ### GAM ###
  gam_mgcv <- gam(presence ~ s(bio15) + s(gdd3) + s(r)+ s(f)+ s(forestaggr)
                  + s(slope)+ s (srad) + s(hillshade),
                  data=calib, family="binomial")
  # prediction on the evaluation data and evaluation using the AUC approach
  Pred_test <- predict(gam_mgcv, eval, type="response")
  Test_results["GAM",i] <- as.numeric(auc(roc(eval$presence,Pred_test)))
  # prediction on the total dataset
  Pred_results[,"GAM",i] <- predict(gam_mgcv, species_data, type="response")
  ### RF ###
  RF_mod = randomForest(x = calib[,c("bio15","gdd3", "f", "r", "forestaggr",
                                     "slope" ,"srad" ,"hillshade" )],
                        y = as.factor(calib$presence),
                        ntree = 1000, importance = TRUE)
  # prediction on the evaluation data and evaluation using the AUC approach
  Pred_test <- predict(RF_mod, eval, type="response")
  Test_results["RF",i] <- as.numeric(auc(roc(eval$presence,
                                             as.numeric(as.character(Pred_test)))))
  # prediction on the total dataset
  Pred_results[,"RF",i] = predict(RF_mod, species_data, type="prob")[,2]
  print(i)
}


# variation in auc between models & CV runs
AUC <- unlist(Test_results)
AUC <- as.data.frame(AUC)
Test_results_ggplot <- cbind(AUC, model=rep(rownames(Test_results), times=nCV))
p <- ggplot(Test_results_ggplot, aes(model, AUC))
p + geom_boxplot() + theme_bw()
# model with the highest AUC - best predictive performance

ddply(Test_results_ggplot,.(model), summarize,
      Mean_AUC = mean(AUC), SD_AUC = sd(AUC),
      Median_AUC = median(AUC))

# a simple glm is best??

# more libraries ----
library(ecospat)
library(mgcv)
library(ade4)
library(terra)
library(fields)
library(dismo)
library(gtools)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(sf)
library(biomod2)
library(randomForest)

# detach gam for mgcv gam to work
if( "package:gam" %in% search()){
  detach("package:gam", unload = TRUE)
}

# missed 1 step - remove coordinates without climatic info. 
# if things not working then maybe do this step

gam_simplified <- mgcv::gam(presence ~ s(bio13,k=4) + s(bio15,k=4) + s(r,k=4),
                  data=species_data, family="binomial")
glm1 <- glm(formula = presence ~ bio13 + I(bio13^2) + bio15 + r,
            family = binomial, data = species_data)
rf1 = randomForest(x = species_data[,c("bio13", "bio15", "r")],
                   y = as.factor(species_data$presence),
                   ntree = 1000)

# but i can just use my models

# Predicted occurence probabilities

# make names unique
names(env_variables) <- make.unique(names(env_variables))
relevant_layers <- env_variables[[c("bio15", "gdd3", "gdd3_centered", "f", "r", "forestaggr", "slope", "srad", "hillshade")]]

fungus.curr.gam <- predict(env_variables, gam_final, type="response")
tmp.fungus.gam.df <- cbind.data.frame(Model = "GAM",
                                    as.data.frame(fungus.curr.gam,
                                                  xy=T, na.rm=T))

fungus.curr.glm <- predict(env_variables, glm_final, type="response")
tmp.fungus.glm.df <- cbind.data.frame(Model = "GLM",
                                    as.data.frame(fungus.curr.glm,
                                                  xy=T, na.rm=T))

fungus.curr.rf <- predict(env_variables, rf_final, type="prob", na.rm=T)[[2]]
tmp.fungus.rf.df <- cbind.data.frame(Model = "RF",
                                   as.data.frame(fungus.curr.rf,
                                                 xy=T, na.rm=T))
colnames(tmp.fungus.rf.df) = c("Model","x","y","lyr1")


# dataset for plotting

PlotData <- rbind.data.frame(tmp.fungus.gam.df,tmp.fungus.glm.df,
                             tmp.fungus.rf.df)

# plot the new spatial patterns
pal1 <- rev(colorRampPalette(brewer.pal(11,"RdYlGn")[11:1])(64))
(predicted_spatial_patterns <- ggplot(PlotData, aes(x=x, y=y, fill=lyr1)) +
  geom_raster(alpha=0.8) +
  facet_wrap(~Model) +
  coord_equal()+
  scale_fill_gradientn("Occurence\nprobability",colours=pal1)+
  theme_void())

# predict SE for GAM (also available for GLM) for the df and plot the predictions

fungus.se <- predict(env_variables,gam1, type="response", se.fit=TRUE)
df.gam <- as.data.frame(fungus.se, xy=T, na.rm=T)
colnames(df.gam) <- c("x","y","GAM.pred","GAM.sd")
gam.pred <- ggplot(df.gam, aes(x=x, y=y, fill=GAM.pred)) +
  geom_raster(alpha=0.8) +
  coord_equal()+
  scale_fill_gradientn(colours=pal1)+
  theme_void()
pal2 <- colorRampPalette(brewer.pal(9,"YlOrRd"))(64)
gam.sd <- ggplot(df.gam, aes(x=x, y=y, fill=GAM.sd)) +
  geom_raster(alpha=0.8) +
  coord_equal()+
  scale_fill_gradientn(colours=pal2)+
  theme_void()
plot_grid(gam.pred, gam.sd, labels=c("A", "B"), ncol=1)


# ensemble modelling ----

# Combine the maps into a stack
fungus.curr <- c(fungus.curr.gam, fungus.curr.glm,fungus.curr.rf)

fungus.curr.mean <- mean(fungus.curr)
tmp.fungus.mean.df <- cbind.data.frame(Ens = "Mean",
                                     as.data.frame(fungus.curr.mean,
                                                   xy=T, na.rm=T))
colnames(tmp.fungus.mean.df) = c( "Ens","x", "y", "Prob")

# Standard deviation
fungus.curr.sd <- stdev(fungus.curr)
tmp.fungus.sd.df <- cbind.data.frame(Ens = "sd",
                                   as.data.frame(fungus.curr.sd,
                                                 xy=T, na.rm=T))
colnames(tmp.fungus.sd.df) = c( "Ens","x", "y", "Prob")
# Weighted Mean
#Here a weighted mean is perform based on the mean AUC value
#obtained during practical 2 based on 20 CV (can slightly change)
fungus.curr.wmean <- weighted.mean(fungus.curr,
                                 w= c(0.6997442 ,0.7613171 , 0.6175831 ))
tmp.fungus.wmean.df <- cbind.data.frame(Ens = "wMean",
                                      as.data.frame(fungus.curr.wmean,
                                                    xy=T, na.rm=T))
colnames(tmp.fungus.wmean.df) = c( "Ens","x", "y", "Prob")

# make dataset for ggplot
PlotData <- rbind.data.frame(tmp.fungus.mean.df,
                             tmp.fungus.sd.df,
                             tmp.fungus.wmean.df)
# Update map names (i.e., factor levels in 'Ens')
PlotData$Ens <- factor(PlotData$Ens,
                       levels = c("Mean", "sd", "wMean"),
                       labels = c("Simple Mean", "Model Disagreement (sd)", "Weighted Mean"))

# Plot the spatial patterns
pal1 <- rev(colorRampPalette(brewer.pal(11,"RdYlGn")[11:1])(64))
(ensemble_plot <- ggplot(PlotData, aes(x=x, y=y, fill=Prob)) +
  geom_raster(alpha=0.8) +
  facet_wrap(~Ens) +
  coord_equal()+
  scale_fill_gradientn("Occurence\nprobability",colours=pal1)+
  theme_void())
ggsave("sdm/report/plots/ensemble_maps.png", dpi= 600)

# predict in a new area - unsure what this is??-----
envData <- as.data.frame(env_variables,xy=T)
tmp1 <- dudi.pca(envData[,c("bio15", "bio13" , "r")],
                 nf=2, scannf=F)
tmp2 <- data.frame(cbind(envData, tmp1$li))
tmp2$ID <- ifelse(envData$x<2570000, 0, 1)
(map.plot <- ggplot(tmp2, aes(x=x, y=y,fill=factor(ID))) +
  geom_raster(alpha=1)+
  coord_equal()+
  scale_fill_manual(values=c("#3090C733","#9F000F4D"),
                    labels = c("Valleys","Mountains"),
                    name="")+
  theme_void()+
  guides(colour=guide_legend(override.aes=list(size=2))))

# wtf is this


# Environmental conditions in the two worlds
pca.plot <- ggplot(tmp2,aes(x=Axis1, y=Axis2,fill=factor(ID)))+ stat_bin2d(alpha=0.5)+
  theme_bw()+
  labs(colour="",x="PCA-Axis 1", y="PCA-Axis 2")+
  scale_fill_manual(values=c("#3090C733","#9F000F4D"),
                    labels = c("Valleys","Mountains"),
                    name="")
plot_grid(map.plot, pca.plot, labels=c("A", "B"), ncol=1)

Valleys <- envData[envData$x<2570000,]

mess.fungus <- ecospat.mess(envData[,c("x", "y",
                                         "bio13", "bio15", "r")],
                          Valleys[,c("x", "y",
                                     "bio13", "bio15", "r")])
mess.fungus <- as.data.frame(mess.fungus)
# Plot the MESS outputs

pal3 <- rev(colorRampPalette(brewer.pal(11,"RdYlGn"))(64))
pal4 <- colorRampPalette(brewer.pal(11,"RdYlGn"))(64)
mess <- ggplot(data=mess.fungus, aes(x=x, y=y, fill=MESS)) +
  geom_raster() +
  scale_fill_gradientn(colours=pal4)+
  coord_equal()+
  theme_void()
mess.w <- ggplot(data=mess.fungus, aes(x=x, y=y, fill=MESSw)) +
  geom_raster() +
  scale_fill_gradientn(colours=pal4)+
  coord_equal()+
  theme_void()
mess.neg <- ggplot(data=mess.fungus, aes(x=x, y=y, fill=MESSneg)) +
  geom_raster() +
  scale_fill_gradientn(colours=pal3)+
  coord_equal()+
  theme_void()
plot_grid(mess, mess.w, mess.neg, labels=c("A", "B", "C"), nrow=1)

# LITERALLY NO IDEA what this does

# projecting in time----

# used rcp 85, also try rcp 45

bio15.fu.rcp45 <- rast("sdm/covariates_future/bio15_20702099_RCP45.tif")
gdd3.fu.rcp45 <- rast("sdm/covariates_future/gdd3_20702099_RCP45.tif")
f.fu <-rast("sdm/Covariate_currentTime/f.tif")
r.fu <- rast("sdm/Covariate_currentTime/r.tif")
forestaggr.fu <- rast("sdm/Covariate_currentTime/forestaggr_100.tif")
slope.fu <- rast("sdm/Covariate_currentTime/slope_mean2m.tif")
srad.fu <- rast ("sdm/Covariate_currentTime/srad_mean_060708_summer.tif")
hillshade.fu <- rast ("sdm/Covariate_currentTime/hillshade_mean2m.tif")

# Load and reproject rasters (if needed)
target_crs <- "EPSG:2056"

bio15.fu.rcp45 <- project(rast("sdm/covariates_Future/bio15_20702099_RCP45.tif"), target_crs)
gdd3.fu.rcp45 <- project(rast("sdm/covariates_Future/gdd3_20702099_RCP45.tif"), target_crs)
f.fu          <- project(rast("sdm/Covariate_currentTime/f.tif"), target_crs)
r.fu          <- project(rast("sdm/Covariate_currentTime/r.tif"), target_crs)
forestaggr.fu <- project(rast("sdm/Covariate_currentTime/forestaggr_100.tif"), target_crs)
slope.fu      <- project(rast("sdm/Covariate_currentTime/slope_mean2m.tif"), target_crs)
srad.fu       <- project(rast("sdm/Covariate_currentTime/srad_mean_060708_summer.tif"), target_crs)
hillshade.fu  <- project(rast("sdm/Covariate_currentTime/hillshade_mean2m.tif"), target_crs)



biostack.fut <- c(bio15.fu.rcp45,gdd3.fu.rcp45,f.fu,r.fu, forestaggr.fu, slope.fu, 
                  srad.fu, hillshade.fu )

names(biostack.fut) <- c("bio15", "gdd3", "r", "f","forestaggr","slope" , "srad" ,
                         "hillshade" )

names(biostack.fut)

fungus.fut.gam <- predict(biostack.fut, gam_simplified, type="response")
fungus.fut.glm <- predict(biostack.fut, glm_final, type="response")
fungus.fut.rf <- predict(biostack.fut, rf_final, type="prob", na.rm=T)[[2]]
## Make an ensemble based on a weighted mean
fungus.fut.wmean <- weighted.mean(c(fungus.fut.gam, fungus.fut.glm,
                                    fungus.fut.rf),
                                  w= c(0.6997442 ,0.7613171 , 0.6175831))

# plot the predictions
fungus.df.cur <- cbind(as.data.frame(fungus.curr.wmean, xy=T, na.rm=T),
                       Map="Current")
fungus.df.fut <- cbind(as.data.frame(fungus.fut.wmean, xy=T, na.rm=T),
                       Map="Future")
tmp.dif <- cbind(as.data.frame(fungus.fut.wmean-fungus.curr.wmean,
                               xy=T, na.rm=T),
                 Map="Differences")
df.all <- rbind(fungus.df.cur, fungus.df.fut, tmp.dif)
(future_map_1 <- ggplot(df.all, aes(x=x, y=y, fill=sum))+
  geom_raster(alpha=0.8)+
  coord_equal()+
  facet_wrap(~Map, nrow=1)+
  scale_fill_gradientn("Values", colours=pal1)+
  theme_void())

# actually a good plot!! area reduces in suitability

# fuuuck

# binarizing the meps: abs & p of suitable habitat
# use the threshold that maximizes TSS

predval <- extract(fungus.curr.wmean, species_data[,1:2], ID = F)
thr <- ecospat.max.tss(Pred = predval, Sp.occ = species_data$presence)$max.threshold
thr # 0.31
# above this probability it is already presence??? that is kinda crazy

fungus.curr.bin <- bm_BinaryTransformation(fungus.curr.wmean,threshold = thr)
fungus.fut.bin <- bm_BinaryTransformation(fungus.fut.wmean,threshold = thr)

fungus.df.curr.bin <- cbind(as.data.frame(fungus.curr.bin, xy=T, na.rm=T),
                         Map="Current")
fungus.df.fut.bin <- cbind(as.data.frame(fungus.fut.bin, xy=T, na.rm=T),
                         Map="Future")

df.all <- rbind(fungus.df.curr.bin,fungus.df.fut.bin)
df.all$sum = as.factor(df.all$sum)

# how is the range expected to change?
(range_change <- ggplot(df.all, aes(x=x, y=y, fill=sum))+
  geom_raster(alpha=0.8)+
  coord_equal()+
  scale_fill_discrete("Habitat\npreference",type=c("#8F8F8F","#8BC8AC"))+
  facet_wrap(~Map, nrow=1)+
  theme_void())
# looks kinda similar??? idk

# in numbers
RangeShift <- BIOMOD_RangeSize(proj.current = fungus.curr.bin,proj.future = fungus.fut.bin)
RangeShift$Compt.By.Models
# some lost some gained, rel. balanced out it seems.

#Plotting the output
tmp <- as.data.frame(RangeShift$Diff.By.Pixel,xy=T)
tmp$values <- ""
tmp$values[which(tmp$sum==0)] <- "Abs"
tmp$values[which(tmp$sum==1)] <- "Gain"
tmp$values[which(tmp$sum==-1)] <- "Stable"
tmp$values[which(tmp$sum==-2)] <- "Lost"

(change_plot <- ggplot(tmp, aes(x=x, y=y, fill=values))+
  geom_raster(alpha=0.8)+
  coord_equal()+
  scale_fill_manual("",values=c("grey90","#33a02c","#b2df8a", "orange"))+
  theme_void() )
# change colors, but otherwise seems ok






# trying with different variables ----
# slope; soil properties (ph and moisture); summer air T, growing degree days, and yearly air T most influential

# try 1 temperature and 1 precip var???

# still doing correlations ....----

# i think this will be final attempt, i am so bored
library(dplyr)
library(caret)

# Step 1: Correlation matrix for environmental predictors
cor_matrix <- cor(env_values, use = "pairwise.complete.obs")

# Step 2: Cluster variables using 1 - abs(correlation)
distance_matrix <- as.dist(1 - abs(cor_matrix))
cor_clust <- hclust(distance_matrix, method = "complete")
cut_height <- 0.3  # corresponds to correlation threshold of 0.7
clusters <- cutree(cor_clust, h = cut_height)

# Step 3: For each cluster, fit univariate logistic regressions and pick best predictor
get_best_predictor <- function(varnames) {
  scores <- sapply(varnames, function(var) {
    df <- data.frame(presence = species_data$presence, x = env_values[[var]])
    model <- glm(presence ~ x, data = df, family = binomial)
    null_model <- glm(presence ~ 1, data = df, family = binomial)
    # McFadden's pseudo-R²
    r2 <- 1 - logLik(model) / logLik(null_model)
    return(as.numeric(r2))
  })
  varnames[which.max(scores)]
}

# Step 4: Loop through clusters and select best predictor from each
selected_vars <- sapply(unique(clusters), function(cluster_id) {
  vars_in_cluster <- names(clusters[clusters == cluster_id])
  get_best_predictor(vars_in_cluster)
})

# Final selection
selected_vars <- sort(selected_vars)
selected_vars

# not doing stepwise anymore - but decide whether quadratic terms in the glm----


selected_predictors <- c("bio15", "gdd3", "f", "r", "forestaggr", 
                         "slope", "twi25ss", "srad", "hillshade")

results <- list()
final_terms <- c()

for (var in selected_predictors) {
  
  # build formulas
  f_linear <- as.formula(paste("presence ~", var))
  f_quad   <- as.formula(paste("presence ~", var, "+ I(", var, "^2)"))
  
  # fit models
  mod_linear <- glm(f_linear, data = species_data, family = binomial)
  mod_quad   <- glm(f_quad, data = species_data, family = binomial)
  
  # compare AIC
  aic_linear <- AIC(mod_linear)
  aic_quad <- AIC(mod_quad)
  
  # store result
  results[[var]] <- data.frame(
    predictor = var,
    AIC_linear = aic_linear,
    AIC_quad = aic_quad,
    better = ifelse(aic_quad + 2 < aic_linear, "quadratic", "linear")
  )
  
  # record best-fitting term
  if (aic_quad + 2 < aic_linear) {
    final_terms <- c(final_terms, var, paste0("I(", var, "^2)"))
  } else {
    final_terms <- c(final_terms, var)
  }
}

# Combine and view model comparison results
model_comparison <- do.call(rbind, results)
print(model_comparison)

# Show final formula
cat("Final model formula:\n")
cat("presence ~", paste(final_terms, collapse = " + "), "\n")



# notes 
# remove twiss - no effect and can say we wanted as few variables as possible due to
# not being able to model their changes in the future
# the rest i can justify

# validation --ORIGINAL ----
# unsure if correct bc some things changed compared to tp script

# Ensure the data is correct
fungus <- data.frame(bio15 = species_data$bio15,
                     bio13 = species_data$bio13,
                     r = species_data$r,
                     presence = as.factor(species_data$presence))

set.seed(123)
cv.error.10 <- rep(0, 10)

for (i in 1:10) {
  glm.fit <- glm(presence ~ poly(bio15 + bio13 + r, i),  # polynomial of the sum
                 family = "binomial",
                 data = fungus)
  
  # Use the same dataset used for model fitting
  cv.error.10[i] <- cv.glm(fungus, glm.fit, K = 10)$delta[1]
}

cv.error.10

# adding complexity beyond degree 2 or 3 does not improve accuracy meaningfully.

# example with rf model
set.seed(123)
fungus$presence<-make.names(fungus$presence)
# idk what this does lol
# define training control
train_control <- trainControl(method = "cv", number = 5, savePredictions = T,
                              summaryFunction = twoClassSummary, classProbs = T)
# train model
model <- train(presence~ ., data = fungus, trControl = train_control, method = "rf")

# plot roc
selectedIndices <- model$pred$mtry == 3
ROC <- roc(as.numeric(model$pred$obs[selectedIndices]),
           as.numeric(model$pred$X0[selectedIndices]), auc = TRUE)
confMat <- caret::confusionMatrix(data = model$pred$obs[selectedIndices],
                                  reference = model$pred$pred[selectedIndices],
                                  mode= "everything")
confMat$overall
par(mfrow=c(1,1))
plot.roc(smooth(ROC))


# k-fold cv: for gam, from chatgpt----
library(caret)
library(pROC)

set.seed(123)
Vero_chama$Veronica_chamaedrys <- make.names(Vero_chama$Veronica_chamaedrys)

# Define train control
train_control <- trainControl(method = "cv", number = 10,
                              summaryFunction = twoClassSummary,
                              classProbs = TRUE, savePredictions = TRUE)

# Train the GAM model
model_gam <- train(Veronica_chamaedrys ~ bio6 + bio7 + bio15,
                   data = Vero_chama,
                   method = "gam",
                   trControl = train_control,
                   metric = "ROC",
                   family = "binomial")

# View results
model_gam
library(pROC)

# Extract predictions for a specific resample
preds <- model_gam$pred
# Use only one resample (caret returns all folds merged)
preds_fold <- preds[preds$Resample == "Fold1", ]

roc_obj <- roc(preds_fold$obs, preds_fold$X1)  # X1 = predicted probability for class '1'

plot.roc(roc_obj)
TSS <- roc_obj$sensitivities + roc_obj$specificities - 1
max(TSS)  # best threshold




# trying so hard to get tss from ncv----


bio15.curr <- env_variables$bio15
gdd3.curr <- env_variables$gdd3
summary(values(bio15.curr))
summary(values(bio15.fu.rcp45))

summary(values(gdd3.curr))
summary(values(gdd3.fu.rcp45))


bio15.fu.clipped <- clamp(bio15.fu.rcp45, lower=min(values(bio15.curr), na.rm=TRUE), 
                          upper=max(values(bio15.curr), na.rm=TRUE))
gdd3.fu.clipped <- clamp(gdd3.fu.rcp45, lower=min(values(gdd3.curr), na.rm=TRUE), 
                         upper=max(values(gdd3.curr), na.rm=TRUE))

# Recreate the biostack with clipped rasters
biostack.fut.clipped <- c(bio15.fu.clipped, gdd3.fu.clipped, f.fu, r.fu, 
                          forestaggr.fu, slope.fu, srad.fu, hillshade.fu)

# Set names for the stack
names(biostack.fut.clipped) <- c("bio15", "gdd3", "r", "f", "forestaggr", "slope", "srad", "hillshade")



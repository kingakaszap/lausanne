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

corrplot.mixed(cor(env_values), number.cex= 0.5, tl.cex = 0.5)
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
# maybe change forestaggr - see which is the best predictor (but 500 focal too correlated)
# bio1 and bio13 too correlated
# keep bio13, drop bio1 - mabye max T?
# srad ok

# try to see if there is relationship?? idk

# modelling prep----
species_data <- na.omit(data.frame(coord, presence=pa_data$resp,
                             extract(env_variables, coord)))
head(species_data)
View(species_data)

# try to start modelling ----
# glm -----
# Simple GLM with linear predictors
glm1 <- glm(presence ~ bio13 + bio15 + f +r+ forestaggr.1 + slope,
            data=species_data, family="binomial")
summary(glm1)
str(species_data$presence)
# Simple GLM with quadratic predictors
glm2 <- glm(presence ~ I(bio13^2)+I(bio15^2)+I(f^2)  +I(r^2)  +I(forestaggr.1^2)+I(slope^2),
            data=species_data, family="binomial")
summary(glm2)
# A bit more complex GLM : both linear and quadratic predictors
glm3 <- glm(Veronica_chamaedrys ~ bio7 + bio6 + bio15 +
              I(bio7^2) + I(bio6^2) + I(bio15^2),
            data=spData, family="binomial")
# idk i think this would be too much, too many predictors

# ok moving on ----
# predictions
var.names <- c("bio13", "bio15", "f", "r", "forestaggr.1", "slope")
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
pred.glm3<-predict(glm3,newdata=new.data,type="response")

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

# plotting rspc for all 3 vars
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

resp1.glm<-ggplot(out1.glm,aes(x=Env.val,y=Occ.prob,color=Algorithm,
                               linetype=Algorithm)) +
  geom_line(lwd=1) +
  scale_color_manual(values=c("#007991","#439A86","#BCD8C1"))+
  facet_wrap(~Var.name,ncol=2,scale="free_x")+
  theme_bw()+
  labs(x="Value",y="Occurence probability")+
  theme_bw()
print(resp1.glm)

ggsave("smd/report/plots/glm_curves_initial.png", dpi = 600)

# slope looks unnecessary
# can keep these minus slope for now, but also play around with other variables later
# also maybe bio13 not as important??
# maybe bio1 instead of bio13

# stepwise parameter selection -----

glmStart <- glm(presence~1, data = species_data, family = binomial)
glm.formula<-bm_MakeFormula("presence",
                            species_data[,c("bio13", "bio15","f", "forestaggr.1", "r")],
                            "quadratic",interaction.level=1)
# stepwise selection with AIC
glmModAIC <- stepAIC(glmStart, glm.formula,
                     direction = "both",
                     trace = TRUE,
                     k = 2)

# only keeps bio13 and r

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
RF = randomForest(x = species_data[,c("bio13", "bio15", "r", "f", "forestaggr.1")],
                  y = as.factor(species_data$presence),
                  ntree = 1000,
                  importance =TRUE)
RF.pred = predict(RF, type = "prob")[,2]

importance(RF)
?importance

# most important is bio13 and r
# bio15 no importance

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
resp.rf<-ggplot(out.rf, aes(x=Env.val,y=Occ.prob,color=Algorithm,
                            linetype=Algorithm)) +
  geom_line(size=1)+
  scale_color_manual(values=c("red"))+
  facet_wrap(~Var.name,ncol=2,scale="free_x")+
  theme_bw()+
  labs(x="Value",y="Occurenceprobability")+
  theme_bw()
print(resp.rf)

# try with different bio variables, and with forest aggr ddifferent scales. 
# r definitely needs to stay
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
var.names <- c("bio13","bio15","f", "forestaggr.1", "r")
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

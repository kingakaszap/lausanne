# libraries ----
library(biomod2) # V.4.2-4
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
library(mgcv)
options(warn=-1) ### to remove warnings within the PDF

# load data, select predictors ----

data("ecospat.testData") #a dataset within the package ecospat that contain our test species : Veronica_chamaedrys
spObs <- ecospat.testData[,c("long","lat","Veronica_chamaedrys")]
head(spObs) # overview of the data frame
# it is presence absence data
coord<-spObs[,1:2]

# gather all the bioclimatic raster files
raster_files <- mixedsort(list.files(path = "sdm/data_tp1/covariates",
                                     pattern = "bio",
                                     full.names = TRUE))
# gather all raster layers
bioclim <-rast(raster_files)
cat(crs(bioclim))# cat() allows to clearly read the string from crs()

# change the names of the rasters
names(bioclim) = sub("_.*", "", names(bioclim))
  
# plot the different climatic maps
# plot function in the package terra

# classic:
plot(bioclim, col =rev(viridis::inferno(100, direction = -1)))

# need to project the species data to the coordinate system of the maps
coord <- sf_project(from = st_crs("EPSG:21781"),
                    to =st_crs("EPSG:2056"), pts = coord)
head(coord)

# extract bioclimatic values
bioclimValues <- na.omit(data.frame(extract(bioclim, coord)))
# coord corresponds to species coordinates, bioclim to bioclim data
head(bioclimValues)

dim(bioclimValues)
str(bioclimValues)

# now we have predictor data points corresponding to
# species data points
# (how ??)
# select which predictor to include in the model
# usually try: relevant for species, not correlates

# prelim analysis of correlation between predictors
# general rule: min 10 occurrences for 1 exp variable
# cor threshold: usually 65-75%

# create correlation matrix
df.cor.bioclim <- data.frame(round(cor(bioclimValues), 3))
head(df.cor.bioclim)

corrplot.mixed(cor(bioclimValues), number.cex= 0.5, tl.cex = 0.5)

# select least correlated variables based on threshold of 0.7

colnames(df.cor.bioclim) [-c(findCorrelation(as.matrix(df.cor.bioclim), cutoff =0.70))]
# bio3, bio6, bio15
# isothermality(mean diurnal range/annual range)
# min T of coldest month
# precipitation seasonality (coef of variation)

# create clustering of correlations
# transform correlation matrix into distance matrix

cor.clust <- hclust(as.dist(1-abs(df.cor.bioclim)))
# plot these distances as a tree
par(mar = c(1,4,3,1))
plot(cor.clust, main = "Cluster of the correlations among variables",
     ylab = "Height as 1-abs(cor)", xlab = "")
# threshold value
abline(h = 0.3, lty = 2, col = "red", lwd = 2)

cor(bioclimValues[,c("bio7", "bio6", "bio15")])
# bio7 - annual T range

# modeling ----

# prep data
# get species data and corresponding env data in the same table
spData <- na.omit(data.frame(coord, Veronica_chamaedrys=sp0bs$Veronica_chamaedrys,
                             extract(bioclim, coord)))
head(spData)

# number of occurrences for our sp
sum(spData[,"Veronica_chamaedrys"])

# prevalence 
sum(spData[,"Veronica_chamaedrys"]) /nrow(spData)

dim(spData)

# 1 - environmental envelope technique---
pred_BIOCLIM_100 <- bm_SRE(resp.var = spData$Veronica_chamaedrys,
                           expl.var = spData[, c("bio7", "bio6", "bio15")],
                          new.env = subset(bioclim, c ("bio7", "bio6", "bio15")),
                          quant = 0.025)

pred_BIOCLIM_950<-bm_SRE(resp.var=spData$Veronica_chamaedrys,
                         expl.var= spData[,c("bio7","bio6","bio15")],
                         new.env= subset(bioclim,c("bio7","bio6","bio15")),
                         quant=0.025)

pred_BIOCLIM_900<-bm_SRE(resp.var=spData$Veronica_chamaedrys,
                         expl.var= spData[,c("bio7","bio6","bio15")],
                         new.env= subset(bioclim,c("bio7","bio6","bio15")),
                         quant=0.05)

# plot to compare distributions
par(mfrow = c(2,2))
plot(pred_BIOCLIM_100, col = c("#8F8F8F","#8BC8AC"),main="BIOCLIM100%")
points(spData[spData$Veronica_chamaedrys==1,1:2],pch=16,cex=0.5)
# presence: plain circle, absence: empty circle
points(spData[spData$Veronica_chamaedrys==0,1:2],cex=0.5)
plot(pred_BIOCLIM_950,col=c("#8F8F8F","#8BC8AC"),main="BIOCLIM95%")
points(spData[spData$Veronica_chamaedrys==1,1:2],pch=16,cex=0.5)
points(spData[spData$Veronica_chamaedrys==0,1:2],cex=0.5)
plot(pred_BIOCLIM_900,col=c("#8F8F8F","#8BC8AC"),main="BIOCLIM90%")
points(spData[spData$Veronica_chamaedrys==1,1:2],pch=16,cex=0.5)
points(spData[spData$Veronica_chamaedrys==0,1:2],cex=0.5)

# predictions dont match too well with data: some locations where its predicted as present but no occurence recorded

# regression based approaches---

# Simple GLM with linear predictors
glm1 <- glm(Veronica_chamaedrys ~ bio7 + bio6 + bio15,
            data=spData, family="binomial")
# Simple GLM with quadratic predictors
glm2 <- glm(Veronica_chamaedrys ~ I(bio7^2)+I(bio6^2)+I(bio15^2),
            data=spData, family="binomial")
# A bit more complex GLM : both linear and quadratic predictors
glm3 <- glm(Veronica_chamaedrys ~ bio7 + bio6 + bio15 +
              I(bio7^2) + I(bio6^2) + I(bio15^2),
            data=spData, family="binomial")

# response curves
# species response along one variable gradient
# represents the changes in PROBABILITY to find that species when
# it is going through all values of a predictor gradient
# all other predictors are fixed!
# as predicted by the model
# plot curve for bio7 (T)
var.names <- c("bio7", "bio6", "bio15")
# fix other values
medians<- apply(spData[, var.names],2,median)
medians_table <- data.frame(sapply(medians, function(x)rep(x, 100)))
summary(medians_table)

# create gradient for bio7
foc.var <- spData[,"bio7"]
var.new <- seq(min(foc.var), max(foc.var), length = 100)
# using 100 points
new.data <- medians_table
new.data[,"bio7"] <- var.new # put the gradient in the table of medians
summary(new.data)
head(new.data)

# ask the 3 glm-s to predict on the new ("fake") dataset
pred.glm1 <-predict(glm1,newdata=new.data,type="response")
pred.glm2<-predict(glm2,newdata=new.data,type="response")
pred.glm3<-predict(glm3,newdata=new.data,type="response")

pred.glob <- c(pred.glm1, pred.glm2, pred.glm3)

tmp1<-cbind(Occ.prob=pred.glob,Env.val=var.new)
tmp2 <-data.frame(cbind(Algorithm = rep(c("GLM1", 'GLM2', "GLM3"), each = 100), 
                        Var.name = "bio7"))
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
for(i in 1:3) {
  foc.var <- spData[,var.names[i]]
  new.data <- medians_table
  var.new <- seq(min(foc.var), max(foc.var), length = 100)
  new.data[,i] <- var.new
  pred.glm1 <- predict (glm1, newdata = new.data, type = "response")
  pred.glm2 <- predict (glm2, newdata = new.data, type = "response")
  pred.glm3 <- predict(glm3, newdata = new.data, type = "response")
  pred.glob <- c(pred.glm1, pred.glm2, pred.glm3)
  tmp1 <- cbind(Occ.prob= pred.glob, Env.val=var.new)
  tmp2 <- data.frame(cbind(Algorithm = rep(c("GLM1", "GLM2", "GLM3"), each = 100),
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

# stepwise parameter selection----
# choose the formula (model) that fits each predictor best
glmStart <- glm(Veronica_chamaedrys~1, data = spData, family = binomial)
# most complex
glm.formula<-bm_MakeFormula("Veronica_chamaedrys",
                            spData[,c("bio7", "bio6","bio15")],
                            "quadratic",interaction.level=1)
#Itwillbeusedasapoolinwhichthestepwiseselectionalgorithmcantakepredictors
#toincorporateinthesimplemodel.
glm.formula

# stepwise selection with AIC
glmModAIC <- stepAIC(glmStart, glm.formula,
                     direction = "both",
                     trace = TRUE,
                     k = 2)

# Stepwise selection with BIC
glmModBIC<-stepAIC(glmStart,
                   glm.formula,
                   direction="both",
                   trace=FALSE,
                   k=log(nrow(spData)))

?stepAIC

out2.glm=NULL
for(i in 1:3){
  foc.var<-spData[,var.names[i]]
  new.data<-medians_table
  var.new<-seq(min(foc.var),max(foc.var),length=100)
  new.data[,i]<-var.new
  pred.glm1<-predict(glm1,newdata=new.data,type="response")
  pred.glm2<-predict(glm2,newdata=new.data,type="response")
  pred.glm3<-predict(glm3,newdata=new.data,type="response")
  pred.glmAIC<-predict(glmModAIC,newdata=new.data,type="response")
  pred.glmBIC<-predict(glmModBIC,newdata=new.data,type="response")
  pred.glob <-c(pred.glm1,pred.glm2,pred.glm3,pred.glmAIC,pred.glmBIC)
  tmp1<-cbind(Occ.prob=pred.glob,Env.val=var.new)
  tmp2<-data.frame(cbind(Algorithm=rep(c("GLM1","GLM2","GLM3","glmAIC","glmBIC"),
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

# classification approaches ----
# cluster analysis

# recursive partitioning ----

RP <- rpart(Veronica_chamaedrys~bio7 + bio6 + bio15,
            data = spData,
            control=rpart.control(xval=1000), method = "class")

par(mfrow=c(1,1), mar=c(0,0,0,0))
plot(RP, uniform=F, margin=0.1, branch=1)
text(RP)

# random forests, baggin appr

RF = randomForest(x = spData[,c("bio7", "bio6", "bio15")],
                  y = as.factor(spData$Veronica_chamaedrys),
                  ntree = 1000,
                  importance =TRUE)
RF.pred = predict(RF, type = "prob")[,2]

importance(RF)

# plot response curves
out.rf <- NULL
for (i in 1:3){
  foc.var <- spData[,var.names[i]]
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

# create GAM model using stepwise selection of variables----
?step.Gam
?predict.Gam
library(gam)

# Start with a minimal GAM model including one predictor
gamStart <- gam(Veronica_chamaedrys ~ s(bio15), data = spData, family = binomial)

# Start with a minimal model with all the predictors
gamModel <- gam(Veronica_chamaedrys ~ s(bio15) + s(bio6) + s(bio7), data = spData, family = binomial)

# Perform stepwise selection with proper scope formatting
gamModAIC <- step.Gam(
  gamStart, gamModel,
  direction = "both", 
  trace = TRUE
)

gamModAIC <- step.Gam(gamStart, 
                      scope = list(upper = gamModel), 
                      direction = "both", trace = TRUE)





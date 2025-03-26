library(unmarked)
data(mallard)
#visualize the data structure
mallard.y[1:5, ]
mallard.site[1:5, ]
mallard.obs$ivel[1:5, ]
mallard.obs$date[1:5, ]
#build the unmarked frame
mallardUMF<-unmarkedFramePCount(y=mallard.y, siteCovs=mallard.site,obsCovs=mallard.obs)
#summaries
summary(mallardUMF)
names(mallard.site)
#build a few models …
(null <- pcount(~1~1, mallardUMF))
(elev<- pcount(~1~elev, mallardUMF))
(date<- pcount(~date~1, mallardUMF))
(forest <- pcount(~1~forest, mallardUMF))
(global <- pcount(~date ~elev+forest, mallardUMF))
#...and compare them
fits1 <- fitList(null=null, elev=elev, date=date, forest=forest, global=global)
(ms1 <- modSel(fits1, nullmod='null'))
# build a more elaborate model
fm.mallard.1<-pcount(~date+I(date^2)~elev+forest,data = mallardUMF, K=50)
summary(fm.mallard.1)
#refine the model
fm.mallard.2<-pcount(~date~elev+forest,data = mallardUMF, K=50)
summary(fm.mallard.2)
#plot the model
par(mfrow = c(1,3))
beta1 <- coef(fm.mallard.2)
plot(function(x) plogis(beta1[1] + beta1[2]*x), 1, 3, xlab="Elevation", ylab="Occupancy probability", ylim = c(0, 1))
plot(function(x) plogis(beta1[1] + beta1[3]*x), -3, 3, xlab="Forest cover", ylab="Occupancy probability", ylim = c(0, 1))
plot(function(x) plogis(beta1[4] + beta1[5]*x), xlab="Date", ylab="Detection probability", ylim = c(0, 1))
#estimate abundance
ranef(fm.mallard.2)
sum(ranef(fm.mallard.2)@post[,2,])

#Data for your exercise

minkFrogData<-read.csv(file.choose(), header=T) 
str(minkFrogData)
minkFrog <- unmarkedFramePCount(y = minkFrogData[,1:5], siteCovs = minkFrogData[,6:8],obsCovs = list(precip = minkFrogData[,9:13]))
summary(minkFrog)
plot(minkFrog)

#Questions:
#Try building several modelswith different detectability and habitat quality covariates (using sensible biological assumptions). What is the best model you can build (lowest AIC)?
#How does precipitation during the surveys affect detectability?
#How does distance to road affect abundance?
#What's the total abundance of Mink Frog in the study system?
# libraries and data----
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

if( "package:gam" %in% search()){
detach("package:gam", unload = TRUE)
 }

data("ecospat.testData")
spObs <- ecospat.testData[,c("long","lat","Veronica_chamaedrys")]

coord <- spObs[,1:2]
#Project the coordinates in the CRS of climate maps
coord <- sf_project(from =st_crs("EPSG:21781"),to= st_crs("EPSG:2056"),pts=coord)

# Load and extract bioclim data
raster_files <- mixedsort(list.files("../data/covariates",pattern="bio",full.names = T))
bioclim <- rast(raster_files)
#Change the names of the rasters
names(bioclim) = sub("_.*","",names(bioclim))
#Remove coordinates that have no climatic information
spData <- na.omit(data.frame(X=coord[,1], Y=coord[,2],
                             Veronica_chamaedrys=spObs$Veronica_chamaedrys,
                             extract(bioclim, coord)))

gam1 <- mgcv::gam(Veronica_chamaedrys ~ s(bio6,k=4) + s(bio7,k=4) + s(bio15,k=4),
                  data=spData, family="binomial")
glm1 <- glm(formula = Veronica_chamaedrys ~ bio7 + I(bio7^2) + I(bio6^2),
            family = binomial, data = spData)
rf1 = randomForest(x = spData[,c("bio7", "bio6", "bio15")],
                   y = as.factor(spData$Veronica_chamaedrys),
                   ntree = 1000)

# predicting in the calibration area----
#For computation time we are going to decrease the resolution of climate data
# For the report, do not do that (except if you are working with insects)

bioclim <- aggregate(bioclim,fact=2,fun="mean",na.rm=T)

# Predicted occurence probabilities
vero.curr.gam <- predict(bioclim, gam1, type="response")
tmp.vero.gam.df <- cbind.data.frame(Model = "GAM",
                                    as.data.frame(vero.curr.gam,
                                                  xy=T, na.rm=T))
vero.curr.glm <- predict(bioclim, glm1, type="response")
tmp.vero.glm.df <- cbind.data.frame(Model = "GLM",
                                    as.data.frame(vero.curr.glm,
                                                  xy=T, na.rm=T))
vero.curr.glm <- predict(bioclim, glm1, type="response")
tmp.vero.glm.df <- cbind.data.frame(Model = "GLM",
                                    as.data.frame(vero.curr.glm,
                                                  xy=T, na.rm=T))
vero.curr.rf <- predict(bioclim, rf1, type="prob", na.rm=T)[[2]]
tmp.vero.rf.df <- cbind.data.frame(Model = "RF",
                                   as.data.frame(vero.curr.rf,
                                                 xy=T, na.rm=T))
colnames(tmp.vero.rf.df) = c("Model","x","y","lyr1")

# dataset for plotting

PlotData <- rbind.data.frame(tmp.vero.gam.df,tmp.vero.glm.df,
                             tmp.vero.rf.df)

# plot the new spatial patterns
pal1 <- rev(colorRampPalette(brewer.pal(11,"RdYlGn")[11:1])(64))
ggplot(PlotData, aes(x=x, y=y, fill=lyr1)) +
  geom_raster(alpha=0.8) +
  facet_wrap(~Model) +
  coord_equal()+
  scale_fill_gradientn("Occurence\nprobability",colours=pal1)+
  theme_void()

# predict SE for GAM (also available for GLM) for the df and plot the predictions

vero.se <- predict(bioclim,gam1, type="response", se.fit=TRUE)
df.gam <- as.data.frame(vero.se, xy=T, na.rm=T)
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
vero.curr <- c(vero.curr.gam, vero.curr.glm,vero.curr.rf)

vero.curr.mean <- mean(vero.curr)
tmp.vero.mean.df <- cbind.data.frame(Ens = "Mean",
                                     as.data.frame(vero.curr.mean,
                                                   xy=T, na.rm=T))
colnames(tmp.vero.mean.df) = c( "Ens","x", "y", "Prob")

# Standard deviation
vero.curr.sd <- stdev(vero.curr)
tmp.vero.sd.df <- cbind.data.frame(Ens = "sd",
                                   as.data.frame(vero.curr.sd,
                                                 xy=T, na.rm=T))
colnames(tmp.vero.sd.df) = c( "Ens","x", "y", "Prob")
# Weighted Mean
#Here a weighted mean is perform based on the mean AUC value
#obtained during practical 2 based on 20 CV (can slightly change)
vero.curr.wmean <- weighted.mean(vero.curr,
                                 w= c(0.84,0.84, 0.69))
tmp.vero.wmean.df <- cbind.data.frame(Ens = "wMean",
                                      as.data.frame(vero.curr.wmean,
                                                    xy=T, na.rm=T))
colnames(tmp.vero.wmean.df) = c( "Ens","x", "y", "Prob")

# make dataset for ggplot
PlotData <- rbind.data.frame(tmp.vero.mean.df,
                             tmp.vero.sd.df,
                             tmp.vero.wmean.df)
# Plot the spatial patterns
pal1 <- rev(colorRampPalette(brewer.pal(11,"RdYlGn")[11:1])(64))
ggplot(PlotData, aes(x=x, y=y, fill=Prob)) +
  geom_raster(alpha=0.8) +
  facet_wrap(~Ens) +
  coord_equal()+
  scale_fill_gradientn("Occurence\nprobability",colours=pal1)+
  theme_void()

# predict in a new area
bioclimData <- as.data.frame(bioclim,xy=T)
tmp1 <- dudi.pca(bioclimData[,c("bio6", "bio7" , "bio15")],
                 nf=2, scannf=F)
tmp2 <- data.frame(cbind(bioclimData, tmp1$li))
tmp2$ID <- ifelse(bioclimData$x<2570000, 0, 1)
map.plot <- ggplot(tmp2, aes(x=x, y=y,fill=factor(ID))) +
  geom_raster(alpha=1)+
  coord_equal()+
  scale_fill_manual(values=c("#3090C733","#9F000F4D"),
                    labels = c("Valleys","Mountains"),
                    name="")+
  theme_void()+
  guides(colour=guide_legend(override.aes=list(size=2)))
# Environmental conditions in the two worlds
pca.plot <- ggplot(tmp2,aes(x=Axis1, y=Axis2,fill=factor(ID)))+ stat_bin2d(alpha=0.5)+
  theme_bw()+
  labs(colour="",x="PCA-Axis 1", y="PCA-Axis 2")+
  scale_fill_manual(values=c("#3090C733","#9F000F4D"),
                    labels = c("Valleys","Mountains"),
                    name="")
plot_grid(map.plot, pca.plot, labels=c("A", "B"), ncol=1)

Valleys <- bioclimData[bioclimData$x<2570000,]

mess.Vero <- ecospat.mess(bioclimData[,c("x", "y",
                                         "bio6", "bio7", "bio15")],
                          Valleys[,c("x", "y",
                                     "bio6", "bio7", "bio15")])
mess.Vero <- as.data.frame(mess.Vero)
# Plot the MESS outputs

pal3 <- rev(colorRampPalette(brewer.pal(11,"RdYlGn"))(64))
pal4 <- colorRampPalette(brewer.pal(11,"RdYlGn"))(64)
mess <- ggplot(data=mess.Vero, aes(x=x, y=y, fill=MESS)) +
  geom_raster() +
  scale_fill_gradientn(colours=pal4)+
  coord_equal()+
  theme_void()
mess.w <- ggplot(data=mess.Vero, aes(x=x, y=y, fill=MESSw)) +
  geom_raster() +
  scale_fill_gradientn(colours=pal4)+
  coord_equal()+
  theme_void()
mess.neg <- ggplot(data=mess.Vero, aes(x=x, y=y, fill=MESSneg)) +
  geom_raster() +
  scale_fill_gradientn(colours=pal3)+
  coord_equal()+
  theme_void()
plot_grid(mess, mess.w, mess.neg, labels=c("A", "B", "C"), nrow=1)

# influence of grain size (resolution) on environment
# Load the yearly average precipitation data
tave0.75km <- rast("sdm/data_tp3/bioclim/current/Europe_tave.tif")
# Average information - aggregare to coarser grain sizes
tave3km <- aggregate(tave0.75km, fact=4, fun=mean,na.rm=T)

tave7.5km <- aggregate(tave0.75km, fact=10, fun=mean,na.rm=T)
tave30km <- aggregate(tave7.5km, fact=4, fun=mean,na.rm=T)
tave75km <- aggregate(tave7.5km, fact=10, fun=mean,na.rm=T)
tave300km <- aggregate(tave75km, fact=6, fun=mean,na.rm=T)
tave750km <- aggregate(tave75km, fact=10, fun=mean,na.rm=T)

### Classical plots - finer patterns disappear
par(mfrow=c(2,2))
plot(tave0.75km, main="Mean Temperature 0.75km")
plot(tave7.5km, main="Mean Temperature 7.5km")
plot(tave75km, main="Mean Temperature 75km")
plot(tave750km, main="Mean Temperature 750km")

# calculating some cell stats
sss <- data.frame(matrix(NA, nrow=7, ncol=4))
names(sss) <- c("GrainSize","min","mean","max")
sss[,1] <- c(0.75,3,7.5,30,75,300,750)
row.names(sss) <- c("0.75km","3km","7.5km","30km","75km","300km","750km")
# Minimum values
sss$min <- c(as.numeric(global(tave0.75km, fun="min",na.rm=T)),
             as.numeric(global(tave3km, fun="min",na.rm=T)),
             as.numeric(global(tave7.5km, fun="min",na.rm=T)),
             as.numeric(global(tave30km, fun="min",na.rm=T)),
             as.numeric(global(tave75km, fun="min",na.rm=T)),
             as.numeric(global(tave300km, fun="min",na.rm=T)),
             as.numeric(global(tave750km, fun="min",na.rm=T)))
print(sss$min) # big difs with dif grain size


# Average values
sss$mean <- c(as.numeric(global(tave0.75km, fun="mean",na.rm=T)),
              as.numeric(global(tave3km, fun="mean",na.rm=T)),
              as.numeric(global(tave7.5km, fun="mean",na.rm=T)),
              as.numeric(global(tave30km, fun="mean",na.rm=T)),
              as.numeric(global(tave75km, fun="mean",na.rm=T)),
              as.numeric(global(tave300km, fun="mean",na.rm=T)),
              as.numeric(global(tave750km, fun="mean",na.rm=T)))
print(sss$mean)

# Maximum values
sss$max <- c(as.numeric(global(tave0.75km, fun="max",na.rm=T)),
             as.numeric(global(tave3km, fun="max",na.rm=T)),
             as.numeric(global(tave7.5km, fun="max",na.rm=T)),
             as.numeric(global(tave30km, fun="max",na.rm=T)),
             as.numeric(global(tave75km, fun="max",na.rm=T)),
             as.numeric(global(tave300km, fun="max",na.rm=T)),
             as.numeric(global(tave750km, fun="max",na.rm=T)))
print(sss$max)

# plot degradation info
tmp.df <- melt(sss, id="GrainSize")
ggplot(tmp.df, aes(x=GrainSize, y=value, color=variable))+
  geom_line()+
  theme_bw()+
  labs(y="Temperature: min/max/mean within area")
# min temp increases, coldest spots get averaged away
# max temp decreases, hottest spots smoothed
# mean stays stable

# projecting in time----
# used rcp 85, also try rcp 45
bio6r.fu <- rast("../data/covariates/covariates_future/bio6_20702099_RCP45.tif")
bio7r.fu <- rast("../data/covariates/covariates_future/bio7_20702099_RCP45.tif")
bio15r.fu <- rast("../data/covariates/covariates_future/bio15_20702099_RCP45.tif")
biostack.fut <- c(bio6r.fu,bio7r.fu,bio15r.fu)
names(biostack.fut) <- c(paste0("bio",c(6,7,15)))

names(biostack.fut)

#For computation time, we are going to decrease the resolution of climate data
## Not to do for your report except if you selected insects
biostack.fut <- aggregate(biostack.fut,2,"mean",na.rm=T)


vero.fut.gam <- predict(biostack.fut, gam1, type="response")
vero.fut.glm <- predict(biostack.fut, glm1, type="response")
vero.fut.rf <- predict(biostack.fut, rf1, type="prob", na.rm=T)[[2]]
## Make an ensemble based on a weighted mean
vero.fut.wmean <- weighted.mean(c(vero.fut.gam, vero.fut.glm,
                                  vero.fut.rf),
                                w= c(0.84,0.84, 0.69))


# plot the predictions
vero.df.cur <- cbind(as.data.frame(vero.curr.wmean, xy=T, na.rm=T),
                     Map="Current")
vero.df.fut <- cbind(as.data.frame(vero.fut.wmean, xy=T, na.rm=T),
                     Map="Future")
tmp.dif <- cbind(as.data.frame(vero.fut.wmean-vero.curr.wmean,
                               xy=T, na.rm=T),
                 Map="Differences")
df.all <- rbind(vero.df.cur, vero.df.fut, tmp.dif)
ggplot(df.all, aes(x=x, y=y, fill=sum))+
  geom_raster(alpha=0.8)+
  coord_equal()+
  facet_wrap(~Map, nrow=1)+
  scale_fill_gradientn("Values", colours=pal1)+
  theme_void()

# actually a good plot!! area reduces in suitability


PredVal <- extract(vero.curr.wmean,spData[,1:2],ID=F)
thr <- ecospat.max.tss(Pred = PredVal, Sp.occ = spData$Veronica_chamaedrys)$max.threshold
thr

fungus.curr.bin <- bm_BinaryTransformation(vero.curr.wmean,threshold = thr)
vero.fut.bin <- bm_BinaryTransformation(vero.fut.wmean,threshold = thr)

vero.df.cur.bin <- cbind(as.data.frame(vero.curr.bin, xy=T, na.rm=T),
                         Map="Current")
vero.df.fut.bin <- cbind(as.data.frame(vero.fut.bin, xy=T, na.rm=T),
                         Map="Future")
df.all <- rbind(vero.df.cur.bin,vero.df.fut.bin)
df.all$sum = as.factor(df.all$sum)
ggplot(df.all, aes(x=x, y=y, fill=sum))+
  geom_raster(alpha=0.8)+
  coord_equal()+
  scale_fill_discrete("Habitat\npreference",type=c("#8F8F8F","#8BC8AC"))+
  facet_wrap(~Map, nrow=1)+
  theme_void()

RangeShif <- BIOMOD_RangeSize(proj.current = vero.curr.bin,proj.future = vero.fut.bin)

RangeShif$Compt.By.Models

#Plotting the output
tmp <- as.data.frame(RangeShif$Diff.By.Pixel,xy=T)
tmp$values <- ""
tmp$values[which(tmp$sum==0)] <- "Abs"
tmp$values[which(tmp$sum==1)] <- "Gain"
tmp$values[which(tmp$sum==-1)] <- "Stable"
tmp$values[which(tmp$sum==-2)] <- "Lost"
ggplot(tmp, aes(x=x, y=y, fill=values))+
  geom_raster(alpha=0.8)+
  coord_equal()+
  scale_fill_manual("",values=c("grey90","#33a02c","#b2df8a"))+
  theme_void() #here only 3 colours bc no lost


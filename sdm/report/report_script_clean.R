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
# choosing variables based on the cor matrix etc -----
df.cor.envdata <- data.frame(round(cor(env_values), 3))
head(df.cor.envdata)
head(env_values)
# remove id column - hope this doesnt cause issues idk
env_values <- env_values %>% dplyr::select(-ID)
View(env_values)
corrplot.mixed(cor(env_values), number.cex= 0.5, tl.cex = 0.5, tl.pos = "lt")
# this looks wrong lol

# and now is the point where I can select the variables, in theory

nrow(pa_data)

cor(env_values[,c("bio1", "bio15", "forestaggr", "slope", "r", "f", "srad", "hillshade")])

# correlation clustering -----
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
# modelling prep----
species_data <- na.omit(data.frame(coord, presence=pa_data$resp,
                                   extract(env_variables, coord)))

# glm -----
# decided not to include quadratic predictors because it
# caused an inflation of the variance inflation factor,
# and since we do ensembling, nonlinear effects might better be captured in the gam.
glm_final <- glm(presence ~ bio15 + gdd3 + f + 
                   r + forestaggr + slope + srad + hillshade,
                 data = species_data, family = binomial)
summary(glm_final)
vif(glm_final)

# predictions of the glm -----
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

ggsave("sdm/report/plots/responses_glm_noquard.png", dpi = 600)

# gam & preidctions----
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

# rf----
rf_final = randomForest(x = species_data[,c("bio15", "f", "r", "forestaggr", "gdd3", "hillshade",  "slope", "srad")],
                        y = as.factor(species_data$presence),
                        ntree = 1000,
                        importance =TRUE)
RF.pred = predict(rf_final, type = "prob")[,2]

importance(rf_final)

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
# more libraries ----
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
# final models, predictions and a new dataset with these----
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


# calibration - extent to which presence is correctly predicted ----

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
# disrcimination metrics -----
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
(plots_kappa_tss<- ggplot(tmp,aes(x=Threshold,y=Value,color=Model))+
  geom_line(linewidth=1.1)+
  scale_color_discrete(type=c("#007991","#439A86","#BCD8C1"))+
  theme_bw()+
  facet_wrap(~Evaluation.Metric))

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
# a better threshold because not affected by prevalence.

tss.max <- ecospat.max.tss(Pred = EvalData$RF, Sp.occ = EvalData$ObsPA)
tss.max$max.TSS
# gam: 0.5415154; glm: 0.4791289; RF: 0.4314882
# a better threshold bc not affected by prevalence.
head(tss.max$table)
# ncv 20-----
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
  
  eval <- species_data[a$validation,]
  ### GLM ###
  # taken out stepwise
  # Use your final GLM model here instead of the stepwise model
  glm_model <- glm(presence ~ bio15 + gdd3  + f + r + forestaggr +
                     slope + srad + hillshade, 
                   data = calib, family = binomial,
                   control = glm.control(maxit = 100))
  
  # prediction on the evaluation data and evaluation using the AUC approach
  Pred_test <- predict(glm_model, eval, type = "response")
  Test_results["GLM", i] <- as.numeric(auc(roc(eval$presence, Pred_test)))
  
  # prediction on the total dataset
  Pred_results[,"GLM",i] <- predict(glm_model, species_data, type="response")
  ### GAM ###
  gam_mgcv <- mgcv::gam(presence ~ s(bio15) + s(gdd3) + s(r)+ s(f)+ s(forestaggr)
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
# predicted occurence probabilities -----
# simplify gam and pray to god it helps
gam_simplified <- mgcv::gam(presence ~  s(bio15, k = 4) + s(gdd3, k = 4) + s(f, k = 4) + s(r, k=4) + s(forestaggr, k=4) + 
                              s(slope, k=4) + s(srad, k=4) + s(hillshade, k=4),
                            data=species_data, family="binomial")
summary(gam_simplified)
rf_final
glm_final

names(env_variables) <- make.unique(names(env_variables))

fungus.curr.gam <- predict(env_variables, gam_simplified, type="response")
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


# ensemble, present projections ----

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

ggsave("sdm/report/plots/ensemble_maps_2_weird.png", dpi= 600)

# future - straight from ensemble -----
bio15.fu <- rast("sdm/Covariates_Fut/bio15_20702099_RCP45.tif")
gdd3.fu <- rast("sdm/Covariates_Fut/gdd3_20702099_RCP45.tif")
f.fu <- rast("sdm/Covariate_currentTime/f.tif")
r.fu <- rast("sdm/Covariate_currentTime/r.tif")
forestaggr.fu <- rast("sdm/Covariate_currentTime/forestaggr_100.tif")
slope.fu <- rast("sdm/Covariate_currentTime/slope_mean2m.tif")
srad.fu <- rast("sdm/Covariate_currentTime/srad_mean_060708_summer.tif")
hillshade.fu  <- rast("sdm/Covariate_currentTime/hillshade_mean2m.tif")

bio15.fu.85 <- rast("sdm/Covariates_Fut/bio15_20702099_RCP85.tif")
gdd3.fu.85 <- rast("sdm/Covariates_Fut/gdd3_20702099_RCP85.tif")

biostack.fut <- c(bio15.fu.85 ,
                  gdd3.fu.85 ,
                  f.fu ,
                  r.fu ,
                  forestaggr.fu ,
                  slope.fu ,
                  srad.fu,
                  hillshade.fu  )
names(biostack.fut) <- c("bio15", "gdd3", "f", "r","forestaggr","slope" , "srad" ,
                         "hillshade")

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
    scale_fill_gradientn(name = "Predicted probability of occurrence",  colours = pal1) +
    theme_void())

ggsave("sdm/report/plots/future_probabilities.png", dpi = 600)

# I will cry
future_map_1


#  extracts the predicted probabilities from the raster 
# at the species sampling locations 
# so that you can compare those predictions 
# to the observed presences/absences and find the best threshold 
# for converting probabilities into binary presence/absence predictions.

predval <- extract(fungus.curr.wmean, species_data[,1:2], ID = F)
thr <- ecospat.max.tss(Pred = predval, Sp.occ = species_data$presence)$max.threshold
thr # 0.26
# above this probability it is already presence??? that is kinda crazy

# binary transformation and maps----

fungus.curr.bin <- bm_BinaryTransformation(fungus.curr.wmean,threshold = thr)
fungus.fut.bin <- bm_BinaryTransformation(fungus.fut.wmean,threshold = thr)

fungus.df.curr.bin <- cbind(as.data.frame(fungus.curr.bin, xy=T, na.rm=T),
                            Map="Current")
fungus.df.fut.bin <- cbind(as.data.frame(fungus.fut.bin, xy=T, na.rm=T),
                           Map="Future")

df.all <- rbind(fungus.df.curr.bin,fungus.df.fut.bin)
df.all$sum = as.factor(df.all$sum)

# how is the range expected to change?
(range_change <- ggplot(df.all, aes(x = x, y = y, fill = sum)) +
  geom_raster(alpha = 0.8) +
  coord_equal() +
  scale_fill_discrete("Habitat\npreference", type = c("#E0E0E0", "#1B9E77")) +
  facet_wrap(~Map, nrow = 1) +
  theme_void() +
  theme(
    strip.text = element_text(size = 16, face = "bold") ))


# in numbers
RangeShift <- BIOMOD_RangeSize(proj.current = fungus.curr.bin,proj.future = fungus.fut.bin)
RangeShift$Compt.By.Models
# WHOOPS

#Plotting the output
tmp <- as.data.frame(RangeShift$Diff.By.Pixel,xy=T)
tmp$values <- ""
tmp$values[which(tmp$sum==0)] <- "Absent"
tmp$values[which(tmp$sum==1)] <- "Gain"
tmp$values[which(tmp$sum==-1)] <- "Stable"
tmp$values[which(tmp$sum==-2)] <- "Lost"

(change_plot <- ggplot(tmp, aes(x=x, y=y, fill=values))+
    geom_raster(alpha=0.8)+
    coord_equal()+
    scale_fill_manual("",values=c("grey90","#D95F02","#1B9E77"))+
    theme_void() )


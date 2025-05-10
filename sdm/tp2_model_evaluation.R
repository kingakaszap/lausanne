# think: which data to use for training and which for validation?

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

ecospat.meva.table(Pred = EvalData$GLM, Sp.occ = EvalData$ObsPA, th = 0.5)

kappa.max <- ecospat.max.kappa(Pred = EvalData$GLM, Sp.occ = EvalData$ObsPA)
head(kappa.max$table)

kappa.max$max.Kappa

tss.max <- ecospat.max.tss(Pred = EvalData$GLM, Sp.occ = EvalData$ObsPA)
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
  facet_wrap(~Evaluation.Metric)

# error threshold plots
data <- EvalData[1:5]
N.models <- ncol(data)- 2
par(oma=c(0,0,0,0), mar=c(4,1,1.5,1), mfrow=c(2,2), cex=0.7, cex.lab=1.4, mgp=c(2, 0.5,0))
for (mod in 1:N.models){
  error.threshold.plot(data,
                       which.model = mod,
                       color = TRUE,
                       add.legend = TRUE,
                       legend.cex = 0.7)
}

# roc and auc
# area under the curve of the roc plot
par(mfrow  = c(1,1))
auc.roc.plot(data, color=c("#007991", "#439A86", "#BCD8C1"), legend.cex=0.7, main="")
# should be towards the top left corner

# calculate optimal threshold using different methods
optimal.thresholds(EvalData, opt.methods = 1:12, req.sens=0.9, req.spec = 0.9, FPC = 2, FNC = 1)


# 2.1.3. presence only data
obs <- EvalData$GLM [which(EvalData$ObsPA ==1)]
boyceplot <- ecospat.boyce(fit = EvalData$GLM, obs, nclass = 0,
                           window.w = "default", res = 100, PEplot = T)
boyceplot$cor
abline(a=0,b=max(boyceplot$F.ratio))
 
# validation and stuff - i didnt understand much.. ----
Vero_chama<-data.frame(bio7=spData$bio7,
                       bio6=spData$bio6,
                       bio15=spData$bio15,
                       Veronica_chamaedrys=spData$Veronica_chamaedrys)
Vero_chama$Veronica_chamaedrys<-as.factor(Vero_chama$Veronica_chamaedrys)

# K FOLD CROSS VALIDATION example w glm - what degree polynomial to choose

set.seed(123)
cv.error.10<-rep(0,10)
for(i in 1:10){
  glm.fit<-glm(Veronica_chamaedrys ~poly(bio7+ bio6 + bio15,i),
               family="binomial",
               data=Vero_chama)
  cv.error.10[i]<-cv.glm(spData, glm.fit,K=10)$delta[1]
}
cv.error.10
# adding complexity beyond degree 2 or 3 does not improve accuracy meaningfully.

# example with rf model
set.seed(123)
Vero_chama$Veronica_chamaedrys<-make.names(Vero_chama$Veronica_chamaedrys)
                      # idk what this does lol
# define training control
train_control <- trainControl(method = "cv", number = 5, savePredictions = T,
                              summaryFunction = twoClassSummary, classProbs = T)
# train model
model <- train(Veronica_chamaedrys~ ., data = Vero_chama, trControl = train_control, method = "rf")

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


# other data partitioning methods ----
Vero_chama$Veronica_chamaedrys <- make.names(Vero_chama$Veronica_chamaedrys)

#-- Usual bootstrap
train_control <- trainControl(method="boot", number=5, savePredictions=T,
                              summaryFunction=twoClassSummary, classProbs=T)
model.boot <- train(Veronica_chamaedrys ~ ., data=Vero_chama, trControl=train_control,
                    method="rf")
ROC.boot <- roc(as.numeric(model.boot$pred$obs[selectedIndices]),
                as.numeric(model.boot$pred$X0[selectedIndices]), auc=TRUE)

#-- the 0.632 bootstrap estimator (Efron 1983)
train_control <- trainControl(method="boot632", number=5, savePredictions=T,
                              summaryFunction=twoClassSummary,classProbs=T)
model.boot632 <- train(Veronica_chamaedrys ~ ., data=Vero_chama,
                       trControl=train_control, method="rf")
ROC.boot632 <- roc(as.numeric(model.boot632$pred$obs[selectedIndices]),
                   as.numeric(model.boot632$pred$X0[selectedIndices]), auc=TRUE)
#-- optimism bootstrap estimator (Efron and Tibshirani, 1994)
train_control <- trainControl(method="optimism_boot", number=5, savePredictions=T,
                              summaryFunction=twoClassSummary,classProbs=T)
model.boot.optim <- train(Veronica_chamaedrys ~ ., data=Vero_chama,
                          trControl=train_control, method="rf")
ROC.boot.optim <- roc(as.numeric(model.boot.optim$pred$obs[selectedIndices]),
                      as.numeric(model.boot.optim$pred$X0[selectedIndices]),auc=TRUE)
#-- 5-fold cross validation
train_control <- trainControl(method="cv", number=5, savePredictions=T,
                              summaryFunction=twoClassSummary, classProbs=T)
model.cv <- train(Veronica_chamaedrys ~ ., data=Vero_chama, trControl=train_control,
                  method="rf")
ROC.cv <- roc(as.numeric(model.cv$pred$obs[selectedIndices]),
              as.numeric(model.cv$pred$X0[selectedIndices]), auc=TRUE)
#-- repeated split-sample CV where 75% of the data is used for calibrating and 25% for the evaluation
# You can change the ratio with the parameter p
train_control <- trainControl(method="LGOCV", number=5, savePredictions=T,
                              summaryFunction=twoClassSummary, classProbs=T, p= 0.75)
model.LGOCV <- train(Veronica_chamaedrys ~ ., data=Vero_chama, trControl=train_control,
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

# comparison of different algorithms----
nCV <- 20 # number of cross validation
nRow <- nrow(spData)
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
  a <- bm_SampleBinaryVector(obs = spData$Veronica_chamaedrys, ratio=0.7)
  calib <- spData[a$calibration,]
  eval <- spData[a$validation,]
  ### GLM ###
  glmStart <- glm(Veronica_chamaedrys~1, data=calib, family=binomial)
  glm.formula <- bm_MakeFormula("Veronica_chamaedrys",
                                spData[,c("bio7", "bio6", "bio15")],
                                "quadratic", interaction.level=1)
  glmModAIC <- stepAIC(glmStart, glm.formula, data = calib,
                       direction = "both", trace = FALSE, k = 2,
                       control=glm.control(maxit=100))
  # prediction on the evaluation data and evaluation using the AUC approach
  Pred_test <- predict(glmModAIC, eval, type="response")
  Test_results["GLM",i] <- as.numeric(auc(roc(eval$Veronica_chamaedrys,Pred_test)))
  # prediction on the total dataset
  Pred_results[,"GLM",i] <- predict(glmModAIC, spData, type="response")
  ### GAM ###
  gam_mgcv <- gam(Veronica_chamaedrys ~ s(bio7) + s(bio6) + s(bio15),
                  data=calib, family="binomial")
  # prediction on the evaluation data and evaluation using the AUC approach
  Pred_test <- predict(gam_mgcv, eval, type="response")
  Test_results["GAM",i] <- as.numeric(auc(roc(eval$Veronica_chamaedrys,Pred_test)))
  # prediction on the total dataset
  Pred_results[,"GAM",i] <- predict(gam_mgcv, spData, type="response")
  ### RF ###
  RF_mod = randomForest(x = calib[,c("bio7", "bio6", "bio15")],
                        y = as.factor(calib$Veronica_chamaedrys),
                        ntree = 500, importance = TRUE)
  # prediction on the evaluation data and evaluation using the AUC approach
  Pred_test <- predict(RF_mod, eval, type="response")
  Test_results["RF",i] <- as.numeric(auc(roc(eval$Veronica_chamaedrys,
                                             as.numeric(as.character(Pred_test)))))
  # prediction on the total dataset
  Pred_results[,"RF",i] = predict(RF_mod, spData, type="prob")[,2]
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

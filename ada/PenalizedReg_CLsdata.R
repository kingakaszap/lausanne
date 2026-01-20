#-------------------------------------------------------------------------------
# EXAMPLE: CANCER DATA
#-------------------------------------------------------------------------------

# Load the data

# load(file="/Users/dsesia/Desktop/Teaching/Penalized regression/data_cancer_dep_CTNNB1.RData")
load("C:/Users/Kinga/Downloads/data_cancer_dep_CTNNB1.RData")

data[1:10,1:10]
# each row is a cell line
# 1st column: dependency score for this gene,
# how essential is the gene for survival.
# - : knocking it out reduces cell viability
# + - knocking it out might increase growth or negligible
# other columns: binary mutation indications for other genes
# ie is this gene mutated in the given cell line
# e.g. which mutations associated with a high dependence on ctnnb1?
y=data$CTNNB1_dep_score
X=as.matrix(data[,-1])
dim(X)

hist(data$CTNNB1_dep_score,breaks = 30,main = "CTNNB1 dependency scores",xlab = "Score",col="#FDCF76")
# a few cell lines have very high dependency, ie function very badly without this gene
#-------------------------------------------------------------------------------
# Linear regression using a subset of the observation (so that n<p or n~p)
#-------------------------------------------------------------------------------

lr=lm(CTNNB1_dep_score~.,data=data[1:100,]) # first 100 rows 
summary(lr)
# NA-s --- these predictors might be completely correlated with others, or have all 0s
# a lot of NA-s. like 1 data point and infinite lines going through it
# unable to estimate all parameters, bc not enough data points

#-------------------------------------------------------------------------------
# Training and Test set --- TEST error
#-------------------------------------------------------------------------------

n=nrow(data)
n.train=floor(0.75*n) # 75% of the observations
set.seed(123)
train=sample(1:n, size = n.train, replace = FALSE)
X.train=X[train,]
y.train=y[train]
X.test=X[-train,]
y.test=y[-train]

#-------------------------------------------------------------------------------
# Standard linear regression (OLS estimates)
#-------------------------------------------------------------------------------

lr=lm(CTNNB1_dep_score~.,data=data[train,])
summary(lr) # mse 0.1889
n.var.lr=length(which(summary(lr)$coeff[-1,1]!=0)) 
# number of variables (excluding the intercept) different from 0 (for summary table)

lr.pred=predict(lr,newdata = data[-train,])
test_error_ols=mean((lr.pred-y.test)^2)
test_error_ols # 0.117
#-------------------------------------------------------------------------------
# RIDGE regression 
#-------------------------------------------------------------------------------

library(glmnet)

# 1) Finding the optimal lambda
set.seed(123)
ridge.cv=cv.glmnet(x = X.train, y = y.train, family = "gaussian", type.measure = "mse", alpha = 0)
# shrink param towards 0

plot(ridge.cv)
bestlam=ridge.cv$lambda.min
bestlam # 0.259

# 2) Fit the model 
ridge.fit=glmnet(x = X.train, y = y.train, family = "gaussian",type.measure = "mse", alpha = 0)
coef(ridge.fit,s=bestlam)  # s is the best lambda we found.
# predict(ridge.fit,type = "coefficients",s=bestlam)# same output of line above
n.var.ridge=length(which(coef(ridge.fit,s=bestlam)[-1]!=0)) # number of variables (excluding the intercept) different from 0 (for summary table)
n.var.ridge # still 168

# plot the coefficient path
annotations=colnames(X)
annotations[-which(abs(coef(ridge.fit,s=bestlam)[-1])>0.1)]=""
plot(ridge.fit, xvar="lambda",xlim=c(-5,5))
# if log lambda is 0 then lamda is 1
# if log lambda is - végtelen then lambda is 0
abline(v=log(bestlam), lty=2,col="black")
text(-5.5, coef(ridge.fit)[-1,length(ridge.cv$lambda)], annotations, pos=4,cex=0.6)

# x: lambda log scale 
# lines represent explanatory variables. what happens to coefs with dif values of lambda?
# vertical line: intersect of the coefs in the min squared error
# the values of coefs estimated with ridge regression
# the lowest line is the beta catine -- but thats also related to the response variable
# (dependency score)
# if you have muts on this gene, you have a higher dependency on the gene in question but a more - dependency score
# if mutations on e.g. mtor -> you will have a more + dependency score
# re-check dependency scores..
# similarly to simple lr, sign shows direction of relationship etc

# 3) Predictions
ridge.pred=predict(ridge.fit, s=bestlam, newx = X.test)
test_error_ridge=mean((ridge.pred-y.test)^2)
test_error_ridge # 0.06225617
#-------------------------------------------------------------------------------
# LASSO regression 
#-------------------------------------------------------------------------------

set.seed(123)
# 1) Finding the optimal lambda
lasso.cv=cv.glmnet(X.train,y.train,family = "gaussian", type.measure = "mse",alpha = 1)
# the alpha 1 is what defines that we do lasso
plot(lasso.cv)
bestlam.min=lasso.cv$lambda.min
bestlam.1se=lasso.cv$lambda.1se

# the two lines - cant say that one is statistically better than the other
# bc both lie within stdev
# 2 models (represented by the 2 lines)
# explain equally well our variable of interest 
# but the line on the right does it with fewer variables (?)
# also interesting how many variables the model has selected

# 2) Fit the model 
lasso.fit=glmnet(X.train,y.train,alpha = 1)
coef(lasso.fit,s=bestlam.1se)[which(coef(lasso.fit,s=bestlam.1se)!=0),,drop=FALSE]
coef(lasso.fit,s=bestlam.min)[which(coef(lasso.fit,s=bestlam.min)!=0),,drop=FALSE]
n.var.lasso=length(which(coef(lasso.fit,s=bestlam.min)[-1]!=0)) # number of variables (excluding the intercept) different from 0 (for summary table)
# make sure understand output - dependency score, etc, negative coef - negative effect on dependency score??
# which means it will reduce dependency score if there is a mutation on that gene.
# and a lower dependency score means the cell lines suffer more if the gene of interest (first col) is knocked out,
n.var.lasso # 5

# plot the coefficient path
labels=colnames(X)
labels=ifelse(coef(lasso.fit,s=bestlam.min)[-1]==0,"",labels)
plot(lasso.fit, xvar="lambda",xlim=c(-10.3,-2))
abline(v=c(log(bestlam.min),log(bestlam.1se)), lty=2,col="grey")
text(rep(-10.5, length(labels)), coef(lasso.fit)[-1,length(lasso.fit$lambda)], labels, pos=4,cex=0.6)
# can see that at the end you have all 0 for coefs except for a few coefs.
# we have 4 vars dif from 0 while all others = 0



# 3) Predictions
lasso.pred=predict(lasso.fit, s=bestlam.min, newx = X.test)
test_error_lasso=mean((lasso.pred-y.test)^2)
test_error_lasso # 0.06318595
#-------------------------------------------------------------------------------
# ELASTIC NET regression 
#-------------------------------------------------------------------------------

# 1) Finding the optimal alpha and lambda
alpha.grid=seq(0.1,0.9,0.1)
list.of.fits=list()
n.a=1
for (a in alpha.grid){
  set.seed(123)
  list.of.fits[[n.a]]=cv.glmnet(X.train, y.train, family = "gaussian", type.measure = "mse", alpha=a)
  names(list.of.fits)[n.a]=paste("alpha",a)
  n.a=n.a+1
}

min.cvm=c()
for (i in 1:length(list.of.fits)){
  min.cvm[i]=min(list.of.fits[[i]]$cvm)
}

bestalpha=alpha.grid[which.min(min.cvm)]
bestlam.min=list.of.fits[[which.min(min.cvm)]]$lambda.min
bestlam.1se=list.of.fits[[which.min(min.cvm)]]$lambda.1se

# 2) Fit the model 
en.fit=glmnet(X.train,y.train,family = "gaussian", type.measure = "mse",alpha=bestalpha)

coef(en.fit,s=bestlam.min)[which(coef(en.fit,s=bestlam.min)!=0),]
n.var.en=length(which(coef(en.fit,s=bestlam.min)[-1]!=0)) # number of variables (excluding the intercept) different from 0 (for summary table)
n.var.en # 5 like lasso
# plot the coefficient path
labels=colnames(X)
labels=ifelse(coef(en.fit,s=bestlam.min)[-1]==0,"",labels)
plot(en.fit, xvar="lambda",xlim=c(-10.5,-2))
abline(v=c(log(bestlam.min),log(bestlam.1se)), lty=2,col="grey")
text(rep(-10.8, length(labels)), coef(en.fit)[-1,length(en.fit$lambda)], labels, pos=4,cex=0.6)

# 3) Predictions
en.pred=predict(en.fit, s=bestlam.min, newx = X.test)
test_error_en=mean((en.pred-y.test)^2)
test_error_en # 0.06332069
#-------------------------------------------------------------------------------
### Summary table - models comparison
#-------------------------------------------------------------------------------

mse=c(test_error_ols,test_error_ridge,test_error_lasso,test_error_en)
n.var=c(n.var.lr,n.var.ridge,n.var.lasso,n.var.en)

table=cbind(mse,n.var)
rownames(table)=c("OLS","ridge","lasso","elastic net")

table
# elastic net basically gives same solution as lasso
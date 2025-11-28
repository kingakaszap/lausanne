#-------------------------------------------------------------------------------
# EXAMPLE: SIMULATED DATA  --- Lasso vs Elastic net --- grouped selection
#-------------------------------------------------------------------------------

n=100
set.seed(1)
Z1=runif(n,0,20)
Z2=runif(n,0,20)
Z3=runif(n,0,20)
y=Z1+0.2*Z2+rnorm(n)

X1=Z1+rnorm(n,sd=0.01)
X2=Z1+rnorm(n,sd=0.01)
X3=-Z1+rnorm(n,sd=0.01)
X4=Z2+rnorm(n,sd=0.01)
X5=-Z2+rnorm(n,sd=0.01)
X6=Z3

X=cbind(X1,X2,X3,X4,X5,X6)
temp=cor(X)

#-------------------------------------------------------------------------------
### Lasso
#-------------------------------------------------------------------------------

library(glmnet)

# 1) Finding the optimal lambda
set.seed(123)
lasso.cv=cv.glmnet(X,y,family = "gaussian", type.measure = "mse",alpha = 1)
bestlam=lasso.cv$lambda.min


# 2) Fit the model 
lasso.fit=glmnet(X,y,alpha = 1)
coef(lasso.fit,s=bestlam)
labels=colnames(X)
plot(lasso.fit, xvar="lambda",xlim=c(-4,2))
abline(v=log(bestlam), lty=2,col="grey")
text(rep(-4.3, length(labels)), coef(lasso.fit)[-1,length(lasso.cv$lambda)], labels, pos=4,cex=1)

#-------------------------------------------------------------------------------
## Elastic net 
#-------------------------------------------------------------------------------

# 1) Finding the optimal lambda
set.seed(123)
alpha=0.5 # arbitary alpha value
en.cv=cv.glmnet(X,y,family = "gaussian", type.measure = "mse",alpha = alpha)
bestlam=en.cv$lambda.min

# 2) Fit the model 
en.fit=glmnet(X,y,family = "gaussian", type.measure = "mse",alpha=alpha)
coef(en.fit,s=bestlam)

plot(en.fit, xvar="lambda",xlim=c(-3.5,2.5))
abline(v=log(bestlam), lty=2,col="grey")
text(rep(-3.7, length(labels)), coef(en.fit)[-1,length(en.cv$lambda)], labels, pos=4,cex=1)
# the ones that are symmetrical to y = 0 plane are correlated?? i think

#-------------------------------------------------------------------------------
#### Ridge
#-------------------------------------------------------------------------------

# 1) Finding the optimal lambda
set.seed(123)
ridge.cv=cv.glmnet(x = X, y = y, family = "gaussian", type.measure = "mse", alpha = 0)
bestlam=ridge.cv$lambda.min

# 2) Fit the model 
ridge.fit=glmnet(x = X, y = y, family = "gaussian", alpha = 0)
coef(ridge.fit,s=bestlam)

annotations=colnames(X)
plot(ridge.fit, xvar="lambda",xlim=c(-1.5,9))
abline(v=log(bestlam), lty=2,col="grey")
text(-1.5, coef(ridge.fit)[-1,length(ridge.cv$lambda)], annotations, pos=4,cex=1)

# ridge in this case more similar to elastic net;
# since not doing variable selection, it can better catch relationships (?)
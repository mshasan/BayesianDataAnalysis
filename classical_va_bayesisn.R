library(LearnBayes)
library(arm) # requires R 3.0 or greater.
library(MCMCpack)
#install.packages("C:\\Users\\Apu\\Downloads\\LaplacesDemon_13.11.17 (1).tar.gz\\")
library(LaplacesDemon)

## Probelm -a ----------------------------------

dat <- read.csv("school.csv",h=T);dat
attach(dat)
logY=log(Y)
summary(Y)
summary(logY)
var(logY)

## some initial plots
#-------------------------
plot(X,logY)
out = (logY > 3.4 & X < -5)
text(X[out], logY[out], label=School[out], pos = 2)
pairs(~logY+X,data=dat, main="Scatterplot Matrix")

## fit classical models(design matrix will be available)
#--------------------------------------------------------
fit=lm(logY~X,data=dat,x=TRUE,y=TRUE)
summary(fit)

modelb <- bayesglm(logY~X,data=dat,prior.mean=3.5,prior.scale=.033, prior.df=19)
display(modelb)



## The algorithm in binreg is based on the decomposition of the joint posterior
#------------------------------------------------------------------------------ 
theta.sample=blinreg(fit$y,fit$x,5000)
S=sum(fit$residual^2)
shape=fit$df.residual/2
rate=S/2
sigma2=rigamma(1,shape,rate)
MSE = sum(fit$residuals^2)/fit$df.residual
vbeta=vcov(fit)/MSE
beta=rmnorm(1,mean=fit$coef,varcov=vbeta*sigma2)

# histograms of the simulated posterior 
# draws of the individual regression coefficients
#--------------------------------------------------------------------------
par(mfrow=c(2,2))
hist(theta.sample$beta[,2],main="Socio economic status",xlab=expression(beta[1]))
hist(theta.sample$sigma,main="ERROR SD",xlab=expression(sigma))

# quantile command to summarize the draws of s
#------------------------------------------------
apply(theta.sample$beta,2,quantile,c(.05,.5,.95))
quantile(theta.sample$sigma,c(.05,.5,.95))

# histograms of simulated expected mean response
#------------------------------------------------------------------------
cov1=c(1,15)
X1=rbind(cov1)
mean.draws=blinregexpected(X1,theta.sample)
c.labels=c("A")
par(mfrow=c(2,2))
hist(mean.draws[,1],main=paste("Covariate set",c.labels[1]),xlab="logY")


# simulated predicted values of the parameters ß and s.
#-----------------------------------------------------------------------
cov1=c(1,15)
X1=rbind(cov1)
pred.draws=blinregpred(X1,theta.sample)
c.labels=c("A")
par(mfrow=c(2,2))
hist(pred.draws[,1],main=paste("Covariate set",c.labels[1]),xlab="logY")



# residuals check
#---------------------------------------------------------------------------
prob.out=bayesresiduals(fit,theta.sample,2)
par(mfrow=c(1,1))
plot(X,prob.out)
out = (prob.out > 0.35)
text(X[out], prob.out[out], label=School[out], pos = 4)


## problem - b-----------------------------
## draw scatter plot with regression line and 95% credible band

modelb <- bayesglm(logY~X,data=dat,prior.mean=3.5,prior.scale=.033, prior.df=19)
display(modelb)
par(mfrow=c(1,1))
fitted <- predict(modelb, interval = "confidence")

# plot the data and the fitted line
plot(X, logY)
lines(X, fitted)

# now the confidence bands
lines(X, lwr, lty = "dotted")
lines(X, upr, lty = "dotted")


## Problem - c-------------------------------------

p# simulated predicted values of the parameters ß and s.
#-----------------------------------------------------------------------
cov1=c(1,-16.04)
X1=rbind(cov1)
pred.draws=blinregpred(X1,theta.sample)
c.labels=c("New school")
par(mfrow=c(2,2))
hist(pred.draws[,1],main=paste("Covariate set",c.labels[1]),xlab="logY")

newdata=data.frame(X=-16.04)
predict(modelb,newdata,interval="prediction",level=0.95)




































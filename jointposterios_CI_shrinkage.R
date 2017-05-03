## problem - 1---------------------------------------------------
library(ggplot2)

#(c)---------------

y <- c(16/74,9/99,10/58,13/70,19/122,20/77,18/104,17/129,35/308,55/119);y
z <- c(12/125,1/19,2/16,4/48,9/217,7/74,9/38,8/162);z

lambda.y1 <- -log(prod(y));lambda.y1
lambda.y2 <- -log(prod(1-y));lambda.y2

alpha.y <- rexp(1000,lambda.y1);alpha.y
beta.y <- rexp(1000,lambda.y2);beta.y

lambda.z1 <- -log(prod(z));lambda.z1
lambda.z2 <- -log(prod(1-z));lambda.z2

alpha.z <- rexp(1000,lambda.z1);alpha.z
beta.z <- rexp(1000,lambda.z2);beta.z

posterior.y <- rbeta(1000,mean(alpha.y),beta.y);posterior.y
posterior.z <- rbeta(1000,mean(alpha.z),beta.z);posterior.z

results <- round(cbind(alpha.y,beta.y,posterior.y,alpha.z,beta.z,posterior.z),5);results




##(d)----------------------

mu.y <- alpha.y/(alpha.y+beta.y);mu.y
mu.z <- alpha.z/(alpha.z+beta.z);mu.z

mu.diff <- mu.y - mu.z; mu.diff
quantile(mu.diff, c(0.05, 0.95))

qplot(mu.diff, main="Mean difference", geom="histogram")



## problem - 2---------------------------------------------------

##(b)draw simulations from the joint posterior distributions----------

n <- c(74,99,58,70,122,77,104,129,308,119)
w <- c(16,9,10,13,19,20,18,17,35,55)
theta2 <- matrix(0,1000,10)
for (i in 1:10){
	theta2[,i] = rbeta(1000, w[i] + alpha, n[i]+beta-w[i])
	}


##(c) comparing raw and simulated proportins-------------

post.prop <- apply(theta2,2,mean)
obs.prop <- w/n
cbind(obs.prop,post.prop)


## (d) 95% posterior interval for average underlying proportions-------

m <- apply(theta2,2,median);m
y <- w/n
plot(y, m)
for (i in 1:10) {
	lami = rbeta(1000, w[i] + alpha, n[i]+beta-w[i])
	probint = quantile(lami, c(0.05, 0.95))
	lines(y[i] * c(1, 1), probint)
	}






## Problem - 3------------------------------------------------

## (b) 
## draw contours

library(LearnBayes)

n <- c(74,99,58,70,122,77,104,129,308,119)
w <- c(16,9,10,13,19,20,18,17,35,55)
dat <- data.frame(n,w);dat

## funtion need to create contour plots
poissgamexch=function (theta, datapar) # default theta=(2,-7)
{
	w = datapar$data[, 2]
	n = datapar$data[, 1]
	alpha = exp(theta[1])
	beta = exp(theta[2])
	
	logf = function(w, n, alpha, beta){
		lgamma(n + alpha) + n*log(beta) - lgamma(alpha) - (n + alpha)*log(1 + beta)
		}
	
	val = sum(logf(w, n, alpha, beta))
	return(val)
}


datapar = list(data = dat)
start=c(3, 4)
fit = laplace(poissgamexch, start, datapar)
fit

## contour plots
par(mfrow = c(1, 1))
mycontour(poissgamexch, c(0, 2.71, 0, 4.1), datapar, xlab="log alpha",ylab="log beta",
		main="Contour plots of hyperparameters")


## Scatter plot with contours from simulation (gibbs sampling)--------------

start = c(3, 4)
fitgibbs = gibbs(poissgamexch, start, 1000, c(1,.15), datapar)
fitgibbs$accept

par(mfrow = c(1, 2))
mycontour(poissgamexch, c(0, 2.71, 0, 4.1), datapar, xlab="log alpha",ylab="log beta")
points(fitgibbs$par[, 1], fitgibbs$par[, 2])
plot(fitgibbs$par[, 1], fitgibbs$par[, 2],xlab="log alpha",ylab="log beta" )


##(c)---------------
# explain from (b)

## (d) integrable density----------------------------
# by trial and error method we choose 

par(mfrow = c(1, 1))
mycontour(poissgamexch, c(-.7, 3.5, 1.5, 6), datapar, xlab="log alpha",ylab="log beta",
		main="Contour plots of hyperparameters")


start = c(3, 4)
fitgibbs = gibbs(poissgamexch, start, 1000, c(1,.15), datapar)
fitgibbs$accept

par(mfrow = c(1, 2))
mycontour(poissgamexch, c(-.7, 3.5, 1.5, 6), datapar, xlab="log alpha",ylab="log beta")
points(fitgibbs$par[, 1], fitgibbs$par[, 2])
plot(fitgibbs$par[, 1], fitgibbs$par[, 2],xlab="log alpha",ylab="log beta" )



## (e) simulation from joint posterior------------------------

n <- c(74,99,58,70,122,77,104,129,308,119)
w <- c(16,9,10,13,19,20,18,17,35,55)


alpha = exp(fitgibbs$par[, 1])
beta = exp(fitgibbs$par[, 2])
theta <- matrix(0,1000,10)
for (i in 1:10){
	theta[,i] = rgamma(1000, n[i] + alpha, beta/(1+beta))
	}


## (f)summaries findings------------------------


## comparing vechicle rate betwenn observed and simulated data
theta.obs <- n
theta.simu <- apply(theta,2,mean);theta.simu
cbind(theta.obs, theta.simu)


## shrinkage

shrink <- function(k) mean(alpha/(alpha + n[k]))
shrinkage=sapply(1:10, shrink)
plot(log(n), shrinkage)

























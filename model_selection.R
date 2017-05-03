# problem 2 r code
#####
# MAKE SURE THE DIRECTORY CONTAINS THE DATASET BANK

bank = read.table("banksalary.txt",head=T)
library(arm) # requires R 3.0 or greater.
library(MCMCpack)

# part a - scatterplot matrix--------------------

pairs(~BegSal+Educ+Exper+Time,data=bank, main="Scatterplot Matrix")

# part b - bayes with noninformative prior----------------

modelb <- bayesglm(BegSal~factor(Sex)+Educ+Exper+Time,data=bank,
	prior.scale=Inf, prior.df=Inf)
display(modelb)

# prediction-----------
newdata=data.frame(X=c(2,16,52.5,30))
predict(modelb,newdata,interval="prediction",level=0.95)


# quantile calculating--------------------
tempquant = function(x)(
	return(quantile(x,probs=c(.025,.975)))
)
apply(coef(sim(modelb)),2,tempquant)

# part c - frequentist analysis---------------------

freq.model = lm(BegSal~factor(Sex)+Educ+Exper+Time,data=bank)
summary(freq.model)

# part d - bayes with informative prior----------------------

modeld <- bayesglm(BegSal~factor(Sex)+Educ+Exper+Time,data=bank,
	prior.mean = c(0,47,0,0),
	prior.scale=Inf, prior.df=Inf)
display(modeld)

apply(coef(sim(modeld)),2,tempquant)

# part e - independence prior-------------------------

modele = MCMCregress(BegSal~factor(Sex)+Educ+Exper+Time,data=bank,
	b0=c(0,0,47,0,0))
summary(modele)

# part g - model selection------------------------------

modelf = MCMCregress(BegSal~factor(Sex)+Educ+Exper+Time,data=bank,
	marginal.likelihood="Chib95",
	b0=c(0,0,47,0,0),B0=c(.001,.001,.001,.001,.001),c0=1,d0=1)

modelf1 = MCMCregress(BegSal~Educ+Exper+Time,data=bank,
	marginal.likelihood="Chib95",
	b0=c(0,47,0,0),B0=c(.001,.001,.001,.001),c0=1,d0=1)
modelf2 = MCMCregress(BegSal~factor(Sex)+Exper+Time,data=bank,
	marginal.likelihood="Chib95",
	b0=c(0,0,0,0),B0=c(.001,.001,.001,.001),c0=1,d0=1)
modelf3 = MCMCregress(BegSal~factor(Sex)+Educ+Time,data=bank,
	marginal.likelihood="Chib95",
	b0=c(0,0,47,0),B0=c(.001,.001,.001,.001),c0=1,d0=1)
modelf4 = MCMCregress(BegSal~factor(Sex)+Educ+Exper,data=bank,
	marginal.likelihood="Chib95",
	b0=c(0,0,47,0),B0=c(.001,.001,.001,.001),c0=1,d0=1)

test = BayesFactor(modelf,modelf1,modelf2,modelf3,modelf4)

#Problem 3
#Part a)

#open the library "coda"
library(coda)


#introduce the data
data=c(65,156,100,134,16,108,121,4,39,143,56,26,1,1,5,65)


#alpha has prior Gamma(a,b), lambda has prior Gamma(c,d)
#Choose a,b,c,d and niter, number of iterations

a=1;
b=3;
c=5;
d=20;
n<-length(data);
niter=40000;

#Calculate the log of joint posterior density
 
logpostfct<-function(alp,lam)
{  return((n+a-1)*log(alp)+(n+c-1)*log(lam)-alp/b-lam*(1/d+sum(data^alp))
           +(alp-1)*sum(log(data)));
}


# This function builds the MCMC chain; it outputs two vectors 
#One vector contains the alphas one the lambdas
#start is the vector of starting values, sig1 and sig2 are the scale
#factors that update alpha and lambda, respectively
#nsim is the number of simulations(length of the chain)



metrhast<-function(start, sig1, sig2, nsim)
{
alphav=rep(start[1], nsim);  #initialize the vectors containing the
lambdav=rep(start[2], nsim);  #alpha and lambda values

for(i in 2:nsim)
{
  alphav[i]<-alphav[i-1]+sig1*rnorm(1);  #update the parameters
  lambdav[i]<-lambdav[i-1]+sig2*rnorm(1); #make sure they both are positive
 while((alphav[i]<=0) | (lambdav[i]<=0))
  { alphav[i]<-alphav[i-1]+sig1*rnorm(1);
    lambdav[i]<-lambdav[i-1]+sig2*rnorm(1);
  }
  
  loglik_ratio=logpostfct(alphav[i], lambdav[i])-
               logpostfct(alphav[i-1], lambdav[i-1]);
  u<-runif(1);
  if(u>min(1, exp(loglik_ratio)))
   {alphav[i]<-alphav[i-1];
    lambdav[i]<-lambdav[i-1];
   }  
}
  newlist<-list("alphavalues"=alphav, "lambdavalues"=lambdav);
  return(newlist);
}




#We  will use these values for part e)
sigma1<-1.0;
sigma2<-0.5;

start1=c(1.0, 1.0);
start2=c(5.0, 3.0);
start3=c(0.5, 0.1);
 

#Part b)


#Find loglikelihood function of posterior of alpha given lambda and data

log.lik<-function(alp,lam)
{
return((n+a-1)*log(alp)+(alp-1)*sum(log(data))-alp/b-lam*sum(data^alp));
}

#This function builds the chains
# startv is the vector of starting values of the chain
#sig is the scale factor used to update alpha
#nsimul is the number of iterations
#datav is the empirical data


Gibbsampling<-function(startv, sig, nsimul, datav)
{
alphavect=rep(startv[1], nsimul);  #initialize the vectors containing the
lambdavect=rep(startv[2], nsimul);  #alpha and lambda values



#Do MH inside Gibbs sampling

for(i in 2:niter)
{
  alphavect[i]<-alphavect[i-1]+sig*rnorm(1);
 while(alphavect[i]<=0)
  { alphavect[i]<-alphavect[i-1]+sig*rnorm(1);}
  
  loglik_ratio=log.lik(alphavect[i], lambdavect[i-1])-
               log.lik(alphavect[i-1], lambdavect[i-1]);
  u<-runif(1);
  if(u>min(1, exp(loglik_ratio)))
   {alphavect[i]<-alphavect[i-1];}

   lambdavect[i]<-rgamma(1, shape=n+c, rate=1/d+sum(data^alphavect[i]));

}
  newliste<-list("alphas"=alphavect, "lambdas"=lambdavect);
  return(newliste);
}

#We will use this sigma value in part e)
sigma=0.5;






#Part e)

#Run the algorithms from part a) and b) 3 times
#Plot the trace of plots and autocorrelation functions

sigma1<-1.0;
sigma2<-0.5;
sigma<-0.5;

start1=c(1.0, 1.0);
start2=c(5.0, 3.0);
start3=c(0.5, 0.1);

res1<-metrhast(start1, sigma1, sigma2, niter);
matr1<-cbind(res1$alphavalues, res1$lambdavalues);
colnames(matr1)<-c("alpha", "lambda");
m1<-as.mcmc(matr1);
plot(m1, trace=TRUE, density=TRUE);
autocorr.plot(m1, auto.layout=FALSE)
xyplot(m1, col="blue");

res2<-metrhast(start2, sigma1, sigma2, niter);
matr2<-cbind(res2$alphavalues, res2$lambdavalues);
colnames(matr2)<-c("alpha", "lambda");
m2<-as.mcmc(matr2);
xyplot(m2, col="red");
autocorr.plot(m2, auto.layout=FALSE);
thinr2<-window(m2, 15001, 40000, 25);
autocorr.plot(thinr2, auto.layout=FALSE);

res3<-metrhast(start3, sigma1, sigma2, niter);
matr3<-cbind(res3$alphavalues, res3$lambdavalues);
colnames(matr3)<-c("alpha", "lambda");
m1<-as.mcmc(matr3);
xyplot(m2, col="green");
autocorr.plot(m3, auto.layout=FALSE);


gres1<-Gibbsampling(start1, 0.5, 40000, data)
gmatr1<-cbind(gres1$alphas, gres1$lambdas);
colnames(gmatr1)<-c("alpha", "lambda");
gm1<-as.mcmc(gmatr1);
autocorr.plot(gm1, auto.layout=FALSE)
xyplot(gm1, col="green");
thing1<-window(gm1, 15001, 40000, 25);
autocorr.plot(thing1, auto.layout=FALSE)


gres2<-Gibbsampling(start2, 0.5, 40000, data)
gmatr2<-cbind(gres2$alphas, gres2$lambdas);
colnames(gmatr2)<-c("alpha", "lambda");
gm2<-as.mcmc(gmatr2);
autocorr.plot(gm2, auto.layout=FALSE)
xyplot(gm2, col="green");
thing2<-window(gm2, 15001, 40000, 25);
autocorr.plot(thing2, auto.layout=FALSE)


gres3<-Gibbsampling(start3, 0.5, 40000, data)
gmatr3<-cbind(gres3$alphas, gres3$lambdas);
colnames(gmatr3)<-c("alpha", "lambda");
gm3<-as.mcmc(gmatr3);
autocorr.plot(gm3, auto.layout=FALSE)
xyplot(gm3, col="green");
thing3<-window(gm3, 15001, 40000, 25);
autocorr.plot(thing3, auto.layout=FALSE)


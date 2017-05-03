dat1<-read.delim("nejm_brca_release.txt",as.is=T)
dat2<-dat1[,4:25]  #get rid of gene info columns
dat2<-as.matrix(dat2)  #convert from data frame to matrix

#make plots comparing ratios to log(ratios):-----------------------------
par(mfrow=c(3,2))
hist(dat2[,1],xlab="expression ratio",main="sample individual #1")
hist(log(dat2[,1]),xlab="log(expression ratio)",main="sample individual #1")

hist(dat2[,4],xlab="expression ratio",main="sample individual #4")
hist(log(dat2[,4]),xlab="log(expression ratio)",main="sample individual #4")

hist(dat2[,6],xlab="expression ratio",main="sample individual #6")
hist(log(dat2[,6]),xlab="log(expression ratio)",main="sample individual #6")


## frequentist clustering--------------------

dat2 <- log(dat2)
HedDist <- dist(dat2)
HedenfalkHclust <- hclust(HedDist)
par(mfrow=c(1,1))
plot(HedenfalkHclust)

## we can zoom in to smaller parts of the tree-------------

h10 <- cutree(HedenfalkHclust,h=10)
HedClusth10 <- hclust(dist(dat2[h10==1,]))
plot(HedClusth10,labels=dat1[h10==1,1])

## again zoom in------------------

k3 <- cutree(HedClusth10,k=3)
HedClustk3 <- hclust(dist(dat2[h10==1,][k3==3,]))
plot(HedClustk3,labels=dat1[h10==1,][k3==3,1])


## bayesian clustering----------------------------------

#install.packages("bayesclust",repos="http://mirrors.nics.utk.edu/cran/")
#install.packages("mclust",repos="http://mirrors.nics.utk.edu/cran/")
library(bayesclust)
library(mclust)

dim(dat2)# 3226   22
# load a workspace into the current session-----------------------
#load("HedFalkBayesClust")


# we are using 'dat2 <- log(dat2)' log of original data
# Search for optimal partitioning of data into 2 clusters---------------

clusterTest3 <- cluster.test(dat2, nsim=100000, p=22, k=3, mcs=0.1,file = "HedFalkBayesClust",replications=4)
summary(clusterTest3)


# Plot the running posterior probabilities to monitor convergence-------
plot(clusterTest)
#plot(cluster.test.reps)


# Generate corresponding null density object.-----------------------
#nulldensK3 <- nulldensity(n=3226, nsim=10000, k=3, mcs=0.1, p=22, prop=0.25)
#hist(nulldensK3, main="Null Density Histogram", xlab="samples")


# Convert EPP to p-value of frequentist p-value--------------------------------------------
#emp2pval(clusterTest, nulldensK3)


# Searching for Optimal Clusters----------------------------
#clusterOptK3 <- cluster.optimal(dat2, nsim=1000, p=22, k=3, mcs=0.1)
#plot(clusterOptK3)


# multiple test for clusterer size-----------------------------
clusterTestK2 <- cluster.test(dat2, nsim=100000, p=22, k=2, mcs=0.1,replications=4)
clusterTestK4 <- cluster.test(dat2, nsim=100000, p=22, k=4, mcs=0.1,replications=4)
clusterTestK5 <- cluster.test(dat2, nsim=100000, p=22, k=5, mcs=0.1,replications=4)
combine(clusterTestK2, clusterTest3, clusterTestK4,clusterTestK5)

#Using the Default Null Distributions----------------------------
data(cutoffs)
tail(cutoffs)

#Comparison with Mclust-------------------------------------------

X <- as.matrix(dat2)
#clusterTestK3 <- cluster.test(X, nsim=10000, p=22, k=3, mcs=0.1,replications=4)
#nulldensK3 <- nulldensity(n=3226, nsim=1000, k=3, mcs=0.1, p=22, prop=0.25)
#emp2pval(clusterTestK3, nulldensK3)
model = Mclust(X, prior = priorControl(functionName="defaultPrior", shrinkage=0.1))
summary(model)




# save the workspace to the file .RData in the cwd 
save.image("p2.RData")

# save specific objects to a file
# if you don't specify the path, the cwd is assumed 
#save(object list,file="myfile.RData")---------------------------------



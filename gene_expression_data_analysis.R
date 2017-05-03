## library needs for this analysis
#---------------------------------
library(bayesclust)
library(BayesFactor)


## data reading
#----------------
dat1<-read.delim("nejm_brca_release.txt",as.is=T)
dat2<-dat1[,4:25]  #get rid of gene info columns
dat2<-as.matrix(dat2)  #convert from data frame to matrix
dat3<-log(dat2)
head(dat2, n=3)

## ---------------------starts data summary-------------------------
# scatterplot matrix for samples 1, 5, 11,16 & 21
#---------------------------------------------------------
summary(dat2)
summary(dat3)
pairs(~s1996+s1252+s1572+s1063,data=dat2, main="Scatterplot Matrix") # given data
pairs(~s1996+s1252+s1572+s1063,data=dat3, main="Scatterplot Matrix") # log-transformed data

##make plots comparing ratios to log(ratios)
#------------------------------------------------------
par(mfrow=c(3,2))
hist(dat2[,1],xlab="expression level",main="sample individual #1")
hist(log(dat2[,1]),xlab="log(expression level)",main="sample individual #1")
hist(dat2[,11],xlab="expression level",main="sample individual #11")
hist(log(dat2[,11]),xlab="log(expression level)",main="sample individual #11")
hist(dat2[,21],xlab="expression ratio",main="sample individual #21")
hist(log(dat2[,21]),xlab="log(expression level)",main="sample individual #21")


## -----------------starts analysis-----------------------------
## for clustering creating a data log transformed data
#---------------------------------------------------------------
dat_log1<-log(dat2)
#write.table(dat_log1, file = "hedenfalk_logfull.txt",sep = "\t",row.names = F,col.names = F)
dat_log2 <- read.delim("nejm_brca_release_logfull.txt",as.is=T)
dat_log3 <- dat_log2[,3:24]
dat_log <- as.matrix(dat_log3)
head(dat_log, n=3)

############################################################################
# Z-cluster part (won't be able to use in personal computer)---------------
############################################################################
# multiple test for clusterer size
#---------------------------------
#clusterTest3 <- cluster.test(dat_log, nsim=5000000, p=22, k=3, mcs=0.1,file = "HedFalkBayesClust",replications=4)
#summary(clusterTest3)
#clusterTestK2 <- cluster.test(dat_log, nsim=5000000, p=22, k=2, mcs=0.1,replications=4)
#clusterTestK4 <- cluster.test(dat_log, nsim=5000000, p=22, k=4, mcs=0.1,replications=4)
#clusterTestK5 <- cluster.test(dat_log, nsim=5000000, p=22, k=5, mcs=0.1,replications=4)
#combine(clusterTestK2, clusterTest3, clusterTestK4,clusterTestK5)
#clusterOptK3 <- cluster.optimal(dat2, nsim=100000, p=22, k=3, mcs=0.1)
#plot(clusterOptK3)
#################################################################################################################

# Clustring---------------------------------------
# frequentists sample hierarchical clustering
#-------------------------------------------------
HedDist <- dist(t(dat_log))
freq.clust <- hclust(HedDist,method = "ward")
par(mfrow=c(1,1))
plot(freq.clust,,xlab="Samples")


# frequentists gene hierarchical clustering
#-------------------------------------------------
HedDist <- dist((dat_log))
freq.clust <- hclust(HedDist,method = "complete")
par(mfrow=c(1,1))
plot(freq.clust,,xlab="Samples")

h10 <- cutree(freq.clust,h=10)
HedClusth10 <- hclust(dist(dat_log[h10==1,]))
plot(HedClusth10,labels=dat1[h10==1,1])

k3 <- cutree(HedClusth10,k=3)
Hedclustk3 <- hclust(dist(dat_log[h10==1,][k3==3,]))
plot(Hedclustk3,labels=dat1[h10==1,][k3==3,1])


## bayesian ANOVA testing by bayesfactor.(BF = alt model/null model)
BF>1 to inf increases significance of the alt. model,10=strong evidence
#----------------------------------------------------------------------
#group1(s1996, s1649, s1320, s1324, s1905, s1542)
#group2(s1822, s1252, s1900, s1486, s1224, s1510, s1321, s1281, s1714, s1787, s1721, s1572) 
#group3(s1816, s1616, s1063, s1936)

group1 <- c(1,12:15,18)
group2 <- c(2:11,16,17)
group3 <- c(19:22)

bayesFactor <- NULL
for(x in 1:nrow(dat_log))
	{
	Y <- c(dat_log[x,group1],dat_log[x,group2],dat_log[x,group3])
	X1 <- c(1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
	X2 <- c(0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0)
	dat <- data.frame(Y,X1,X2)
	bf = lmBF(Y ~ X1*X2,data=dat,progress=T)
	bayesFactor <- c(bayesFactor,extractBF(bf)[1,1])
	}
(BayesSignGene <- sum(bayesFactor > 10)) # 398 significant genes.
signGenebyBF <- dat1[,1][order(bayesFactor)][1:BayesSignGene] # order gives position of the values of sort()



## frequentists ANOVA testing
#--------------------------------
freq.pvals = c()
for(x in 1:nrow(dat_log))
	{
	Y <- c(dat_log[x,group1],dat_log[x,group2],dat_log[x,group3])
	X1 <- c(1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
	X2 <- c(0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0)
	freq.anova <- lm(Y~X1*X2)
	a <- summary(freq.anova)
	pval <- 1-pf(a$fstatistic[1],a$fstatistic[2],a$fstatistic[3])
	freq.pvals = c(freq.pvals,pval)
	}

length(which(freq.pvals < .05)) # gives 1000 significant genes
length(which(freq.pvals < .01)) # gives 516 significant genes


#Benjemini and Hochberg FDR to get # of significant genes
# it is valid only for frequentist analysis
#---------------------------------------------------------
# FDR gives 438 significant genes 
FDR=0.05
g<-nrow(dat_log) 		#number of genes
rval<-(1:g)/g*FDR
spval<-sort(freq.pvals)
(numsig<-sum(cumsum(spval<rval)==(1:g))) #This gives the number of significant genes at FDR=0.05
sigGenebyFDR <- dat1[,1][order(freq.pvals)][1:numsig] #this gives the significant genes. 

#signGenebyBF[sapply(signGenebyBF ,function(x) any(grepl(x,sigGenebyFDR)))]



## getting significant genes and corresponding rows and columns
# from bayesina ANOVA (we created a data set of significant genes)
#---------------------------------------------------------------------
datSigGenes <- read.delim("signGenesBaselineG1.txt",as.is=T)
attach(data.frame(datSigGenes))

x <- as.vector(datSigGenes)				# converting "SZ.RLV" column as vector
y <- strsplit(x, '[ ]')				# seperating numerical and character values
y2 <- sapply(y, '[[', 1)			# keeping only numerical values for further calculations
S.size <- as.numeric(y2)




Gene Label
datSigGenes <- as.matrix(datSigGenes)[,1]



significant genes at FDR=0.05
sigGenebyFDR <- dat1[,1][order(freq.pvals)][1:numsig] #this gives the significant genes. 

signGenebyBF[sapply(signGenebyBF ,function(x) any(grepl(x,sigGenebyFDR)))]

dat1[,1] %in% sigGenebyFDR == 1 # t/f of the sig genes
which(dat1[,1] %in% sigGenebyFDR == 1) # which rows are the significant genes
dat1[which(dat1[,1] %in% sigGenebyFDR == 1),] # the dataset with only selects the rows of sig genes








# two sample t-test among 3 group
#-----------------------------------------------------
grp1 <- c(19:22)
grp2 <- c(1,12:15,18)
grp3 <- c(2:11,16,17)

pval_12 <- NULL
pval_13 <- NULL
pval_23 <- NULL
for(x in 1:nrow(dat2))
	{
	pval_12<-c(pval_12,t.test(dat_log[x,grp1],dat_log[x,grp2])$p.value)
	pval_13<-c(pval_13,t.test(dat_log[x,grp1],dat_log[x,grp3])$p.value)
	pval_23<-c(pval_23,t.test(dat_log[x,grp2],dat_log[x,grp3])$p.value)
	}

#Benjemini and Hochberg FDR
#----------------------------
FDR = 0.05
#FDR = 0.20
g<-nrow(dat_log) #number of genes
rval_12<-(1:g)/g*FDR
rval_13<-(1:g)/g*FDR
rval_23<-(1:g)/g*FDR

spval_12<-sort(pval_12)
spval_13<-sort(pval_13)
spval_23<-sort(pval_23)

(numsig1<-sum(cumsum(spval_12<rval_12)==(1:g))) #gives the number of significant genes at FDR=0.05
(numsig2<-sum(cumsum(spval_13<rval_13)==(1:g)))
(numsig3<-sum(cumsum(spval_23<rval_23)==(1:g)))

dat1[,1][order(pval_12)][1:numsig1] 	#this gives the significant genes. 
dat1[,1][order(pval_13)][1:numsig2]
dat1[,1][order(pval_23)][1:numsig3]











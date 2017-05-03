######################################################
##### problem 2

##########
# inputting the known data

# number of bikes
w = c(16,9,10,13,19,20,18,17,35,55)

# number of non bikes
oth = c(58,90,48,57,103,57,86,112,273,64)

# total number of vehicles
n = w+oth

# proportion of bikes
y = w/n

##########
# part b.)

#####
# making the hyperprior grid of alpha and betas

alpha_grid = seq(0.01,10,0.01)
beta_grid = seq(0.02,20,0.02)

#####
# evaluating the hyperprior density at each (alpha,beta) gridpoint

# we'll write a quick function that will evaluate a
# given (alpha,beta) pair on the hyperprior distribution

hyperprior = function(alpha,beta){

	holder = c(lgamma(alpha+beta)-lgamma(alpha)-lgamma(beta))

	# iterating over the 10 sampled streets
	for(i in 1:10){
	
		#evaluating the density factor at street i
		density_i = lgamma(alpha+w[i])+lgamma(n[i]-w[i]+beta)-lgamma(alpha+beta+n[i])	

		# adding that value to a holding vector
		holder = c(holder,density_i)

	}

	# multiplying the density factors and outputting their product
	output = exp(sum(holder))
	return(output)

}

# using the above function, we'll run through all possible pairs of
# alpha and beta grid values. note that we'll have to get a little cute
# and have to do the exp(ln(x)) trick for computational reasons

ab_grid = matrix(0,ncol=length(alpha_grid),nrow=length(beta_grid))

# then by row, we'll fill in the matrix ab_grid with the densities
for(i in 1:length(beta_grid)){
	for(j in 1:length(alpha_grid)){
		ab_grid[i,j] = hyperprior(alpha=alpha_grid[j],beta=beta_grid[i])
	}
}

# and to create a proper distribution, we'll normalize the densities
# by their sum

ab_grid = ab_grid/sum(ab_grid)

#####
# marginal distributions of the hyperparameters

# the marginal distribution of the alphas are the column sums
# the marginal distribution of the betas are the row sums

par(mfrow=c(2,1))
plot(alpha_grid,colSums(ab_grid),type="l",xlim=c(0,4),
	main="Marginal distribution of hyperparameter alpha")
plot(beta_grid,rowSums(ab_grid),type="l",xlim=c(0,4),
	main="Marginal distribution of hyperparameter beta")
par(mfrow=c(1,1))

#####
# random samples from the joint hyperprior

# we want to get 1000 random samples from the joint hyperprior
# we'll pick alpha from its marginal distribution, and then conditional
# on that alpha value, sample a beta. lastly, we'll add a little jitter
# to each alpha/beta pair

a_sample = c()
b_sample = c()

for(i in 1:1000){

	# sampling an alpha and beta value
	a_temp = sample(alpha_grid,1,prob=colSums(ab_grid))
	b_temp = sample(beta_grid,1,prob=ab_grid[,which(a_temp==alpha_grid)])
	
	# adding some noise
	a_sample = c(a_sample,a_temp+runif(1,min=-0.005,max=0.005))
	b_sample = c(b_sample,b_temp+runif(1,min=-0.01,max=0.01))

}

########## 
# part c.)

# now we want to use the drawn hyper parameter sample in part b to draw
# a sample of the theta_1 to theta_10

# we'll write a quick function, since we need to do this 10 times,
# once for every theta

thetasamp = function(which,name,posteriorprob){
	
	temp = c()
	
	# for every (alpha, beta) pair, we draw and make a theta_i sample
	for(i in 1:length(a_sample)){
		
		temp = c(temp,rbeta(1,a_sample[i]+w[which],b_sample[i]+n[which]-w[which]))
		
	}
	
	# creating a (global) variable named "theta[i]" to hold all the theta[i] samples
	assign(paste("theta",name,sep=""),temp,envir=.GlobalEnv)
	
	# plotting the histogram of the theta
	breaks = seq(0,0.8,0.025)
	hist(temp,breaks=breaks,main=paste("Histogram of theta",name,sep=""))
	abline(v=y[which],lty=2,lwd=2)
	
	# returning a posterior probability interval for the sampled theta
	low = (1-posteriorprob)/2
	high = 1-low
	ci = c(quantile(temp,low),quantile(temp,high))
	abline(v=ci[1],lwd=2)
	abline(v=ci[2],lwd=2)	
	return(ci)
	
}

par(mfrow=c(5,2))
thetasamp(1,"1",.95)
thetasamp(2,"2",.95)
thetasamp(3,"3",.95)
thetasamp(4,"4",.95)
thetasamp(5,"5",.95)
thetasamp(6,"6",.95)
thetasamp(7,"7",.95)
thetasamp(8,"8",.95)
thetasamp(9,"9",.95)
thetasamp(10,"10",.95)
par(mfrow=c(1,1))

########## 
# part d.)

# note that the posterior interval for the average underlying prop of
# bicycle traffic is alpha/(alpha+beta). again we use the hyper parameter
# sample we drew in part b to get our estimates of the underlying theta

# calculating a sample of the underlying theta
underlying = a_sample/(a_sample+b_sample)
hist(underlying,main="Histogram of the average underlying proportion of bike traffic")

low = (1-0.95)/2
high = 1-low
ci = c(quantile(underlying,low),quantile(underlying,high))
abline(v=ci[1],lwd=2)
abline(v=ci[2],lwd=2)	

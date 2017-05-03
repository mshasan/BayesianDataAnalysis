library(foreign)
example.1 <- read.spss("http://www.unt.edu/rss/class/Jon/R_SC/Module3/ExampleData1.sav",
   use.value.labels=TRUE, max.value.labels=Inf, to.data.frame=TRUE)
summary(example.1)

library(Rcmdr)
library(abind)

# NOTE: The package 'BayesFactorPCL' is not available (currently; Mar. 3, 2011) at CRAN and must 
# be installed from R-Forge; so, use the following line of code in the R console to install it 
# directly from R-Forge:

install.packages("BayesFactorPCL", repos="http://R-Forge.R-project.org")

# Then load the library into the current session. 

library(BayesFactorPCL)

# For a description of the package:

help(BayesFactorPCL)

#### T-test example.

# Summary for types of Candy.

numSummary(example.1$Recall1 , groups=example.1$Candy,statistics=c("mean", "sd"))

boxplot(example.1$Recall1 ~ example.1$Candy, col = "lightgreen")

# Levene's test of Homogeneity of Variances ("var").

tapply(example.1$Recall1, example.1$Candy, var, na.rm=TRUE)
leveneTest(example.1$Recall1, example.1$Candy, center=median)

# First, conduct the traditional t-test.

t.t1 <- t.test(Recall1~Candy, alternative="less", conf.level=.95, var.equal=TRUE, data=example.1)
t.t1

attach(example.1)
x1 <- split(Recall1, Candy)
x1
detach(example.1)

library(MBESS)
smd(x1$Skittles, x1$None)

# Second, one may want to conduct a Bayesian version of the Levene's test for homogeneity of 
# variances; using the 'BayesFactorPCL' library function 'eqVariance.Gibbs'; which requires a 
# matrix of data with each group as a column and each row a case. In the output, the "$BF" is 
# the Bayes Factor, which if greater than one indicates the groups' variances are equal. 

x2 <- cbind(x1$Skittles, x1$None)
is.matrix(x2)

eqVariance.Gibbs(x2, iterations = 1000, whichModel = 2, M2.metrop.scale = 2)
 
# Third, conduct the Bayes Factor t-test; which returns "a scalar giving the Bayes factor IN FAVOR  
# of the NULL HYPOTHESIS that the effect size is 0.." (Rouder & Morey, ttest.Quad function help document).

t.b1 <- ttest.Quad(t = -7.7566, n1 = 27, n2 = 27, rscale = 1, prior.cauchy = TRUE)
t.b1

t.b2 <- ttest.Quad(t = -7.7566, n1 = 27, n2 = 27, rscale = 1, prior.cauchy = FALSE)
t.b2


help(ttest.Quad)

#### One Way ANOVA example.

# Summary for types of distraction.

numSummary(example.1$Recall1 , groups=example.1$Distraction,statistics=c("mean", "sd"))

# Levene's test of Homogeneity of Variances ("var").

tapply(example.1$Recall1, example.1$Distraction, var, na.rm=TRUE)
leveneTest(example.1$Recall1, example.1$Distraction, center=median)

boxplot(example.1$Recall1 ~ example.1$Distraction, col = "lightgreen")

# First conduct the traditional ANOVA (Distraction has 3 groups, each with 18 cases).

aov.t1 <- aov(Recall1 ~ Distraction, data=example.1)
summary(aov.t1)

# Second, conduct the Bayes Factor analysis (only works with balanced designs right now: Mar. 2011); 
# which returns "a scalar giving the Bayes factor in favor of the null hypothesis that the 
# standardized effect size is 0" (Morey, oneWayAOV function help documentation). 

aov.b1 <- oneWayAOV.Quad(F = 2.1164, N = 18, J = 3, rscale = 1)
aov.b1

help(oneWayAOV.Quad)

## Note on interpretation of Bayes Factors: Jeffreys (1961) recommends that odds greater than 3 be considered 
## some evidence, odds greater than 10 be considered strong evidence, odds greater than 30 be considered very 
## strong evidence for one hypothesis over another. In the one way ANOVA example above, the Bayes Factor value 
## is 3.234, which indicates that the NULL hypothesis is 3.234 times more probable than the ALTERNATIVE 
## hypothesis, given the data. Kass and Raftery (1995) offer a slightly different strategy for interpreting 
## Bayes Factors: 1 to 3.2 not worth mentioning, 3.2 to 10 substantial, 10 to 100 strong, and greater than 
## 100 decisive. 

# Jeffreys, H. (1961). Theory of probability (3rd ed.). Oxford: Oxford University Press.

# Kass, R. E., & Raftery, A. E. (1995). Bayes Factors. Journal of the American Statistical 
#     Association, 90, 773 - 795. 

# NOTE: the ‘LearnBayes’ package, which is a companion for the book Bayesian Computation with 
# R, both of which authored by Jim Albert (2010, 2007); also contains functions for computing 
# Bayes Factors. 

# Albert, J. (2007). Bayesian computation with R. New York: Springer Science+Business Media, LLC.

# Albert, J. (2010). Package ‘LearnBayes’. Available at CRAN: 
# http://cran.r-project.org/web/packages/LearnBayes/index.html


###### Links for some references/resources

# http://www.socsci.uci.edu/~mdlee/WetzelsEtAl2010.pdf

# http://cran.r-project.org/web/packages/mcmc/vignettes/bfst.pdf

# http://www.stat.cmu.edu/~kass/papers/bayesfactors.pdf

# https://r-forge.r-project.org/projects/bayesfactorpcl/




# End: Last updated, Mar. 10, 2011; (added Bayes Factor Levene's test).
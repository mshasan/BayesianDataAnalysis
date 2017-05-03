#install the package MCMCpack to do a bayesian factor analysis using MCMC methods

library(MCMCpack)


#Read in the whole data and transform it into a matrix
#Whole data
data1<-read.delim("hedenfalk_full.txt",as.is=T, header=FALSE)
data2<-log(data1)
data2<-as.matrix(data2, ncol=22)
reardat<-cbind(data2[,c(1,12:15,18)], data2[,c(2:11, 16, 17)], data2[,c(19:22)])
colnames(reardat)<-c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11","x12",
"x13", "x14","x15","x16","x17","x18","x19","x20","x21","x22")
whdata<-data.frame(reardat)


#Group1 corresponds tovars x1-x6, group2 corresponds to vars x7-x18 and group3
#corresponds to x19-x22



#Do an usual factor analysis on the whole data using 3 factors
whdfan.1<-factanal(~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20+x21+x22, factors=3, rotation="varimax", data=whdata)
whdfan.1

#Results
Loadings:
    Factor1 Factor2 Factor3
x1  0.476   0.262   0.500  
x2  0.405   0.461   0.611  
x3  0.363   0.319   0.685  
x4  0.342   0.311   0.802  
x5  0.414   0.323   0.420  
x6  0.427   0.237   0.560  
x7  0.739   0.257   0.347  
x8  0.572   0.304   0.320  
x9  0.744   0.313   0.282  
x10 0.691   0.302   0.401  
x11 0.706   0.340   0.264  
x12 0.677   0.314   0.433  
x13 0.552   0.385   0.400  
x14 0.457   0.361   0.481  
x15 0.684   0.370   0.319  
x16 0.556   0.440   0.209  
x17 0.687   0.186   0.307  
x18 0.673   0.224   0.268  
x19 0.330   0.730   0.415  
x20 0.277   0.738   0.220  
x21 0.346   0.789   0.224  
x22 0.296   0.738   0.451  

               Factor1 Factor2 Factor3
SS loadings      6.463   4.140   4.123
Proportion Var   0.294   0.188   0.187
Cumulative Var   0.294   0.482   0.669

Test of the hypothesis that 3 factors are sufficient.
The chi square statistic is 4956.65 on 168 degrees of freedom.
The p-value is 0 

#Group1 corresponds to factor3, group2 to factor1 and group3 to factor2


#Do a Bayesian factor analysis of the whole data using MCMC method
 
whdfan.2<-MCMCfactanal(~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20+x21+x22, factors=3, data=whdata, 
lambda.constraints=list(x1=c(1,0), x1=c(2,0), x2=c(1,0), x2=c(2,0), x3=c(1,0), x3=c(2,0),
x4=c(1,0), x4=c(2,0), x5=c(1,0), x5=c(2,0), x6=c(1,0), x6=c(2,0), x7=c(2,0), x7=c(3,0),
x8=c(2,0), x8=c(3,0), x9=c(2,0), x9=c(3,0), x10=c(2,0), x10=c(3,0), x11=c(2,0), x11=c(3,0), 
x12=c(2,0), x12=c(3,0), x13=c(2,0), x13=c(3,0),  x14=c(2,0), x14=c(3,0), x15=c(2,0), x15=c(3,0), x16=c(2,0), x16=c(3,0),
x17=c(2,0), x17=c(3,0), x18=c(2,0), x18=c(3,0), x19=c(1,0), x19=c(3,0), x20=c(1,0), x20=c(3,0),
x21=c(1,0), x21=c(3,0), x22=c(1,0), x22=c(3,0)), burnin=20000, mcmc=80000, thin=10, verbose=0,
seed=NA, lambda.start=NA, psi.start=NA, l0=0, L0=0, a0=0.001, b0=0.001, store.scores=FALSE, std.var=TRUE)


summary(whdfan.2)


#Results

Iterations = 20001:99991
Thinning interval = 10 
Number of chains = 1 
Sample size per chain = 8000 

1. Empirical mean and standard deviation for each variable,
   plus standard error of the mean:

               Mean       SD  Naive SE Time-series SE
Lambdax1_3  -0.7445 0.015494 1.732e-04      1.886e-04
Lambdax2_3  -0.8599 0.014494 1.620e-04      1.908e-04
Lambdax3_3  -0.8401 0.014464 1.617e-04      1.845e-04
Lambdax4_3  -0.8942 0.013936 1.558e-04      1.886e-04
Lambdax5_3  -0.6837 0.015781 1.764e-04      1.952e-04
Lambdax6_3  -0.7501 0.015529 1.736e-04      1.938e-04
Lambdax7_1   0.8499 0.014210 1.589e-04      2.186e-04




#Check if stationarity of MC was achieved

heidel.diag(whdfan.2)



#Results


            Stationarity start     p-value
            test         iteration        
Lambdax1_3  passed          1      0.4529 
Lambdax2_3  passed          1      0.4853 
Lambdax3_3  passed          1      0.7601 
Lambdax4_3  passed          1      0.6475 
Lambdax5_3  passed          1      0.3770 
Lambdax6_3  passed          1      0.2505 
Lambdax7_1  passed       2401      0.1764 
Lambdax8_1  passed          1      0.7125 
Lambdax9_1  passed          1      0.0973 
Lambdax10_1 passed        801      0.0875 
Lambdax11_1 failed         NA      0.0290 
Lambdax12_1 passed          1      0.2099 
Lambdax13_1 passed       3201      0.0759 
Lambdax14_1 passed          1      0.0624 
Lambdax15_1 passed          1      0.0727 
Lambdax16_1 passed          1      0.1285 



#See if the densities of the factor loadings are normally distributed

plot(whdfan.2)




#Factor analysis for reduced data, i.e. significant genes
datfana<-read.delim("signGenesBaselineG1.txt", as.is=T, header=TRUE)
dat2fana<-as.matrix(datfana, ncol=22)


#Do a usual factor analysis
fan.1<-factanal(~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20+x21+x22, factors=3, rotation="varimax", data=datfana)
fan.1

#Results

Call:
factanal(x = ~x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 +     x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + x20 +     x21 + x22, factors = 3, data = datfana, rotation = "varimax")

Uniquenesses:
   x1    x2    x3    x4    x5    x6    x7    x8    x9   x10   x11   x12   x13 
0.289 0.220 0.245 0.163 0.415 0.324 0.190 0.367 0.191 0.220 0.214 0.197 0.383 
  x14   x15   x16   x17   x18   x19   x20   x21   x22 
0.479 0.186 0.284 0.314 0.345 0.127 0.230 0.156 0.147 

Loadings:
    Factor1 Factor2 Factor3
x1  0.304   0.768   0.167  
x2  0.352   0.737   0.336  
x3  0.269   0.795   0.225  
x4  0.251   0.857   0.200  
x5  0.319   0.666   0.199  
x6  0.230   0.778   0.134  
x7  0.825   0.303   0.192  
x8  0.713   0.274   0.223  
x9  0.839   0.227   0.232  
x10 0.786   0.326   0.235  
x11 0.803   0.278   0.253  
x12 0.832   0.254   0.216  
x13 0.685   0.293   0.249  
x14 0.572   0.367   0.243  
x15 0.830   0.207   0.286  
x16 0.755   0.167   0.343  
x17 0.787   0.246          
x18 0.738   0.286   0.170  
x19 0.312   0.339   0.813  
x20 0.255   0.154   0.825  
x21 0.318   0.186   0.842  
x22 0.254   0.362   0.811  

               Factor1 Factor2 Factor3
SS loadings      7.901   4.753   3.658
Proportion Var   0.359   0.216   0.166
Cumulative Var   0.359   0.575   0.741

Test of the hypothesis that 3 factors are sufficient.
The chi square statistic is 1244.36 on 168 degrees of freedom.
The p-value is 1.42e-163 


#The factor analysis results confirm the groups given by cluster analysis
#Group1 is factor2, group2 is factor1 and group3 is factor3
#First 6 variables x1-x6 form factor2, vars x7-x18 form factor 1 and x19-x22 form factor 3


#Do a Bayesian factor analysis using MCMC method
 
fan.2<-MCMCfactanal(~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20+x21+x22, factors=3, data=datfana, 
lambda.constraints=list(x1=c(1,0), x1=c(3,0), x2=c(1,0), x2=c(3,0), x3=c(1,0), x3=c(3,0),
x4=c(1,0), x4=c(3,0), x5=c(1,0), x5=c(3,0), x6=c(1,0), x6=c(3,0), x7=c(2,0), x7=c(3,0),
x8=c(2,0), x8=c(3,0), x9=c(2,0), x9=c(3,0), x10=c(2,0), x10=c(3,0), x11=c(2,0), x11=c(3,0), 
x12=c(2,0), x12=c(3,0), x13=c(2,0), x13=c(3,0),  x14=c(2,0), x14=c(3,0), x15=c(2,0), x15=c(3,0), x16=c(2,0), x16=c(3,0),
x17=c(2,0), x17=c(3,0), x18=c(2,0), x18=c(3,0), x19=c(1,0), x19=c(2,0), x20=c(1,0), x20=c(2,0),
x21=c(1,0), x21=c(2,0), x22=c(1,0), x22=c(2,0)), burnin=15000, mcmc=50000, thin=10, verbose=0,
seed=NA, lambda.start=NA, psi.start=NA, l0=0, L0=0, a0=0.001, b0=0.001, store.scores=FALSE, std.var=TRUE)

summary(fan.2)

#Results

Iterations = 15001:64991
Thinning interval = 10 
Number of chains = 1 
Sample size per chain = 5000 

1. Empirical mean and standard deviation for each variable,
   plus standard error of the mean:

               Mean      SD  Naive SE Time-series SE
Lambdax1_2  -0.8510 0.03007 0.0004253      0.0005411
Lambdax2_2  -0.8689 0.02980 0.0004215      0.0005494
Lambdax3_2  -0.8707 0.02969 0.0004198      0.0005505
Lambdax4_2  -0.9140 0.02862 0.0004048      0.0005609
Lambdax5_2  -0.7754 0.03152 0.0004458      0.0005289
Lambdax6_2  -0.8240 0.03050 0.0004314      0.0005430
Lambdax7_1  -0.9051 0.02851 0.0004031      0.0006711
Lambdax8_1  -0.8024 0.03063 0.0004332      0.0006451
Lambdax9_1  -0.9059 0.02852 0.0004034      0.0006719
Lambdax10_1 -0.8867 0.02957 0.0004182      0.0006882
Lambdax11_1 -0.8906 0.02877 0.0004068      0.0007073
Lambdax12_1 -0.9029 0.02865 0.0004052      0.0006788
Lambdax13_1 -0.7917 0.03106 0.0004392      0.0006464
Lambdax14_1 -0.7083 0.03261 0.0004612      0.0006272
Lambdax15_1 -0.9026 0.02832 0.0004004      0.0006694
Lambdax16_1 -0.8351 0.02975 0.0004207      0.0006336
Lambdax17_1 -0.8172 0.03012 0.0004260      0.0006264
Lambdax18_1 -0.8102 0.03055 0.0004320      0.0006477
Lambdax19_3  0.9345 0.02842 0.0004020      0.0005700
Lambdax20_3  0.8692 0.02993 0.0004233      0.0005797
Lambdax21_3  0.9137 0.02877 0.0004069      0.0005939
Lambdax22_3  0.9189 0.02862 0.0004048      0.0005711
Psix1        0.2841 0.01728 0.0002443      0.0002618
Psix2        0.2534 0.01581 0.0002236      0.0002206
Psix3        0.2505 0.01579 0.0002233      0.0002233
Psix4        0.1743 0.01267 0.0001792      0.0001792
Psix5        0.4067 0.02321 0.0003283      0.0003283
Psix6        0.3279 0.01954 0.0002763      0.0002763
Psix7        0.1928 0.01145 0.0001619      0.0001619
Psix8        0.3660 0.02018 0.0002854      0.0002786
Psix9        0.1919 0.01133 0.0001602      0.0001602
Psix10       0.2247 0.01305 0.0001845      0.0001845
Psix11       0.2187 0.01254 0.0001773      0.0001773
Psix12       0.1968 0.01170 0.0001654      0.0001654
Psix13       0.3829 0.02041 0.0002887      0.0002954
Psix14       0.5073 0.02670 0.0003776      0.0003776
Psix15       0.1970 0.01178 0.0001665      0.0001627
Psix16       0.3128 0.01703 0.0002408      0.0002356
Psix17       0.3410 0.01860 0.0002631      0.0002631
Psix18       0.3530 0.01926 0.0002724      0.0002724
Psix19       0.1315 0.01051 0.0001486      0.0001530
Psix20       0.2484 0.01530 0.0002164      0.0002213
Psix21       0.1695 0.01216 0.0001720      0.0001769
Psix22       0.1607 0.01181 0.0001670      0.0001670



#Check if stationarity of MC was achieved

heidel.diag(fan.2)

#Results

            Stationarity start     p-value
            test         iteration        
Lambdax1_2  passed          1      0.6260 
Lambdax2_2  passed          1      0.3754 
Lambdax3_2  passed          1      0.1676 
Lambdax4_2  passed          1      0.4186 
Lambdax5_2  passed          1      0.4105 
Lambdax6_2  passed          1      0.5269 
Lambdax7_1  passed          1      0.5699 
Lambdax8_1  passed          1      0.2551 
Lambdax9_1  passed          1      0.2061 
Lambdax10_1 passed          1      0.0545 
Lambdax11_1 passed          1      0.2233 
Lambdax12_1 passed          1      0.0954 
Lambdax13_1 passed          1      0.8908 
Lambdax14_1 passed          1      0.2277 
Lambdax15_1 passed          1      0.2989 
Lambdax16_1 passed          1      0.5387 
Lambdax17_1 passed          1      0.4305 
Lambdax18_1 passed          1      0.2273 
Lambdax19_3 passed          1      0.4893 
Lambdax20_3 passed          1      0.5351 
Lambdax21_3 passed          1      0.3302 
Lambdax22_3 passed          1      0.3806 
Psix1       passed          1      0.3420 
Psix2       passed          1      0.6259 
Psix3       passed          1      0.0757 
Psix4       passed          1      0.6760 
Psix5       passed       2001      0.1220 
Psix6       passed          1      0.3922 
Psix7       passed          1      0.3686 
Psix8       passed          1      0.2840 
Psix9       failed         NA      0.0494 
Psix10      passed          1      0.2189 
Psix11      passed          1      0.4640 
Psix12      passed          1      0.2293 
Psix13      passed       1501      0.0521 
Psix14      passed        501      0.0891 
Psix15      passed          1      0.1131 
Psix16      passed          1      0.4042 
Psix17      passed          1      0.8547 
Psix18      passed          1      0.5718 
Psix19      passed          1      0.9045 
Psix20      passed          1      0.8294 
Psix21      passed          1      0.8385 
Psix22      passed          1      0.7041 
                                      
            Halfwidth Mean   Halfwidth
            test                      
Lambdax1_2  passed    -0.851 0.001061 
Lambdax2_2  passed    -0.869 0.001077 
Lambdax3_2  passed    -0.871 0.001079 
Lambdax4_2  passed    -0.914 0.001099 
Lambdax5_2  passed    -0.775 0.001037 
Lambdax6_2  passed    -0.824 0.001064 
Lambdax7_1  passed    -0.905 0.001315 
Lambdax8_1  passed    -0.802 0.001264 
Lambdax9_1  passed    -0.906 0.001317 
Lambdax10_1 passed    -0.887 0.001349 
Lambdax11_1 passed    -0.891 0.001386 
Lambdax12_1 passed    -0.903 0.001331 
Lambdax13_1 passed    -0.792 0.001267 
Lambdax14_1 passed    -0.708 0.001229 
Lambdax15_1 passed    -0.903 0.001312 
Lambdax16_1 passed    -0.835 0.001242 
Lambdax17_1 passed    -0.817 0.001228 
Lambdax18_1 passed    -0.810 0.001270 
Lambdax19_3 passed     0.934 0.001117 
Lambdax20_3 passed     0.869 0.001136 
Lambdax21_3 passed     0.914 0.001164 
Lambdax22_3 passed     0.919 0.001119 
Psix1       passed     0.284 0.000513 
Psix2       passed     0.253 0.000432 
Psix3       passed     0.251 0.000438 
Psix4       passed     0.174 0.000351 
Psix5       passed     0.407 0.000849 
Psix6       passed     0.328 0.000542 
Psix7       passed     0.193 0.000317 
Psix8       passed     0.366 0.000546 
Psix9       <NA>          NA       NA 
Psix10      passed     0.225 0.000362 
Psix11      passed     0.219 0.000348 
Psix12      passed     0.197 0.000324 
Psix13      passed     0.383 0.000699 
Psix14      passed     0.507 0.000777 
Psix15      passed     0.197 0.000319 
Psix16      passed     0.313 0.000462 
Psix17      passed     0.341 0.000516 
Psix18      passed     0.353 0.000534 
Psix19      passed     0.132 0.000300 
Psix20      passed     0.248 0.000434 
Psix21      passed     0.169 0.000347 
Psix22      passed     0.161 0.000327


#looks like there is a problem with psix9, the residual variance of 
# x9, but the trace and density of it shows that the chain 
#converged and variance is normally distributed. The p-value is 0.0494
#so, it's okay.

#See if the densities of the factor loadings are normally distributed

plot(fan.2)

#See in the word file some plots, everything looks okay about the model

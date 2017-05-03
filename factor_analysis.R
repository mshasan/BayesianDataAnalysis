#install the package MCMCpack to do a bayesian factor analysis using MCMC methods

library(MCMCpack)


#Read in the data and transform it into a matrix

datfana<-read.delim("signGenesBaselineG1.txt", as.is=T, header=TRUE)
dat2fana<-as.matrix(datfana, ncol=22)

#genevar<-t(dat2fana)

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

2. Quantiles for each variable:

               2.5%     25%     50%     75%   97.5%
Lambdax1_2  -0.9118 -0.8711 -0.8500 -0.8299 -0.7942
Lambdax2_2  -0.9289 -0.8888 -0.8681 -0.8480 -0.8131
Lambdax3_2  -0.9300 -0.8904 -0.8703 -0.8505 -0.8139
Lambdax4_2  -0.9718 -0.9327 -0.9129 -0.8949 -0.8586
Lambdax5_2  -0.8366 -0.7963 -0.7758 -0.7542 -0.7136
Lambdax6_2  -0.8867 -0.8445 -0.8235 -0.8028 -0.7663
Lambdax7_1  -0.9621 -0.9239 -0.9043 -0.8854 -0.8519
Lambdax8_1  -0.8644 -0.8225 -0.8016 -0.7821 -0.7440
Lambdax9_1  -0.9624 -0.9248 -0.9056 -0.8867 -0.8519
Lambdax10_1 -0.9458 -0.9064 -0.8859 -0.8661 -0.8304
Lambdax11_1 -0.9481 -0.9096 -0.8899 -0.8713 -0.8353
Lambdax12_1 -0.9593 -0.9220 -0.9025 -0.8834 -0.8470
Lambdax13_1 -0.8546 -0.8121 -0.7912 -0.7706 -0.7331
Lambdax14_1 -0.7721 -0.7300 -0.7081 -0.6858 -0.6457
Lambdax15_1 -0.9577 -0.9220 -0.9024 -0.8832 -0.8483
Lambdax16_1 -0.8935 -0.8548 -0.8343 -0.8143 -0.7792
Lambdax17_1 -0.8791 -0.8374 -0.8162 -0.7969 -0.7597
Lambdax18_1 -0.8728 -0.8305 -0.8096 -0.7895 -0.7513
Lambdax19_3  0.8796  0.9151  0.9343  0.9533  0.9919
Lambdax20_3  0.8112  0.8486  0.8688  0.8893  0.9291
Lambdax21_3  0.8603  0.8937  0.9129  0.9331  0.9711
Lambdax22_3  0.8641  0.8990  0.9189  0.9375  0.9759
Psix1        0.2508  0.2724  0.2840  0.2952  0.3185
Psix2        0.2239  0.2422  0.2527  0.2636  0.2862
Psix3        0.2211  0.2397  0.2499  0.2607  0.2830
Psix4        0.1506  0.1654  0.1738  0.1825  0.2001
Psix5        0.3626  0.3910  0.4058  0.4220  0.4534
Psix6        0.2918  0.3143  0.3270  0.3405  0.3690
Psix7        0.1714  0.1851  0.1925  0.2002  0.2166
Psix8        0.3293  0.3523  0.3655  0.3786  0.4077
Psix9        0.1708  0.1841  0.1918  0.1991  0.2155
Psix10       0.2001  0.2160  0.2243  0.2335  0.2512
Psix11       0.1956  0.2101  0.2183  0.2269  0.2445
Psix12       0.1755  0.1888  0.1963  0.2043  0.2211
Psix13       0.3455  0.3688  0.3821  0.3964  0.4242
Psix14       0.4563  0.4889  0.5070  0.5252  0.5602
Psix15       0.1748  0.1887  0.1966  0.2046  0.2211
Psix16       0.2806  0.3012  0.3124  0.3241  0.3472
Psix17       0.3068  0.3281  0.3401  0.3531  0.3793
Psix18       0.3172  0.3394  0.3522  0.3659  0.3922
Psix19       0.1121  0.1243  0.1311  0.1385  0.1528
Psix20       0.2191  0.2377  0.2480  0.2586  0.2794
Psix21       0.1467  0.1610  0.1690  0.1778  0.1945
Psix22       0.1377  0.1526  0.1603  0.1684  0.1850

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


#looks like there is a problem with psix9, the variance of the factor loading 
#corresponding to x9, but the trace and density of it shows that the chain 
#converged and variance is normally distributed

#See if the densities of the factor loadings are normally distributed

plot(fan.2)

#See in the word file some plots, everything looks okay about the model

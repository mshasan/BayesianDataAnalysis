## Problem 2 ----------------

library(ggplot2)
library(grid)
library(gridExtra)

x <- seq(0, 1, length = 100)
a <- dbeta(x, .5, .5)
b <- dbeta(x, 10.2, 1.5)
c <- dbeta(x, 1.5, 10.2)
d <- dbeta(x, 100, 62)


p1 <- qplot(x,a, main="Beta (0.5, 0.5)",geom="line")
p2 <- qplot(x,b, main="Beta (10.2, 1.5)",geom="line") 
p3 <- qplot(x,c, main="Beta (1.5, 10.2)",geom="line")
p4 <- qplot(x,d, main="Beta (100, 62)",geom="line") 
grid.arrange(p1, p2, p3, p4, nrow=2, ncol = 2, main = "Different Beta Distribution")



## How to create direct pdf in R

pdf("Different Beta distribution.pdf", width = 8, height = 6)
grid.arrange(p1, p2, p3, p4, nrow=2, ncol = 2, main = "Different Beta Distribution")



## Problem 5--------------------

y <- c(46,58,40,47,47,54,51,50,52,50,53,43,48,50,55,49,50,52,56,49);y
sum((y-51)^2)







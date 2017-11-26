#### Crude comparisons Poisson versus binomial 

xs= seq(0,10,1)
ys=dbinom(x = xs, size = 1000,p=.001)
y2s=dpois(x = xs,lambda = 1)
sum(ys)
sum(y2s)
plot(xs,ys, main ="Comparing Poisson vs Binomial ", type="h", ylab="P(x)", xlab="x", col="red", lwd=3)
points(xs,y2s, type="p", col="blue", pch=20, cex=1)
(ys-y2s)/ys #relative differences in probabilities

xs= seq(0,10,1)
ys=dbinom(x = xs, size = 100,p=.01)
y2s=dpois(x = xs,lambda = 1)
sum(ys)
sum(y2s)
plot(xs,ys, main ="Comparing Poisson vs Binomial ", type="h", ylab="P(x)", xlab="x", col="red", lwd=3)
points(xs,y2s, type="p", col="blue", pch=20, cex=1)
(ys-y2s)/ys

xs= seq(0,10,1)
ys=dbinom(x = xs, size = 20,p=.05)
y2s=dpois(x = xs,lambda = 1)
sum(ys)
sum(y2s)
plot(xs,ys, main ="Comparing Poisson vs Binomial ", type="h", ylab="P(x)", xlab="x", col="red", lwd=3)
points(xs,y2s, type="p", col="blue", pch=20, cex=1)
(ys-y2s)/ys


xs= seq(0,10,1)
ys=dbinom(x = xs, size = 2,p=.5)
y2s=dpois(x = xs,lambda = 1)
sum(ys)
sum(y2s)
plot(xs,ys, main ="Comparing Poisson vs Binomial ", type="h", ylab="P(x)", xlab="x", col="red", lwd=3)
points(xs,y2s, type="p", col="blue", pch=20, cex=1)

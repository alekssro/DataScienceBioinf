#Lecture of Data Science | Week02 | 5/9/17
#Checking how the sample size influences the mean and the sd
d1 <- rnorm(10, 300, 500)
cat("mean:", mean(d1), "standard deviation:", sd(d1))

d2 <- rnorm(100000, 300, 500)
cat("mean:", mean(d2), "standard deviation:", sd(d2))

#sample size only affect the precission of the standard deviation/mean

#Find 3 positive values that will give you the result above (mean = 300, sd = 500)
p <- c(5, 150, 755)
cat("mean:", mean(p), "standard deviation:", sd(p))

#Exercise
x = c(70:90, 1000, 1100, 1200)
mean(x)
median(x)
sd(x)
quantile(x)
mad(x)

# coefficients of variation: we use this to compare different ratio scales 
# (for example, when we measure the sd of the weight of a mouse and the same with an elephant,
#   we won't be able to compare them unless we take the mean into account')

c <- c(0, 10, 20, 30, 40)
f <- c(32, 50, 68, 86, 104)
sd(c)/mean(c)*100
sd(f)/mean(f)*100

kelvin = c + 273
rankine = f + 459.67
sd(kelvin)/mean(kelvin)*100
sd(rankine)/mean(rankine)*100

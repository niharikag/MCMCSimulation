initValue = rep(1,3)
m1 = rep(0,3)
sigma1.inv = matrix(c(1,.6,0,.6,1,.4,0,.4,1),3,3)
x1 = mvnLangevin(20000, sigma1.inv, initValue)
z1 = mvrnorm(10000,m1,solve(sigma1.inv))
par(mfrow=c(1,2))
x = x1[2000:20000,]
plot(x)
plot(z1)
#compare mean of first component
mean(x[,1:1])
mean(z1[,1:1])

initValue = c(0,0)
m2 = c(0,0)
sigma2.inv = matrix(c(1,.8,.8,1),2,2)
x2 = mvnLangevin(12000, sigma2.inv, initValue)
z2 = mvrnorm(10000,m2, solve(sigma2.inv))
x = x2[2000:12000,]
par(mfrow=c(1,2))
plot(x)
plot(z2)
#compare mean
mean(x[,1:1])
mean(z2[,1:1])

#Example 3: bivariate with non centralized data
initValue = c(1,1)
m3 = c(-1,3)
sigma3.inv = matrix(c(1,.8,.8,1),2,2)
x3 = mvnLangevin(12000,sigma3.inv,initValue,m3)
z3 = mvrnorm(10000,m3, solve(sigma3.inv))
x = x3[2000:12000,]
par(mfrow=c(1,2))
plot(x)
plot(z3)
#compare mean
mean(x[,1:1])
mean(z3[,1:1])

#now, call mcmc function for 50 variables with zero mean
hd = 1000 #high dimension
initValue = rnorm(hd)
m2 = rep(0,hd)
sigma2.inv = diag(hd)
x2 = mvnMCMC(50000, sigma2.inv, initValue)
z2 = mvrnorm(10000,m2, solve(sigma2.inv))
x = x2[20000:50000,]
par(mfrow=c(1,2))
plot(x)
plot(z2)
mean(x[,1:1])
mean(z2[,1:1])
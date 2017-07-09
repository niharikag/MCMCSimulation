#to be deleted
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


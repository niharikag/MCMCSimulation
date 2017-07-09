
dimension = 7
iteration = 50000
initValue = rep(0,dimension) #initValue = runif(dim4)
meanValue = rep(0,dimension)
#generate graph
g1 = generate.randomGraph(dimension,.4) #plot.igraph(g1)
precisionMatrix = graph2cov(g1,dimension)
jumpFactor = .1

#Common code for Testing  
#Call simulation method and exact sampling method
data1 = mvnMCMC(iteration, precisionMatrix, initValue, jumpFactor, meanValue)
data2= mvnMCMC(iteration, precisionMatrix, initValue, jumpFactor, meanValue)
z1 = mvrnorm(1000,meanValue,solve(precisionMatrix))
#trace plot and density plot
chain1 = mcmc(data1)
#plot(chain1)
#acfplot(chain1)
chain2 = mcmc(data2)
chain = mcmc.list(chain1,chain2)
gelman.plot(chain)
# discard first half of the iterations
index = seq((iteration/2),iteration,10)
x = data1[index,]
plot(x)
plot(z1)
#compare mean of first component
mean(x[,1:1])
mean(z1[,1:1])

#common code for graph plot
cx = cov(x)
gls = glasso(cx,.003)

par(mfrow=c(1,2))
plot.igraph(g1)
g2 = cov2graph(gls$wi)
plot.igraph(g2)

round(gls$wi,1)
round(precisionMatrix,1)


h1 = huge(cx,.005, method="glasso")
par(mfrow=c(1,2))
plot.igraph(g1)
plot(h1)

h2 = huge(cx)
plot(h2)


#final sequesnce
dimension = 10
iteration = 100000
initValue = rep(0,dimension) #initValue = runif(dim4)
meanValue = rep(0,dimension)
#generate graph
g1 = generate.randomGraph(dimension,.3) #plot.igraph(g1)
precisionMatrix = graph2cov(g1,dimension)
jumpFactor = .15

data1 = mvnMCMC(iteration, precisionMatrix, initValue, jumpFactor, meanValue)
index = seq((iteration/2),iteration,10)
x = data1[index,]
cx = cov(x)
h1 = huge(cx,.005, method="glasso")
par(mfrow=c(1,2))
plot.igraph(g1)
plot(h1)
h2 = huge(cx)
plot(h2)
round(solve(cx),2)
floor(precisionMatrix)
z1 = mvrnorm(50000,meanValue,solve(precisionMatrix))
c1 = cov(z1)
round(solve(c1),2)


#glasso#
cx = cov(x)
gls = glasso(cx,.003)
par(mfrow=c(1,2))
plot.igraph(g1)
round(gls$wi,2)
g2 = cov2graph(gls$wi)
plot.igraph(g2)

round(gls$wi,2)

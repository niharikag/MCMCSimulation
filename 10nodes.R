#final sequesnce
dimension = 10
iteration = 100000
initValue = rep(0,dimension) #initValue = runif(dim4)
meanValue = rep(0,dimension)
#generate graph
hg = huge.generator(100,dimension,"hub",g=3)
g1 = getAdjMatrix(hg$theta)
#g1 = generate.randomGraph(dimension,.3) #plot.igraph(g1)
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
gls = glasso(cx,.0055)
par(mfrow=c(1,2))
plot.igraph(g1,sub="A Randomly generated graph",frame=TRUE)
g2 = cov2graph(gls$wi)
plot.igraph(g2,sub="A Graph obtained through GLASSO",frame=TRUE)

data2= mvnMCMC(iteration, precisionMatrix, initValue, jumpFactor, meanValue)

chain1 = mcmc(data1)

chain2 = mcmc(data2)
chain = mcmc.list(chain1,chain2)
gelman.plot(chain)
acfplot(chain1)

#Call Langevin Algorithm
Ldata1 = mvnLangevin(iteration, precisionMatrix, initValue, meanValue, jumpFactor)
Ldata2= mvnLangevin(iteration, precisionMatrix, initValue, meanValue, jumpFactor)
Lchain1 = mcmc(Ldata1)
Lchain2 = mcmc(Ldata2)
Lchain = mcmc.list(Lchain1,Lchain2)
gelman.plot(Lchain)
acfplot(Lchain1)

index = seq((iteration/2),iteration,10)
x = Ldata1[index,]

#glasso#
cx = cov(x)
gls = glasso(cx,.0009)
par(mfrow=c(1,2))
plot.igraph(g1,sub="A Randomly generated graph",frame=TRUE)
g2 = cov2graph(gls$wi)
plot.igraph(g2,sub="A Graph obtained through GLASSO",frame=TRUE)

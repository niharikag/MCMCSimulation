#testing mvnMH
#TEST:1, call mcmc function for 3 variables with zero mean
#initializations
dimension = 3
iteration = 10000
initValue = rep(0,dimension)
meanValue = rep(0, dimension)
precisionMatrix = matrix(c(1,.6,0,.6,1,.4,0,.4,1),3,3)
jumpFactor = .9

#Common code for Testing  
#Call simulation method and exact sampling method
data1 = mvnMCMC(iteration, precisionMatrix, initValue, jumpFactor, meanValue)
data2= mvnMCMC(iteration, precisionMatrix, initValue, jumpFactor, meanValue)
z1 = mvrnorm(1000,meanValue,solve(precisionMatrix))
#trace plot and density plot
chain1 = mcmc(data1)
plot(chain1)
#acfplot(chain1)
chain2 = mcmc(data2)
chain = mcmc.list(chain1,chain2)
gelman.plot(chain)
# discard first half of the iterations
par(mfrow=c(1,2))

index = seq((iteration/2):iteration,10)
x = data1[index,]
plot(x)
plot(z1)
#compare mean of first component
mean(x[,1:1])
mean(z1[,1:1])

#common code for graph plot
#cx = cov(z1)
cx = cov(x)
gls = glasso(cx,.007)
round(gls$wi,2)
round(precisionMatrix)
par(mfrow=c(1,2))
plot.igraph(g1)
g2 = cov2graph(gls$wi)
plot.igraph(g2)

#TEST 2: bivariate normal
#initializations
dimension = 2
iteration = 10000
initValue = rep(0,dimension)
meanValue = rep(0, dimension)
precisionMatrix = matrix(c(1,.8,.8,1),2,2)
jumpFactor = .9


#TEST 3:  bivariate non centralized 
#initializations
dimension = 2
iteration = 10000
initValue = rep(0,dimension)
meanValue = c(-1,3)
precisionMatrix = matrix(c(1,.8,.8,1),2,2)
jumpFactor = 1


#TEST 3: call mcmc function for 100 variables with zero mean
#initializations
dimension = 100
iteration = 50000
initValue = rep(0,dimension) #initValue = runif(dim4)
meanValue = rep(0,dimension)
precisionMatrix = diag(dimension)
jumpFactor = 2.38/sqrt(dimension)




#TEST 4: simulate data for randomly generated graph
#initializations
dimension = 7
iteration = 50000
initValue = rep(0,dimension) #initValue = runif(dim4)
meanValue = rep(0,dimension)
#generate graph
g1 = generate.randomGraph(dimension,.4) #plot.igraph(g1)
precisionMatrix = graph2cov(g1,dimension)
jumpFactor = .1


#common code for graph plot
cx = cov(x)
gls = glasso(cx,.01)
round(gls$wi,2)
round(precisionMatrix)
par(mfrow=c(1,2))
plot.igraph(g1)
g2 = cov2graph(gls$wi)
plot.igraph(g2)
#test 2, 
#generate random graph, inverse cov matrix
#generate data from inv. cov matrix (burn-in 2000, thining 10)
#and run glasso, to check if glasso returns almost same inverse cov matrix
dm = 7
g1 = generate.sparseGraph(dm,.4)
iCv = g1$iCovMat
#iCr = cov2cor(iCv)

initMean =rep(0,dm)
initSigma = .01*diag(dm)

x1 = mvnMCMC(50000, iCv, initMean, .1)

s = seq(from=10000,to=50000,by = 10)
x = x1[s,]


cx = cov(x)
#cr = cor(cx)
gls = glasso(cx,.01)
round(gls$wi,2)
round(precisionMatrix)
iCv
gls$wi
round(gls$wi,2)

h1 = huge.glasso(cr,.1)
h1$icov
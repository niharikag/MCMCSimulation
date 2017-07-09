# simulate data samples from multivariate normal distribution 
# given concentration/matrix

#MCMC algorithm
#step1: write a target function from where we want to have data samples

#declare global variables
gPrecisionMatrix = NULL
gMean  = NULL #target Mean
gSD  = NULL #standard deviation for proposal
gDim  = NULL #dimension


targetFunc= function(x){
  #constant Terms will get canceled
  #print(x)
  y = exp(-0.5*(t(x-gMean)%*% gPrecisionMatrix%*% (x-gMean)))
  return(y)
}


#generate sample from proposalproposal density
#it takes currValue and dimension as input
proposalFunc= function(currValue){
  
  y = currValue+gSD*rnorm(gDim)
  #y = currValue+rnorm(dm,0,2.38/sqrt(dm))
  return(y)
}

#step2: write Metropolis-Hastings algorithm, to sample 
#from a distribution that is proportional to the target function
#sc: sigma
mvnMCMC = function(numIter, precisionMatrix, startValue=NULL, sd =NULL, m=NULL){
  
  assign("gDim", nrow(precisionMatrix), envir = .GlobalEnv)
  
  #gDim  <<- nrow(precisionMatrix)
  if(is.null(m)){
    m = rep(0,gDim)
  }
  
  assign("gPrecisionMatrix", precisionMatrix, envir = .GlobalEnv)
  assign("gSD", sd, envir = .GlobalEnv)
  assign("gMean", m, envir = .GlobalEnv)
  #print(gMean)
  #print(gPrecisionMatrix)
  
  
  if(is.null(gSD)){
    #gSD = 2.38/sqrt(gDim)
    assign("gSD", 2.38/sqrt(gDim), envir = .GlobalEnv)
  }
  
  
  #initialize the data matrix with zero
  x = matrix(rep(0,gDim*numIter), numIter, gDim)
  
  #initialize the first observation based startValue
  if(is.null(startValue)){
    startValue = proposalFunc(m)
  }
  
  
  
  x[1,] = startValue
  nAccepted = 0;
  
  for(i in 2:numIter){
    x.curr = x[i-1,]
    #x.proposed = mvrnorm(1,x.curr,initSigma)
    x.proposed = proposalFunc(x.curr)
    mcmcRatio  =  targetFunc(x.proposed)/ targetFunc(x.curr)
    u = runif(1)
    if(is.na(mcmcRatio)){
      x[i,] = x.curr
    }
    else{
      if( u< mcmcRatio){
        #accept the proposed value with probability min(1,mcmcRatio)
        nAccepted = nAccepted+1;
        x[i,] = x.proposed
      }
      else{
        x[i,] = x.curr
      }
    }
  }
  cat("Acceptance Rate:", nAccepted/numIter*100 ,"%")
  return(x)
}



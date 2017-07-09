#MCMC Langevin algorithm
#step1: write a target function from where we want to have data samples

#declare global variables
gPrecisionMatrix = NULL
gMean  = NULL #target Mean
gDS  = NULL #discretization size
gDim  = NULL #dimension


tarLangFunc= function(x){
  #constant Terms are not needed
  y = exp(-0.5*(t(x-gMean)%*% gPrecisionMatrix%*% (x-gMean)))
  return(y)
}

gradLangFunc= function(x){
  #constant Terms will get canceled
  y = -t(x-gMean)%*% gPrecisionMatrix
  return(as.vector(t(y)))
}

#generate sample from proposalproposal density
#it takes currValue and dimension as input
propLangFunc= function(currValue){
  y = currValue+ (.5*(gDS)*gradLangFunc(currValue))+gDS*rnorm(gDim,0,1)
  return(y)
}

mcmcLangRatio = function(x,y){
  term1 = tarLangFunc(y)/ tarLangFunc(x)
  exp1 = (y-x - (.5*(gDS^2)*gradLangFunc(x)))
  exp2 = (x-y-(.5*(gDS^2)*gradLangFunc(y)))
  term2 = exp((sum(exp2^2) - sum(exp1^2))*.5*(1/gDS^2))
  ratio = min(1,term1*term2)
  return(ratio)
}
#step2: write Metropolis-Hastings algorithm, to sample 
#from a distribution that is proportional to the target function
mvnLangevin = function(numIter, precisionMatrix, initValue=NULL, m=NULL,ds = 1){
  # to be changed
  
  assign("gDim", nrow(precisionMatrix), envir = .GlobalEnv)
  
  if(is.null(m)){
    m = rep(0,dimension)
  }
  
  #initialize the data matrix with zero
  x = matrix(rep(0,dimension*numIter), numIter,dimension)
  
  #initialize the first observation based initValue
  if(is.null(initValue)){
    initValue = m
  }
  
  assign("gPrecisionMatrix", precisionMatrix, envir = .GlobalEnv)
  assign("gDS", ds, envir = .GlobalEnv)
  assign("gMean", m, envir = .GlobalEnv)
  
  
  x[1,] = initValue
  nAccepted = 0;
  
  for(i in 2:numIter){
    x.curr = x[i-1,]
    x.proposed = propLangFunc(x.curr)
    mcmcRatio  =  mcmcLangRatio(x.curr,x.proposed)
    u = runif(1)
    
    if( u< mcmcRatio){
      #accept the proposed value with probability min(1,mcmcRatio)
      nAccepted = nAccepted+1;
      x[i,] = x.proposed
    }
    else{
      x[i,] = x.curr
    }
  }
  cat("Acceptance Rate:", nAccepted/numIter*100 ,"%")
  return(x)
}
library(MASS)
library(Matrix)
library(glasso)
library(igraph)
library(coda)
library(huge)
#this method generates a symmetric positive definite matrix of size n
generate.psdMat = function(n){
  #generate a random n by n matrix
  m = matrix(runif(n*n),n,n)
  
  #make is symmetric
  m = m + t(m)
  
  #make it positive definite by adding n*I
  m = m+n*diag(n)
  
  return(n)
}

#this method generates a sparse symmetric positive definite matrix of size n

generate.sparse.spdMat = function(n){
  #generate a sparse random n by n matrix
  m = matrix(0, n, n, sparse=TRUE)
  
  #make is symmetric
  m = m + t(m)
  
  #make it positive definite by adding n*I
  m = m+n*diag(n)
  
  return(n)
}

#this method generates a sparse graph of n nodes 
# and translates its adjacency matrix into inverse covariance matrix
library(igraph)
generate.sparseGraph = function(n,d){
  #generate a random graph
  g = sample_gnp(n,d)
  adj = as_adjacency_matrix(g)
  
  #E(g)$weight <- runif(ecount(g))
  #E(g)$weight <- rnorm(ecount(g))
  #E(g)$weight <- runif(ecount(g),-1,1)
  #E(g)$weight <- runif(ecount(g),0.5,0.8) #to be changed
  E(g)$weight <- runif(ecount(g),0.4,0.5) #to be changed
  
  m1 = as_adjacency_matrix(g, attr="weight")
  
  m = apply(m1, 1:2, function(x) { if(is.na(x)){x=0}else{x}})
  
  #convert this adjacency matrix into a SPD inverse covariance matrix
  #make is symmetric
  #m = m + t(m)
  
  #make it positive definite by adding n*I
  m = m+n*diag(n)
  
  ls = list(adjMat = adj, iCovMat = m)
  return(ls)
}
#this method generates a sparse graph of n nodes 
# and returns its adjacency matrix
library(igraph)
generate.randomGraph = function(n,d){
  #generate a random graph
  g = sample_gnp(n,d)
  
  #m = as_adjacency_matrix(g)
  
  return(g)
}

#this method translates a graph object into a covariance matrix
library(igraph)
graph2cov = function(g,n){
  #put weight for each edge
  E(g)$weight <- round(runif(ecount(g),-0.8,-0.5),2) #to be changes
  
  m1 = as_adjacency_matrix(g, attr="weight")
  
  m = apply(m1, 1:2, function(x) { if(is.na(x)){x=0}else{x}})
  
  #convert this adjacency matrix into a SPD inverse covariance matrix
  #make is symmetric
  
  #make it positive definite by adding n*I
  m = m+n*diag(n)
  
  return(m)
}
#convert cov matrix to graph
cov2graph = function(m){
  m = apply(m, 1:2, function(x) { if(is.null(x) || abs(x)< 0.05){x=0}else{1}})
  m = m +-1*diag(nrow(m))
  m
  g = graph.adjacency(m,mode="undirected")
  return(g)
}



getAdjMatrix = function(m){
  m = as.matrix(m)
  g = graph.adjacency(m,mode="undirected")
  return(g)
}


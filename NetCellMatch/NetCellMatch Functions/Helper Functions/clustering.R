#' A clustering function
#'
#'This function creates the clustering results from the pipeline 
#'Run this function after the graph wavelets
#'@param pd number of nodes in Laplacian
#'@param scales number of scales to consider 
#'@param Wave is the wavelet parameter (run from graph.wavelet) 
#'
#'

clustering=function(pd,scales,Wave){
pd=pd
Dist = vector("list",scales) #initialize the Distance matrices
T <- 100
set.seed(8773)
R <- matrix(rnorm(pd*T,mean=0,sd=1), pd, T)
#Get approximate correlation distance between wavelets

FWT <- vector("list",scales)
for(i in 1:scales){
  #print(i)
  FWT[[i]] <- t(Wave[[i]])%*%R
  Dist[[i]] <- 1 - cor(t(FWT[[i]]))
  cat("Distance matrix at Scale:",i," ")
  
}

result.cluster = hclust.Dist(Dist = Dist,scales)

clustmem <- list()
nclust <- rep(0,scales)
#Do it by number of scales
for(ii in 1:scales){
  h <- min(cut_height(ii,pd)$min_val_cut) - 0.01
  nclust[ii] <- max(cutree(result.cluster[[ii]], h = h))
  clustmem[[ii]] = cutree(result.cluster[[ii]],h = h)
  print(clustmem[[ii]])
  cat("Scale Completed: ", ii)
}
return(clustmem)
}




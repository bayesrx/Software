#'Graph Wavelet Function
#'
#'This function calculates graphical wavelets at different scales of the laplacian
#'@param dimension is number of nodes in graph
#'@param Lp is Laplacian
#'@param nscales is number of scales 
#'

Graph.wavelet <- function(dimension,Lp,nscales){
  p=dimension
  ev <- eigen(Lp)
  ev.val <- ev$value[p:1] #sort from lowest to highest
  ev.vec <- ev$vector[,p:1] 
  
  #filter matrix
  Gs <- list()
  x1 <- 1
  x2 <- x1/ev.val[2]
  alpha <- 2
  beta <- 1/log10(ev.val[3]/ev.val[2])  #Potential Problem Here  
  #(For manhattan- second and 3rd eigenvalues identical)
  s.min = x1/ev.val[2]
  s.max = x1/ev.val[2]^2
  library(sfsmisc)
  s = lseq(from = s.min, to = s.max, length=nscales)   #very miniscule scale determined by the eigenvalues 
  for(j in 1:length(s)){
    Gs[[j]] <- matrix(0,p,p)
    #print(Gs[[j]])
    for(i in 1:p){
      Gs[[j]][i,i] = gfunc(s[j]*ev.val[i],alpha,beta,x1,x2)
    }
  }
  
  Wave <- list()
  for(j in 1:length(s)){
    Wave[[j]] <- ev.vec %*% Gs[[j]] %*% t(ev.vec)
    print(Wave[[j]])
    cat("Scale Completed: ", j)
  }
  
  #define the distance matrix at each scale
  res <- list()
  res <- NULL
  res$Wave <- Wave
  res$s <- s
  res$eval <- ev.val
  return(res)
}
